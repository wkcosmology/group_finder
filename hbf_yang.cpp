/**
 * File            : hbf_yang.cpp
 * Author          : Kai Wang <wkcosmology@gmail.com>
 * Date            : 10.18.2022
 * Last Modified By: Kai Wang <wkcosmology@gmail.com>
 *
 * Implementation of Yang's halo-based group finder.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include <hippcntl.h>
#include <hippio.h>

#include "nbr_finder.hpp"

using namespace std;

using HIPP::IO::H5File;

class LinearInterp
{
   public:
    LinearInterp() {}
    LinearInterp(const vector<double> &xp, const vector<double> &yp, bool sorted = false)
    {
        Init(xp, yp);
        if (!sorted)
            SortFp();
    }
    void Init(const vector<double> &xp, const vector<double> &yp, bool sorted = false)
    {
        if (xp.size() != yp.size()) {
            std::cout << "xp and yp have different size" << std::endl;
            exit(EXIT_FAILURE);
        }
        fp_.resize(xp.size());
        for (size_t i = 0; i < xp.size(); ++i)
            fp_[i] = std::make_pair(xp[i], yp[i]);
        if (!sorted)
            SortFp();
    };
    double Interp(double xp) const
    {
        if (fp_.size() == 0) {
            std::cout << "Linear Interp input illegal" << std::endl;
            exit(EXIT_FAILURE);
        }

        if ((xp < fp_.front().first) || std::isnan(xp))
            return fp_.front().second;
        else if (xp >= fp_.back().first)
            return fp_.back().second;

        size_t i = 0;
        while (xp >= fp_[i + 1].first)
            ++i;
        double x1(fp_[i].first), y1(fp_[i].second), x2(fp_[i + 1].first),
            y2(fp_[i + 1].second);
        return y1 + (xp - x1) / (x2 - x1) * (y2 - y1);
    }

   private:
    void SortFp()
    {
        std::sort(
            fp_.begin(),
            fp_.end(),
            [](const std::pair<double, double> &x, const std::pair<double, double> &y) {
                return x.first < y.first;
            });
    }
    vector<std::pair<double, double>> fp_;
};

class HaloMassFunc
{
   public:
    HaloMassFunc(const string &filename, double vol) : kVol_(vol)
    {
        auto f = H5File(filename, "r");
        f.open_dataset("HaloMass").read(hms_);
        f.open_dataset("CHMF").read(chmf_);
    }

    vector<double> GenerateHMSamp(int n_samp)
    {
        vector<double> hm_samp(n_samp);

        for (auto &v: chmf_)
            v *= kVol_;

        if (*chmf_.end() < n_samp)
            cout << "Maximum # of sample: " << *chmf_.end() << ", required #: " << n_samp
                 << endl;

        auto f = LinearInterp(chmf_, hms_);
        for (int i = 0; i < n_samp; ++i)
            hm_samp[i] = f.Interp(i + 1);

        return hm_samp;
    }

   private:
    vector<double> hms_, chmf_;
    const double kVol_;
};

class YangFinder
{
   private:
    struct Group {
        double z, tot_sm, hm;
        double r_vir, sig_v;  // for neighbor searching
        set<int> mems;
        array<double, 3> pos;
    };

   public:
    /**
     * @brief Init the halo-based group finder
     *
     * @param pos [positions, unit (Mpc, Mpc, redshift), z-direction is the
     * third one]
     * @param z [redshift for group searching]
     * @param sm [log of stellar mass]
     * @param hmf [halo mass array, generated from halo mass function]
     * @param period [box size for periodic box]
     */
    YangFinder(
        const vector<array<double, 3>> &pos,
        const vector<double> &z,
        const vector<double> &sm,
        const vector<double> &hmf,
        double period,
        double p_bg)
        : gal_pos_(pos), gal_sm_(sm), hmf_(hmf), gal_z_(z), period_(period), p_bg_(p_bg)
    {
        // input validity check
        n_gal_ = pos.size();
        cout << "# of input galaxies: " << n_gal_ << endl;
        if (gal_sm_.size() != n_gal_)
            cout << "# of stellar mass is wrong: " << gal_sm_.size() << endl;
        if (hmf_.size() != n_gal_)
            cout << "# of halo shape is wrong: " << hmf_.size() << endl;
        if (gal_z_.size() != n_gal_)
            cout << "# of redshift is wrong: " << gal_z_.size() << endl;
        host_hm_.assign(n_gal_, -1);
        if_mm_.assign(n_gal_, -1);
        gal_prob_.assign(n_gal_, -1);
        gal_igrp_.assign(n_gal_, -1);
        // check position
        if (period_ > 0) {
            for (auto &p: gal_pos_) {
                for (auto &x: p) {
                    if (x < 0 || x > period_) {
                        cout << "Position out of range" << x << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    void Run(int n_iter)
    {
        vector<Group> grps;
        // initial the group catalog
        for (int i = 0; i < n_gal_; ++i) {
            Group grp_tmp;
            grp_tmp.mems.insert(i);
            grp_tmp.z = gal_z_[i];
            grp_tmp.pos = gal_pos_[i];
        }
        gal_igrp_ = RegulateGrp(grps);

        // itaration
        for (int it = 0; it < n_iter; ++it) {
            // reset all the probability
            gal_prob_.assign(n_gal_, -1);

            // finding neighbor galaxies around groups
            vector<array<double, 3>> grp_pos;
            vector<array<double, 2>> grp_search;  // maximum searching radius
            for (auto &grp: grps) {
                grp_pos.push_back(grp.pos);
                grp_search.push_back(
                    {f_rvir_ * grp.r_vir, f_sigv_ * grp.sig_v / kH_ / 100});
            }
            auto nbr_finder = NbrFinder3d(
                gal_pos_,
                array<int, 3>{100, 100, 100},
                array<double, 3>{0, 0, 0},
                array<double, 3>{period_, period_, period_});
            auto nbrs = nbr_finder.get_nbr_ids(grp_pos, grp_search);

            // loop all the groups
            for (size_t igrp = 0; igrp < grps.size(); ++igrp) {
                auto &grp = grps[igrp];

                // skip groups with no members
                // (members were assigned to more massive groups)
                if (grp.mems.size() == 0)
                    continue;

                // iteration until no more changes for group members
                set<int> prev_mems;
                int it_mem = 0;
                double grp_hm = grp.hm;
                double grp_z = grp.z;
                array<double, 3> grp_pos = grp.pos;
                while (prev_mems != grp.mems && it_mem < max_mem_iter_) {
                    prev_mems = grp.mems;
                    ++it_mem;
                    // loop all the neighbor galaxies
                    for (auto &nbr: nbrs[igrp]) {
                        int i_gal = nbr.first;
                        double rp = Dist(grp_pos, gal_pos_[nbr.first]);
                        double p = Prob(grp_hm, rp, grp_z, gal_z_[nbr.first]);
                        // for beighbors above the threshold
                        if (p > p_bg_ && p > gal_prob_[i_gal]) {
                            // 1. delete this gal from previous group
                            // 2. insert it to current group
                            // 3. update gal_prob_ and gal_igrp_
                            grps[gal_igrp_[i_gal]].mems.erase(i_gal);
                            grp.mems.insert(i_gal);
                            gal_prob_[i_gal] = p;
                            gal_igrp_[i_gal] = igrp;
                        }
                    }
                    // update halo information using new members
                    UpdateHalo(grp.mems, grp_hm, grp_pos, grp_z);
                }
            }
            gal_igrp_ = RegulateGrp(grps);
        }
    }

    // Update halo mass, pos, z using member galaxies
    void UpdateHalo(const set<int> &mems, double &hm, array<double, 3> &pos, double &z)
    {
        if (mems.size() == 0)
            return;
        double tot_mass = 0;
        pos = {0, 0, 0};
        z = 0;
        for (auto i: mems) {
            double mass = pow(10.0, gal_sm_[i]);
            tot_mass += mass;
            for (int j = 0; j < 3; ++j)
                pos[j] += gal_pos_[i][j] * mass;
            z += gal_z_[i] * mass;
        }
        z /= tot_mass;
        for (int j = 0; j < 3; ++j)
            pos[j] /= tot_mass;
        hm = totsm2hm_.Interp(log10(tot_mass + 1e-9));
    }

    // Regulate groups
    // 0. input grps must have valid mem information (only)
    // 1. eliminate groups with zero members
    // 2. update halo properties
    // 3. sort groups according to halo mass
    // 4. update gal_igrp
    vector<int> RegulateGrp(vector<Group> &grps)
    {
        // eliminate invalid groups
        // calculate total stellar mass
        vector<Group> grps_new;
        vector<double> arr_totsm;
        for (auto &grp: grps) {
            if (grp.mems.size() == 0)
                continue;
            double tot_sm = 0;
            for (auto i: grp.mems)
                tot_sm += pow(10.0, gal_sm_[i]);
            grp.tot_sm = log10(tot_sm + 1e-9);
            arr_totsm.push_back(grp.tot_sm);
            grps_new.push_back(grp);
        }
        grps = grps_new;

        // re-calibrate M_{*, tot} - M_h relation
        CaliTotSm2Hm(arr_totsm);

        // update halo mass, halo pos, and halo z
        // update halo properties
        for (int igrp = 0; igrp < grps.size(); ++igrp) {
            auto &grp = grps[igrp];
            UpdateHalo(grp.mems, grp.hm, grp.pos, grp.z);
            grp.r_vir = Mass2Rvir(grp.hm, grp.z);
            grp.sig_v = Mass2SigV(grp.hm, grp.z);
        }
        // sort grps according to halo mass
        sort(grps.begin(), grps.end(), [](Group a, Group b) { return a.hm > b.hm; });
        // assign gal_igrp
        vector<int> gal_igrp(n_gal_, -1);
        for (int i = 0; i < grps.size(); ++i)
            for (auto i_mem: grps[i].mems)
                gal_igrp_[i_mem] = i;

        return gal_igrp;
    }

    // calibrate the M_{*, tot} - M_h relation
    // output a interpolator
    LinearInterp CaliTotSm2Hm(const vector<double> &tot_sm)
    {
        auto hm = SHAM(tot_sm);
        auto f1 = LinearInterp(tot_sm, hm);
        vector<double> totsm_samp(np_totsm2hm_), hm_samp(np_totsm2hm_);
        int step = int(floor(tot_sm.size() / np_totsm2hm_)) - 1;
        for (int i = 0; i < np_totsm2hm_; ++i) {
            totsm_samp[i] = tot_sm[i * step];
            hm_samp[i] = hm[i * step];
        }
        *totsm_samp.end() = *tot_sm.end();
        *hm_samp.end() = *hm.end();
        return LinearInterp(totsm_samp, hm_samp);
    }

    // abundance match proxy with cumulative halo mass function
    // output order is the same as the rank of proxy
    vector<double> SHAM(const vector<double> &proxy)
    {
        vector<double> hm(proxy.size(), 0);
        vector<int> i_sort(proxy.size(), 0);
        iota(i_sort.begin(), i_sort.end(), 1);

        sort(i_sort.begin(), i_sort.end(), [proxy](int i, int j) {
            return proxy[i] > proxy[j];
        });
        for (size_t i = 0; i < hm.size(); ++i)
            hm[i_sort[i]] = hmf_[i];
        return hm;
    }

    double Prob(double hm, double rp, double z_grp, double z_gal)
    {
        double p;
        double conc = Mass2Conc(hm, z_grp);
        double sig_v = Mass2SigV(hm, z_grp);
        double r_vir = Mass2Rvir(hm, z_grp);
        p = kH_ * 100 / kC_ * SigmaR(rp, r_vir, conc) * GaussianZ(z_grp, z_gal, sig_v);
        return p;
    }

    // TODO:
    inline double Mass2Conc(double hm, double z)
    {
        double c = 10;
        return c;
    }

    inline double Mass2Rvir(double hm, double z)
    {
        // This formula comes from Yang et al. 2021, corrected for r_200
        hm = pow(10.0, hm);
        double r = 0.781 / kH_ * pow(hm / 1e14 * kH_ / kOmega_m_ * 180 / 200, 1.0 / 3);
        r /= (1 + z);
        return r;
    }

    inline double Mass2SigV(double hm, double z)
    {
        hm = pow(10.0, hm);
        double sig_v = 632 * pow(hm * kOmega_m_ * kH_ / 1e14, 0.3224);
        return sig_v;
    }

    inline double SigmaR(double r_p, double r_vir, double conc)
    {
        if (r_p < 1e-3)
            return 1e10;
        double delta = 60 * conc * conc * conc / (log(1 + conc) - conc / (1 + conc));
        double x = r_p * conc / r_vir;
        double x2m = x * x - 1;
        double f;
        if (x == 1)
            f = 1.0 / 3.0;
        else if (x < 1)
            f = (1.0 - log((1.0 + sqrt(-x2m)) / x) / sqrt(-x2m)) / x2m;
        else
            f = (1 - atan(sqrt(x2m)) / sqrt(x2m)) / x2m;

        double sigma = 2.0 * r_vir / conc * delta * f;
        return sigma;
    }

    inline double GaussianZ(double z_grp, double z_gal, double sig_v)
    {
        double zp1 = z_grp + 1;
        double dz = z_gal - z_grp;
        double p = kC_ / sig_v / zp1 * exp(-0.5 * pow(kC_ * dz / sig_v / zp1, 2));
        return p;
    }

    inline double Dist(const array<double, 3> &x, const array<double, 3> &y)
    {
        double d = 0;
        for (int i = 0; i < 3; ++i)
            d += (x[i] - y[i]) * (x[i] - y[i]);
        return sqrt(d);
    }

   private:
    // input
    const vector<array<double, 3>> &gal_pos_;
    const vector<double> &gal_sm_;
    const vector<double> &gal_z_;
    const vector<double> &hmf_;
    const double period_;
    const double p_bg_;
    // output
    vector<double> host_hm_;
    vector<int> if_mm_;
    // ancillary
    vector<double> gal_prob_;
    vector<int> gal_igrp_;
    int n_gal_;
    LinearInterp totsm2hm_;
    // constants: light speed, Hubble constant
    const double kC_ = 3e5;
    const double kH_ = 0.67;
    const double kOmega_m_ = 0.3;
    // number of points for interpolating totsm and hm relation
    const int np_totsm2hm_ = 100;
    // maximum number of iterations in group member assignment
    const int max_mem_iter_ = 20;
    // neighbor searching range
    const double f_rvir_ = 3;
    const double f_sigv_ = 3;
};
