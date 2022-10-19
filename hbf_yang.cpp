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
#include <vector>

#include "nbr_finder.hpp"

using namespace std;

class YangFinder
{
   private:
    struct Group {
        double z, tot_sm, hm, conc, r_vir, sig_v;
        set<int> mems;
        array<double, 3> pos;
    };

   public:
    /**
     * @brief Init the halo-based group finder
     *
     * @param pos [positions, unit (Mpc, Mpc, redshift), z-direction is the
     * third one]
     * @param log_sm [log of stellar mass]
     * @param hmf [halo mass array, generated from halo mass function]
     * @param z [redshift for group searching]
     * @param period [box size for periodic box]
     * @param rp_max [maximum searching radius on the projection direction,
     * unit: Mpc]
     * @param pi_max [maximum searching radius on the los direction, unit: Mpc]
     */
    YangFinder(
        const vector<array<double, 3>> &pos,
        const vector<double> &z,
        const vector<double> &sm,
        const vector<double> &hmf,
        double period,
        double rp_max,
        double pi_max,
        double p_bg)
        : gal_pos_(pos),
          gal_sm_(sm),
          hmf_(hmf),
          gal_z_(z),
          period_(period),
          rp_max_(rp_max),
          pi_max_(pi_max),
          p_bg_(p_bg)
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
                for (auto &i: p) {
                    if (i < 0 || i > period_) {
                        cout << "Position out of range" << i << endl;
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
            gal_prob_.assign(n_gal_, -1);
            vector<array<double, 3>> grp_pos;
            vector<array<double, 2>> grp_search;  // maximum searching radius
            for (auto &grp: grps) {
                grp_pos.push_back(grp.pos);
                grp_search.push_back({grp.r_vir, grp.sig_v / H0_});
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

                // loop all the neighbor galaxies around given group
                for (auto &nbr: nbrs[igrp]) {
                    int i_gal = nbr.first;
                    double p =
                        Prob(grp.hm, nbr.second[0], grp.z, gal_z_[nbr.first]);
                    // skip when the probability is too low
                    if (p < p_bg_ || p < gal_prob_[i_gal])
                        continue;
                    // 1. delete this gal from previous group
                    // 2. insert it to current group
                    // 3. update gal_prob_ and gal_igrp_
                    grps[gal_igrp_[i_gal]].mems.erase(i_gal);
                    grp.mems.insert(i_gal);
                    gal_prob_[i_gal] = p;
                    gal_igrp_[i_gal] = igrp;
                }

                // !!! Since we expect very large prob for central, for group
                // with only one member, we should make it vulnerable to be
                // included by other groups. Otherwise, these centrals will
                // never be assigned to other groups
                if (grp.mems.size() == 1)
                    gal_prob_[*grp.mems.begin()] = -1;
            }
            gal_igrp_ = RegulateGrp(grps);
        }
    }

    // Regulate groups
    // 0. input grps must have valid mem information (only)
    // 1. eliminate unreliable group
    // 2. sort groups according to total stellar mass or halo mass
    // 3. calculate halo properties
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
        // performance SHAM
        // assign halo mass
        // calculate halo related quantities
        auto arr_hm = SHAM(arr_totsm);
        for (int igrp = 0; igrp < grps.size(); ++igrp) {
            auto &grp = grps[igrp];
            double hm = arr_hm[igrp];
            double z = grp.z;
            grp.hm = hm;
            grp.conc = Mass2Conc(hm, z);
            grp.r_vir = Mass2Rvir(hm, z);
            grp.sig_v = Mass2SigV(hm, z);
            // TODO: update pos and redshift
            array<double, 3> pos_tmp;
            double z_tmp = 0;
            double tot_sm_tmp = 0;
            for (auto i_mem: grp.mems) {
                double sm_tmp = pow(10.0, gal_sm_[i_mem] - 10);
                tot_sm_tmp += sm_tmp;
                for (int i = 0; i < 3; ++i)
                    pos_tmp[i] += sm_tmp * gal_pos_[i_mem][i];
                z_tmp += sm_tmp * gal_z_[i_mem];
            }
            for (int i = 0; i < 3; ++i)
                grp.pos[i] = pos_tmp[i] / tot_sm_tmp;
            grp.z = z_tmp / tot_sm_tmp;
        }
        // sort grps according to halo mass
        sort(grps.begin(), grps.end(), [](Group a, Group b) {
            return a.hm > b.hm;
        });
        // assign gal_igrp
        vector<int> gal_igrp(n_gal_, -1);
        for (int i = 0; i < grps.size(); ++i)
            for (auto i_mem: grps[i].mems)
                gal_igrp_[i_mem] = i;

        return gal_igrp;
    }

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
        p = H0_ / c_ * SigmaR(rp, r_vir, conc) * GaussianZ(z_grp, z_gal, sig_v);
        return p;
    }

    // TODO:
    inline double Mass2Conc(double hm, double z)
    {
        double c = 10;
        return c;
    }

    // TODO:
    inline double Mass2Rvir(double hm, double z)
    {
        double r;
        return r;
    }

    // TODO:
    inline double Mass2SigV(double hm, double z)
    {
        double sig_v;
        return sig_v;
    }

    inline double SigmaR(double r_p, double r_vir, double conc)
    {
        if (r_p < 1e-3)
            return 1e10;
        double delta =
            60 * conc * conc * conc / (log(1 + conc) - conc / (1 + conc));
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
        double p = c_ / sig_v / zp1 * exp(-0.5 * pow(c_ * dz / sig_v / zp1, 2));
        return p;
    }

   private:
    // input
    const vector<array<double, 3>> &gal_pos_;
    const vector<double> &gal_sm_;
    const vector<double> &gal_z_;
    const vector<double> &hmf_;
    const double period_;
    const double rp_max_;
    const double pi_max_;
    const double p_bg_;
    // output
    vector<double> host_hm_;
    vector<int> if_mm_;
    // ancillary
    vector<double> gal_prob_;
    vector<int> gal_igrp_;
    int n_gal_;
    const double c_ = 3e5;  // light speed, unit: km/s
    const double H0_ = 67;  // Hubble constant, unit: km/s/Mpc
};
