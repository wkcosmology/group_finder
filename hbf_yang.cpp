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
#include <vector>

#include "nbr_finder.hpp"

using namespace std;

class YangFinder
{
   public:
    /**
     * @brief Init the halo-based group finder
     *
     * @param pos [positions, unit (Mpc, Mpc, redshift), z-direction is the third one]
     * @param log_sm [log of stellar mass]
     * @param hmf [halo mass array, generated from halo mass function]
     * @param z [redshift for group searching]
     * @param period [box size for periodic box]
     * @param rp_max [maximum searching radius on the projection direction, unit: Mpc]
     * @param pi_max [maximum searching radius on the los direction, unit: Mpc]
     */
    YangFinder(
        const vector<array<double, 3>> &pos,
        const vector<double> z,
        const vector<double> &log_sm,
        const vector<double> &hmf,
        double period,
        double rp_max,
        double pi_max)
        : pos_(pos),
          log_sm_(log_sm),
          hmf_(hmf),
          z_(z),
          period_(period),
          rp_max_(rp_max),
          pi_max_(pi_max)
    {
        // input validity check
        n_gal_ = pos.size();
        cout << "# of input galaxies: " << n_gal_ << endl;
        if (log_sm_.size() != n_gal_)
            cout << "# of stellar mass is wrong: " << log_sm_.size() << endl;
        if (hmf_.size() != n_gal_)
            cout << "# of halo shape is wrong: " << hmf_.size() << endl;
        if (z_.size() != n_gal_)
            cout << "# of redshift is wrong: " << z_.size() << endl;
        grp_id_.assign(n_gal_, -1);
        host_hm_.assign(n_gal_, -1);
        if_mm_.assign(n_gal_, -1);
        // check position
        if (period_ > 0) {
            for (auto &p: pos_) {
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
        auto nbr_finder = NbrFinder3d(
            pos_,
            array<int, 3>{100, 100, 100},
            array<double, 3>{0, 0, 0},
            array<double, 3>{period_, period_, period_});

        vector<array<double, 3>> grp_pos(pos_);
        for (int it = 0; it < n_iter; ++it) {
            auto nbr = nbr_finder.get_nbr_ids(pos_, array<double, 2>{rp_max_, pi_max_});
            vector<int> cen_id;
            vector<double> host_hm;
        }
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

    double Prob(double log_hm, double rp, double z_grp, double z_gal)
    {
        double p;
        double conc = Mass2Conc(log_hm, z_grp);
        double sig_v = Mass2SigV(log_hm, z_grp);
        double r_vir = Mass2Rvir(log_hm, z_grp);
        p = H0_ / c_ * SigmaR(rp, r_vir, conc) * GaussianZ(z_grp, z_gal, sig_v);
        return p;
    }

    inline double Mass2Conc(double log_hm, double z)
    {
        double c = 10;
        return c;
    }

    inline double Mass2Rvir(double log_hm, double z)
    {
        double r;
        return r;
    }

    inline double Mass2SigV(double log_hm, double z)
    {
        double sig_v;
        return sig_v;
    }

    inline double SigmaR(double r_p, double r_vir, double conc)
    {
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
        double p = c_ / sig_v / zp1 * exp(-0.5 * pow(c_ * dz / sig_v / zp1, 2));
        return p;
    }

   private:
    // input
    const vector<array<double, 3>> &pos_;
    const vector<double> &log_sm_;
    const vector<double> &hmf_;
    const vector<double> z_;
    const double period_;
    const double rp_max_;
    const double pi_max_;
    // output
    vector<int> grp_id_;
    vector<double> host_hm_;
    vector<int> if_mm_;
    // ancillary
    int n_gal_;
    const double c_ = 3e5;  // light speed, unit: km/s
    const double H0_ = 67;  // Hubble constant, unit: km/s/Mpc
};
