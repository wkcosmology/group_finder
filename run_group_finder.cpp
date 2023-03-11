#include <hippcntl.h>
#include <hippio.h>
#include <iostream>
#include <vector>

#include "hbf_yang.hpp"

int main(void)
{
    const string gal_file = "";
    const string hmf_file = "";
    const string out_file = "";
    double period = 205;
    double vol = period * period * period;
    double p_bg = 10;

    // load galaxy catalog
    vector<array<double, 3>> pos;
    vector<double> sm, z, arr_pos;
    auto f = H5File(gal_file, "r");
    f.open_dataset("Pos").read(arr_pos);
    f.open_dataset("z").read(z);
    f.open_dataset("sm").read(sm);
    pos.resize(z.size());
    for (int i = 0; i < z.size(); ++i)
        pos.push_back({arr_pos[3 * i], arr_pos[3 * i + 1], arr_pos[3 * i + 2]});
    size_t n_gal = z.size();

    // load halo mass function
    auto hmf = HaloMassFunc(hmf_file, vol);

    // run group finder
    auto hmf_vals = hmf.GenerateHMSamp(pos.size());
    auto finder = YangFinder(pos, z, sm, hmf_vals, period, p_bg);

    // output results
    auto f_out = H5File(out_file, "w");
    f_out.create_dataset<double>("GroupID", {n_gal}).write(finder.GetGrpID().data());
    f_out.create_dataset<double>("HostHM", {n_gal}).write(finder.GetHaloMass().data());
    f_out.create_dataset<double>("IfMM", {n_gal}).write(finder.GetIfMM().data());
    f_out.create_dataset<double>("GroupPos", {n_gal, 3}).write(finder.GetGrpPos().data());
    f_out.create_dataset<double>("Z", {n_gal}).write(finder.GetGrpZ().data());

    return 0;
}
