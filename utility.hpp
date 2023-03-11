#ifndef UTILITY_CPP_
#define UTILITY_CPP_

#include <hippcntl.h>
#include <hippio.h>
#include <iostream>
#include <numeric>
#include <vector>

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
#endif /* end of include guard: UTILITY_CPP_ */
