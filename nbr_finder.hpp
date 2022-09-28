#ifndef NBR_FINDER_HPP_
#define NBR_FINDER_HPP_

#include <algorithm>
#include <array>
#include <boost/multi_array.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>

using std::array;
using std::ostream;
using std::pair;
using std::vector;

/**
 * @brief Finding the neighbour points in a reference sample, in 3-D space
 * The class is instantiated with the reference sample.
 * The query of points around certain position is achieved though class
 * method NbrFinder3d::get_nbr_ids.
 *
 * Support query methods:
 * 1-D: sqrt(dr0 * dr0 + dr1 * dr1 + dr2 * dr2) < r_max
 * 2-D: sqrt(dr0 * dr0 + dr1 * dr1) < r_max[0] && |dr2| < r_max[1]
 * 3-D: |dr0| < r_max[0] && |dr1| < r_max[1] && |dr2| < r_max[2]
 * Each r_max could vary with query positions
 */
class NbrFinder3d
{
   public:
    typedef pair<int, double> IdDistPair1d;
    typedef pair<int, array<double, 2>> IdDistPair2d;
    typedef pair<int, array<double, 3>> IdDistPair3d;
    typedef boost::multi_array<int, 3> ArrInt3d;
    typedef array<double, 3> Pos;
    typedef array<int, 3> CellId;

    /**
     * @brief Constructor for periodic box
     */
    NbrFinder3d(
        const vector<Pos> &ps,
        array<int, 3> num_cell,
        array<double, 3> edge_l,
        array<double, 3> edge_u,
        bool verbose = false);
    /**
     * @brief Constructor for non-periodic box
     */
    NbrFinder3d(
        const vector<Pos> &ps, array<int, 3> num_cell, bool verbose = false);

    ostream &summary(ostream &out);

    vector<vector<IdDistPair1d>> get_nbr_ids(
        const vector<Pos> &ps, double max_r)
    {
        return get_nbr_ids_helper(ps, vector<double>{max_r}, true);
    }
    vector<vector<IdDistPair2d>> get_nbr_ids(
        const vector<Pos> &ps, array<double, 2> max_r)
    {
        return get_nbr_ids_helper(ps, vector<array<double, 2>>{max_r}, true);
    }
    vector<vector<IdDistPair3d>> get_nbr_ids(
        const vector<Pos> &ps, array<double, 3> max_r)
    {
        return get_nbr_ids_helper(ps, vector<array<double, 3>>{max_r}, true);
    }

    vector<vector<IdDistPair1d>> get_nbr_ids(
        const vector<Pos> &ps, const vector<double> &max_r)
    {
        return get_nbr_ids_helper(ps, max_r, false);
    }
    vector<vector<IdDistPair2d>> get_nbr_ids(
        const vector<Pos> &ps, const vector<array<double, 2>> &max_r)
    {
        return get_nbr_ids_helper(ps, max_r, false);
    }
    vector<vector<IdDistPair3d>> get_nbr_ids(
        const vector<Pos> &ps, const vector<array<double, 3>> &max_r)
    {
        return get_nbr_ids_helper(ps, max_r, false);
    }

   private:
    bool verbose_;
    int size_;
    bool if_periodic_;

    vector<array<double, 3>> pos_;
    array<double, 3> edge_l_;
    array<double, 3> edge_u_;
    array<double, 3> box_size_;
    array<int, 3> num_cell_;
    array<double, 3> cell_len_;
    vector<int> index_;
    ArrInt3d first_id_;

    vector<vector<IdDistPair1d>> get_nbr_ids_helper(
        const vector<Pos> &ps, const vector<double> &max_r, bool if_uniform);
    vector<vector<IdDistPair2d>> get_nbr_ids_helper(
        const vector<Pos> &ps,
        const vector<array<double, 2>> &max_r,
        bool if_uniform);
    vector<vector<IdDistPair3d>> get_nbr_ids_helper(
        const vector<Pos> &ps,
        const vector<array<double, 3>> &max_r,
        bool if_uniform);
    void traverse_cell(
        Pos p_pos,
        const vector<int> &ps_in_cell,
        double max_r,
        vector<IdDistPair1d> &res_p,
        bool if_periodic);
    void traverse_cell(
        Pos p_pos,
        const vector<int> &ps_in_cell,
        array<double, 2> max_r,
        vector<IdDistPair2d> &res_p,
        bool if_periodic);
    void traverse_cell(
        Pos p_pos,
        const vector<int> &ps_in_cell,
        array<double, 3> max_r,
        vector<IdDistPair3d> &res_p,
        bool if_periodic);

    vector<CellId> pos2cell_id_vec(const vector<Pos> &ps);
    void build_box();
    inline vector<int> box_id2p_id(CellId cell_id);
    inline double dist(const vector<double> &p1, const vector<double> &p2);
    inline double norm(const vector<double> &p1);
};

vector<int> NbrFinder3d::box_id2p_id(NbrFinder3d::CellId cell_id)
{
    vector<int> p_id;
    int id = first_id_[cell_id[0]][cell_id[1]][cell_id[2]];
    while (id != -1) {
        p_id.push_back(id);
        id = index_[id];
    }
    return p_id;
}

double NbrFinder3d::norm(const vector<double> &p1)
{
    double dist = 0;
    for (int i = 0; i < p1.size(); ++i)
        dist += p1[i] * p1[i];
    return sqrt(dist);
}

#endif /* end of include guard: NBR_FINDER_HPP_ */
