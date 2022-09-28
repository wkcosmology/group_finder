#include "nbr_finder.hpp"

#include <array>
#include <boost/multi_array.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

using std::abs;
using std::array;
using std::ceil;
using std::cout;
using std::endl;
using std::floor;
using std::max_element;
using std::ostream;
using std::vector;

NbrFinder3d::NbrFinder3d(
    const vector<NbrFinder3d::Pos> &ps, array<int, 3> num_cell, bool verbose)
    : verbose_(verbose),
      size_(ps.size()),
      if_periodic_(false),
      pos_(ps),
      num_cell_(num_cell)
{
    double inf = std::numeric_limits<double>::infinity();
    edge_l_.fill(inf);
    edge_u_.fill(-inf);
    for (int i_p = 0; i_p < ps.size(); ++i_p) {
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
            if (ps[i_p][i_dim] < edge_l_[i_dim])
                edge_l_[i_dim] = ps[i_p][i_dim] - 1e-10;
            if (ps[i_p][i_dim] > edge_u_[i_dim])
                edge_u_[i_dim] = ps[i_p][i_dim] + 1e-10;
        }
    }
    for (int i = 0; i < 3; ++i) {
        cell_len_[i] = (edge_u_[i] - edge_l_[i]) / num_cell_[i];
        box_size_[i] = edge_u_[i] - edge_l_[i];
    }
    if (verbose_)
        summary(cout);
    build_box();
}

NbrFinder3d::NbrFinder3d(
    const vector<NbrFinder3d::Pos> &ps,
    array<int, 3> num_cell,
    array<double, 3> edge_l,
    array<double, 3> edge_u,
    bool verbose)
    : verbose_(verbose),
      size_(ps.size()),
      if_periodic_(true),
      pos_(ps),
      edge_l_(edge_l),
      edge_u_(edge_u),
      num_cell_(num_cell)
{
    for (int i = 0; i < 3; ++i) {
        if (edge_l_[i] >= edge_u_[i]) {
            cout << "lower limit >= upper limit, exit" << endl;
            exit(EXIT_FAILURE);
        }
        cell_len_[i] = (edge_u_[i] - edge_l_[i]) / num_cell_[i];
        box_size_[i] = edge_u_[i] - edge_l_[i];
    }
    if (verbose_)
        summary(cout);
    build_box();
}

void NbrFinder3d::build_box()
{
    index_.assign(size_, -1);
    first_id_.resize(boost::extents[num_cell_[0]][num_cell_[1]][num_cell_[2]]);
    for (int i = 0; i < num_cell_[0]; ++i)
        for (int j = 0; j < num_cell_[1]; ++j)
            for (int k = 0; k < num_cell_[2]; ++k)
                first_id_[i][j][k] = -1;
    // build the box
    if (verbose_)
        cout << ">>> NbrFinder3d: build the index box" << endl;
    auto cell_ids = pos2cell_id_vec(pos_);
    for (int i = 0; i < pos_.size(); ++i) {
        auto &box_id = cell_ids[i];
        int &first_id = first_id_[box_id[0]][box_id[1]][box_id[2]];
        index_[i] = first_id;
        first_id = i;
    }
    if (verbose_)
        cout << "<<< NbrFinder3d: Finish building the index box" << endl;
}

ostream &NbrFinder3d::summary(ostream &out)
{
    out << "======Summary of NbrFinder3d======\n";
    out << "Num of input points: " << pos_.size() << "\n";
    out << "Box edge: ";
    for (int i = 0; i < 3; ++i)
        out << "(" << edge_l_[i] << ", " << edge_u_[i] << ") x ";
    out << "\n";
    out << "Num of cells: ";
    for (int i = 0; i < 3; ++i)
        out << num_cell_[i] << " x ";
    out << "\n";
    out << "cell side length: ";
    for (int i = 0; i < 3; ++i)
        out << cell_len_[i] << " x ";
    out << "\n";
    out << "======End of Summary======" << endl;

    return out;
}

vector<vector<NbrFinder3d::IdDistPair1d>> NbrFinder3d::get_nbr_ids_helper(
    const vector<NbrFinder3d::Pos> &ps,
    const vector<double> &max_r,
    bool if_uniform)
{
    vector<vector<IdDistPair1d>> res(ps.size());
    if (max_r.size() == 0)
        return res;

    array<int, 3> d_cell;
    if (if_uniform) {
        d_cell[0] = (int)ceil(max_r[0] / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r[0] / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r[0] / cell_len_[2]);
    } else {
        double max_r0 = *max_element(max_r.cbegin(), max_r.cend());
        d_cell[0] = (int)ceil(max_r0 / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r0 / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r0 / cell_len_[2]);
    }
    auto ps_cell_ids = pos2cell_id_vec(ps);
    for (int i_p = 0; i_p < ps.size(); ++i_p) {
        for (int i_x = -d_cell[0]; i_x < d_cell[0] + 1; ++i_x)
            for (int i_y = -d_cell[1]; i_y < d_cell[1] + 1; ++i_y)
                for (int i_z = -d_cell[2]; i_z < d_cell[2] + 1; ++i_z) {
                    double max_r0 = if_uniform ? max_r[0] : max_r[i_p];
                    if (max_r0
                        < norm(
                              {i_x * cell_len_[0],
                               i_y * cell_len_[1],
                               i_z * cell_len_[2]})
                              - norm(
                                  {cell_len_[0], cell_len_[1], cell_len_[2]}))
                        continue;
                    int x_id = ps_cell_ids[i_p][0] + i_x;
                    int y_id = ps_cell_ids[i_p][1] + i_y;
                    int z_id = ps_cell_ids[i_p][2] + i_z;
                    if (if_periodic_) {
                        x_id = x_id < 0 ? x_id + num_cell_[0]
                                        : x_id % num_cell_[0];
                        y_id = y_id < 0 ? y_id + num_cell_[1]
                                        : y_id % num_cell_[1];
                        z_id = z_id < 0 ? z_id + num_cell_[2]
                                        : z_id % num_cell_[2];
                    } else {
                        if (x_id < 0 || x_id >= num_cell_[0] || y_id < 0
                            || y_id >= num_cell_[1] || z_id < 0
                            || z_id >= num_cell_[2])
                            continue;
                    }
                    auto ps_in_cell = box_id2p_id({x_id, y_id, z_id});
                    traverse_cell(
                        ps[i_p], ps_in_cell, max_r0, res[i_p], if_periodic_);
                }
    }
    return res;
}

vector<vector<NbrFinder3d::IdDistPair2d>> NbrFinder3d::get_nbr_ids_helper(
    const vector<NbrFinder3d::Pos> &ps,
    const vector<array<double, 2>> &max_r,
    bool if_uniform)
{
    vector<vector<IdDistPair2d>> res(ps.size());
    if (max_r.size() == 0)
        return res;
    array<int, 3> d_cell;
    if (if_uniform) {
        d_cell[0] = (int)ceil(max_r[0][0] / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r[0][0] / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r[0][1] / cell_len_[2]);
    } else {
        double max_r0 = 0, max_r1 = 0;
        for (int i = 0; i < max_r.size(); ++i) {
            if (max_r[i][0] > max_r0)
                max_r0 = max_r[i][0];
            if (max_r[i][1] > max_r1)
                max_r1 = max_r[i][1];
        }
        d_cell[0] = (int)ceil(max_r0 / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r0 / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r1 / cell_len_[2]);
    }
    auto ps_cell_ids = pos2cell_id_vec(ps);
    for (int i_p = 0; i_p < ps.size(); ++i_p) {
        for (int i_x = -d_cell[0]; i_x < d_cell[0] + 1; ++i_x)
            for (int i_y = -d_cell[1]; i_y < d_cell[1] + 1; ++i_y)
                for (int i_z = -d_cell[2]; i_z < d_cell[2] + 1; ++i_z) {
                    double max_r0 = if_uniform ? max_r[0][0] : max_r[i_p][0];
                    double max_r1 = if_uniform ? max_r[0][1] : max_r[i_p][1];
                    if (max_r0 < norm({i_x * cell_len_[0], i_y * cell_len_[1]})
                                     - norm({cell_len_[0], cell_len_[1]}))
                        continue;
                    int x_id = ps_cell_ids[i_p][0] + i_x;
                    int y_id = ps_cell_ids[i_p][1] + i_y;
                    int z_id = ps_cell_ids[i_p][2] + i_z;
                    if (if_periodic_) {
                        x_id = x_id < 0 ? x_id + num_cell_[0]
                                        : x_id % num_cell_[0];
                        y_id = y_id < 0 ? y_id + num_cell_[1]
                                        : y_id % num_cell_[1];
                        z_id = z_id < 0 ? z_id + num_cell_[2]
                                        : z_id % num_cell_[2];
                    } else {
                        if (x_id < 0 || x_id >= num_cell_[0] || y_id < 0
                            || y_id >= num_cell_[1] || z_id < 0
                            || z_id >= num_cell_[2])
                            continue;
                    }
                    auto ps_in_cell = box_id2p_id({x_id, y_id, z_id});
                    traverse_cell(
                        ps[i_p],
                        ps_in_cell,
                        {max_r0, max_r1},
                        res[i_p],
                        if_periodic_);
                }
    }
    return res;
}

vector<vector<NbrFinder3d::IdDistPair3d>> NbrFinder3d::get_nbr_ids_helper(
    const vector<NbrFinder3d::Pos> &ps,
    const vector<array<double, 3>> &max_r,
    bool if_uniform)
{
    vector<vector<IdDistPair3d>> res(ps.size());
    if (max_r.size() == 0)
        return res;
    array<int, 3> d_cell;
    if (if_uniform) {
        d_cell[0] = (int)ceil(max_r[0][0] / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r[0][1] / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r[0][2] / cell_len_[2]);
    } else {
        double max_r0 = 0, max_r1 = 0, max_r2 = 0;
        for (int i = 0; i < max_r.size(); ++i) {
            if (max_r[i][0] > max_r0)
                max_r0 = max_r[i][0];
            if (max_r[i][1] > max_r1)
                max_r1 = max_r[i][1];
            if (max_r[i][2] > max_r2)
                max_r2 = max_r[i][2];
        }
        d_cell[0] = (int)ceil(max_r0 / cell_len_[0]);
        d_cell[1] = (int)ceil(max_r1 / cell_len_[1]);
        d_cell[2] = (int)ceil(max_r2 / cell_len_[2]);
    }
    auto ps_cell_ids = pos2cell_id_vec(ps);
    for (int i_p = 0; i_p < ps.size(); ++i_p) {
        for (int i_x = -d_cell[0]; i_x < d_cell[0] + 1; ++i_x)
            for (int i_y = -d_cell[1]; i_y < d_cell[1] + 1; ++i_y)
                for (int i_z = -d_cell[2]; i_z < d_cell[2] + 1; ++i_z) {
                    double max_r0 = if_uniform ? max_r[0][0] : max_r[i_p][0];
                    double max_r1 = if_uniform ? max_r[0][1] : max_r[i_p][1];
                    double max_r2 = if_uniform ? max_r[0][2] : max_r[i_p][2];
                    int x_id = ps_cell_ids[i_p][0] + i_x;
                    int y_id = ps_cell_ids[i_p][1] + i_y;
                    int z_id = ps_cell_ids[i_p][2] + i_z;
                    if (if_periodic_) {
                        x_id = x_id < 0 ? x_id + num_cell_[0]
                                        : x_id % num_cell_[0];
                        y_id = y_id < 0 ? y_id + num_cell_[1]
                                        : y_id % num_cell_[1];
                        z_id = z_id < 0 ? z_id + num_cell_[2]
                                        : z_id % num_cell_[2];
                    } else {
                        if (x_id < 0 || x_id >= num_cell_[0] || y_id < 0
                            || y_id >= num_cell_[1] || z_id < 0
                            || z_id >= num_cell_[2])
                            continue;
                    }
                    auto ps_in_cell = box_id2p_id({x_id, y_id, z_id});
                    traverse_cell(
                        ps[i_p],
                        ps_in_cell,
                        {max_r0, max_r1, max_r2},
                        res[i_p],
                        if_periodic_);
                }
    }
    return res;
}

vector<NbrFinder3d::CellId> NbrFinder3d::pos2cell_id_vec(
    const vector<NbrFinder3d::Pos> &ps)
{
    vector<CellId> cell_ids(ps.size());
    for (int i = 0; i < ps.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            cell_ids[i][j] = (int)floor((ps[i][j] - edge_l_[j]) / cell_len_[j]);
    }
    return cell_ids;
}

void NbrFinder3d::traverse_cell(
    Pos p_pos,
    const vector<int> &ps_in_cell,
    double max_r,
    vector<IdDistPair1d> &res_p,
    bool if_periodic)
{
    double dist0;
    array<double, 3> delta;
    if (if_periodic) {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i) {
                delta[i] = abs(pos_[p][i] - p_pos[i]);
                delta[i] = delta[i] < box_size_[i] / 2
                               ? delta[i]
                               : box_size_[i] - delta[i];
            }
            dist0 = sqrt(
                delta[0] * delta[0] + delta[1] * delta[1]
                + delta[2] * delta[2]);
            if (dist0 < max_r)
                res_p.push_back({p, dist0});
        }
    } else {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i)
                delta[i] = abs(pos_[p][i] - p_pos[i]);
            dist0 = sqrt(
                delta[0] * delta[0] + delta[1] * delta[1]
                + delta[2] * delta[2]);
            if (dist0 < max_r)
                res_p.push_back({p, dist0});
        }
    }
}

void NbrFinder3d::traverse_cell(
    Pos p_pos,
    const vector<int> &ps_in_cell,
    array<double, 2> max_r,
    vector<IdDistPair2d> &res_p,
    bool if_periodic)
{
    double dist0, dist1;
    array<double, 3> delta;
    if (if_periodic) {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i) {
                delta[i] = abs(pos_[p][i] - p_pos[i]);
                delta[i] = delta[i] < box_size_[i] / 2
                               ? delta[i]
                               : box_size_[i] - delta[i];
            }
            dist0 = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
            dist1 = delta[2];
            if (dist0 < max_r[0] && dist1 < max_r[1])
                res_p.push_back({p, {dist0, dist1}});
        }
    } else {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i)
                delta[i] = abs(pos_[p][i] - p_pos[i]);
            dist0 = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
            dist1 = delta[2];
            if (dist0 < max_r[0] && abs(pos_[p][2] - p_pos[2]) < max_r[1])
                res_p.push_back({p, {dist0, dist1}});
        }
    }
}

void NbrFinder3d::traverse_cell(
    Pos p_pos,
    const vector<int> &ps_in_cell,
    array<double, 3> max_r,
    vector<IdDistPair3d> &res_p,
    bool if_periodic)
{
    array<double, 3> delta;
    if (if_periodic) {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i) {
                delta[i] = abs(pos_[p][i] - p_pos[i]);
                delta[i] = delta[i] < box_size_[i] / 2
                               ? delta[i]
                               : box_size_[i] - delta[i];
            }
            if (delta[0] < max_r[0] && delta[1] < max_r[1]
                && delta[2] < max_r[2])
                res_p.push_back({p, delta});
        }
    } else {
        for (auto p: ps_in_cell) {
            for (int i = 0; i < 3; ++i)
                delta[i] = abs(pos_[p][i] - p_pos[i]);
            if (delta[0] < max_r[0] && delta[1] < max_r[1]
                && delta[2] < max_r[2])
                res_p.push_back({p, delta});
        }
    }
}
