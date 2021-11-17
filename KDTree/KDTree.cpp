/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 * https://rosettacode.org/wiki/K-d_tree
 *
 * It is a reimplementation of the C code using C++.  It also includes a few
 * more queries than the original, namely finding all points at a distance
 * smaller than some given distance to a point.
 *
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>
#include <queue>

#include "KDTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h> 
#include <pybind11/chrono.h>

KDNode::KDNode() = default;

KDNode::KDNode(const point_t &pt, const size_t &idx_, const KDNodePtr &left_,
               const KDNodePtr &right_) {
    x = pt;
    index = idx_;
    left = left_;
    right = right_;
}

KDNode::KDNode(const pointIndex &pi, const KDNodePtr &left_,
               const KDNodePtr &right_) {
    x = pi.first;
    index = pi.second;
    left = left_;
    right = right_;
}

KDNode::~KDNode() = default;

double KDNode::coord(const size_t &idx) { return x.at(idx); }
KDNode::operator bool() { return (!x.empty()); }
KDNode::operator point_t() { return x; }
KDNode::operator size_t() { return index; }
KDNode::operator pointIndex() { return pointIndex(x, index); }

KDNodePtr NewKDNodePtr() {
    KDNodePtr mynode = std::make_shared< KDNode >();
    return mynode;
}

inline double dist2(const point_t &a, const point_t &b) {
    double distc = 0;
    for (size_t i = 0; i < a.size(); i++) {
        double di = a.at(i) - b.at(i);
        distc += di * di;
    }
    return distc;
}

inline double dist2(const KDNodePtr &a, const KDNodePtr &b) {
    return dist2(a->x, b->x);
}

inline double dist(const point_t &a, const point_t &b) {
    return std::sqrt(dist2(a, b));
}

inline double dist(const KDNodePtr &a, const KDNodePtr &b) {
    return std::sqrt(dist2(a, b));
}

comparer::comparer(size_t idx_) : idx{idx_} {};

inline bool comparer::compare_idx(const pointIndex &a,  //
                                  const pointIndex &b   //
) {
    return (a.first.at(idx) < b.first.at(idx));  //
}

inline void sort_on_idx(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        size_t idx) {
    comparer comp(idx);
    comp.idx = idx;

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::nth_element(begin, begin + std::distance(begin, end) / 2,
                     end, std::bind(&comparer::compare_idx, comp, _1, _2));
}

void queue_insert(pointQueue &queue, pointDist &item, int k) {
    int s = queue.size();
    if(s < k) {
        queue.push(item);
    } else {
        auto t = queue.top();
        if(t.second > item.second) {
            queue.pop();
            queue.push(item);
        }
    }
}

using pointVec = std::vector< point_t >;

KDNodePtr KDTree::make_tree(const pointIndexArr::iterator &begin,  //
                            const pointIndexArr::iterator &end,    //
                            const size_t &length,                  //
                            const size_t &level                    //
) {
    if (begin == end) {
        return NewKDNodePtr();  // empty tree
    }

    size_t dim = begin->first.size();

    if (length > 1) {
        sort_on_idx(begin, end, level);
    }

    auto middle = begin + (length / 2);

    auto l_begin = begin;
    auto l_end = middle;
    auto r_begin = middle + 1;
    auto r_end = end;

    size_t l_len = length / 2;
    size_t r_len = length - l_len - 1;

    KDNodePtr left;
    if (l_len > 0 && dim > 0) {
        left = make_tree(l_begin, l_end, l_len, (level + 1) % dim);
    } else {
        left = leaf;
    }
    KDNodePtr right;
    if (r_len > 0 && dim > 0) {
        right = make_tree(r_begin, r_end, r_len, (level + 1) % dim);
    } else {
        right = leaf;
    }
    auto node = std::make_shared< KDNode >(*middle, left, right);
    node->level = level;
    left->parent = node;
    left->level = (level + 1) % dim;
    right->parent = node;
    right->level = (level + 1) % dim;
    // KDNode result = KDNode();
    return node;
}

KDTree::KDTree(pointVec point_array) {
    leaf = std::make_shared< KDNode >();
    // iterators
    pointIndexArr arr;
    for (size_t i = 0; i < point_array.size(); i++) {
        arr.push_back(pointIndex(point_array.at(i), i));
    }

    auto begin = arr.begin();
    auto end = arr.end();

    size_t length = arr.size();
    size_t level = 0;  // starting

    root = KDTree::make_tree(begin, end, length, level);
}

KDNodePtr KDTree::nearest_(   //
    const KDNodePtr &branch,  //
    const point_t &pt,        //
    const size_t &level,      //
    const KDNodePtr &best,    //
    const double &best_dist,   //
    pointQueue &nearest_set,
    const int topk=1
) {
    double d, dx, dx2;

    if (!bool(*branch)) {
        return NewKDNodePtr();  // basically, null
    }

    point_t branch_pt(*branch);
    size_t dim = branch_pt.size();

    d = dist2(branch_pt, pt);
    dx = branch_pt.at(level) - pt.at(level);
    dx2 = dx * dx;

    KDNodePtr best_l = best;
    double best_dist_l = best_dist;
    // pointDist pd = pointDist(best, best_dist_l);
    // queue_insert(nearest_set, pd, topk);
    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
        pointDist pd_l = pointDist(best_l, best_dist_l);
        queue_insert(nearest_set, pd_l, topk);
    }

    size_t next_lv = (level + 1) % dim;
    KDNodePtr section;
    KDNodePtr other;

    // select which branch makes sense to check
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    // keep nearest neighbor from further down the tree
    KDNodePtr further = nearest_(section, pt, next_lv, best_l, best_dist_l, nearest_set, topk);
    if (!further->x.empty()) {
        double dl = dist2(further->x, pt);
        if (dl < best_dist_l) {
            best_dist_l = dl;
            best_l = further;
            pointDist pd_l = pointDist(best_l, best_dist_l);
            queue_insert(nearest_set, pd_l, topk);
        }
    }
    // only check the other branch if it makes sense to do so
    if (dx2 < best_dist_l) {
        further = nearest_(other, pt, next_lv, best_l, best_dist_l, nearest_set, topk);
        if (!further->x.empty()) {
            double dl = dist2(further->x, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
                pointDist pd_l = pointDist(best_l, best_dist_l);
                queue_insert(nearest_set, pd_l, topk);
            }
        }
    }

    return best_l;
};



// default caller
KDNodePtr KDTree::nearest_(const point_t &pt, const int topk=1) {
    size_t level = 0;
    // KDNodePtr best = branch;
    double branch_dist = dist2(point_t(*root), pt);
   
    pointQueue nearest_set;
    return nearest_(root,          // beginning of tree
                    pt,            // point we are querying
                    level,         // start from level 0
                    root,          // best is the root
                    branch_dist,
                    nearest_set,
                    topk);  // best_dist = branch_dist
};

bool node_visited(const KDNodePtr &node,  std::set<size_t> visited) {
    bool f = true;
    if(bool(*(node))) {
        if (visited.find(node->index) == visited.end()) {
            f = false;
        }
    }
    return f;
}

void KDTree::view(KDNodePtr node) {
    if(!bool(*node)) {
        return;
    }
    
		std::queue<KDNodePtr> q = std::queue<KDNodePtr>();
		q.push(root);
		while (!q.empty()) {
			auto n = q.front();
            if (n!=node)

            
			if (bool(*n->left)) {
				q.push(n->left);
			}
			if (bool(*n->right)) {
				q.push(n->right);
			}
            q.pop();
        }
}

// the better way to implement this
void KDTree::nearest_(const point_t &pt, pointQueue &nearest_set, 
                        const KDNodePtr &node, const int topk){
    if(!bool(*node)) return;
    auto level = node->level;
    KDNodePtr sibling;
    if(pt.at(level) <= point_t(*node).at(level)){
        nearest_(pt, nearest_set, node->left, topk);
        sibling = node->right;
    } else {
        nearest_(pt, nearest_set, node->right, topk);
        sibling = node->left;
    }
    double d = dist2(point_t(*node), pt);
    auto pd = pointDist(node, d);
    queue_insert(nearest_set, pd, topk);
    int dx = point_t(*node).at(level) - pt.at(level);
    if(nearest_set.size() < topk || dx * dx < nearest_set.top().second) {
        nearest_(pt, nearest_set, sibling, topk);
    }
}

// the ugly implementation of knn
// konw the mind. Do not use it.
void KDTree::nearest_(const point_t &pt, pointQueue &nearest_set, const KDNodePtr &node, const int topk,  
                        const int level, std::set<size_t> visited){
    if(!bool(*node)) return;
    size_t dim = pt.size();
    auto left = node->left;
    auto right = node->right;
    double d, dx;
    if (!bool(*left) && bool(*right)) {
        nearest_(pt, nearest_set, right, topk, (level + 1) % dim, visited);

    } else if (bool(*left) && !bool(*right))
    {
        nearest_(pt, nearest_set, left, topk, (level + 1) % dim, visited);
    } else if(bool(*left) && bool(*right)) {
        if (point_t(*node).at(level) - pt.at(level) > 0) {
            nearest_(pt, nearest_set, left, topk, (level + 1) % dim, visited);
        } else {
            nearest_(pt, nearest_set, right, topk, (level + 1) % dim, visited);
        }
    } else {
        d = dist2(point_t(*node), pt);
        auto pd = pointDist(node, d);
        queue_insert(nearest_set, pd, topk);
        auto res_queue = nearest_set;
        while (!res_queue.empty()) 
        {
            res_queue.pop();
        } 

        visited.insert(node->index);
        auto pre = node;
        bool flag = false;
        while(bool(*pre) && root != pre && (visited.find(pre->index) != visited.end() || flag)){
            pre = pre->parent;
            bool left_visited = node_visited(pre->left, visited);
            bool right_visited = node_visited(pre->right, visited);
            flag = false;
            if(left_visited && right_visited && visited.find(pre->index) == visited.end()) {
                flag = true;
                d = dist2(point_t(*pre), pt);
                auto pd = pointDist(pre, d);
                queue_insert(nearest_set, pd, topk);
                visited.insert(pre->index);
            }
        }

        
        if(bool(*pre) && visited.find(pre->index) == visited.end()) {
            visited.insert(pre->index);
            bool left_visited = node_visited(pre->left, visited);
            bool right_visited = node_visited(pre->right, visited);
            d = dist2(point_t(*pre), pt);
            auto pd = pointDist(pre, d);
            queue_insert(nearest_set, pd, topk);
            auto l = pre->level;
            dx = point_t(*pre).at(l) - pt.at(l);
        auto res_queue = nearest_set;

            while (!res_queue.empty()) 
        {
            res_queue.pop();
        } 
            if (topk > nearest_set.size() || dx * dx < nearest_set.top().second) {
                if (!left_visited) {
                    nearest_(pt, nearest_set, pre->left, topk, (l + 1) % dim, visited);
                } else  {
                    nearest_(pt, nearest_set, pre->right, topk, (l + 1) % dim, visited);
                }
            }
        }    
    }
    
}

// for c++ using
// example in tests directory
void KDTree::nearest_topk_point(const point_t &pt, pointQueue &nearest_set, const int topk=1) {
    nearest_(pt,  nearest_set, root, topk);
}

// for python module
// nearest_topk_point_array(List, int)
pointVec KDTree::nearest_topk_point_array(const point_t &pt, const int topk=1){
    pointQueue nearest_set;
    nearest_topk_point(pt, nearest_set, topk);
    pointVec res;
    while(!nearest_set.empty()){
        res.push_back(point_t(*nearest_set.top().first));
        nearest_set.pop();
    }
    return res;   
}

// for python module
// nearest_topk_index(List, int)
std::vector<int> KDTree::nearest_topk_index(const point_t &pt, const int topk=1){
    pointQueue nearest_set;
    nearest_topk_point(pt, nearest_set, topk);
    std::vector<int> res;
    while(!nearest_set.empty()){
        res.push_back(nearest_set.top().first->index);
        nearest_set.pop();
    }
    return res;   
}

point_t KDTree::nearest_point(const point_t &pt) {
    return point_t(*nearest_(pt));
};
size_t KDTree::nearest_index(const point_t &pt) {
    return size_t(*nearest_(pt));
};


pointIndex KDTree::nearest_pointIndex(const point_t &pt) {
    KDNodePtr Nearest = nearest_(pt);
    return pointIndex(point_t(*Nearest), size_t(*Nearest));
}

pointIndexArr KDTree::neighborhood_(  //
    const KDNodePtr &branch,          //
    const point_t &pt,                //
    const double &rad,                //
    const size_t &level               //
) {
    double d, dx, dx2;

    if (!bool(*branch)) {
        // branch has no point, means it is a leaf,
        // no points to add
        return pointIndexArr();
    }

    size_t dim = pt.size();

    double r2 = rad * rad;

    d = dist2(point_t(*branch), pt);
    dx = point_t(*branch).at(level) - pt.at(level);
    dx2 = dx * dx;

    pointIndexArr nbh, nbh_s, nbh_o;
    if (d <= r2) {
        nbh.push_back(pointIndex(*branch));
    }

    //
    KDNodePtr section;
    KDNodePtr other;
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
    if (dx2 < r2) {
        nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
};

pointIndexArr KDTree::neighborhood(  //
    const point_t &pt,               //
    const double &rad) {
    size_t level = 0;
    return neighborhood_(root, pt, rad, level);
}

pointVec KDTree::neighborhood_points(  //
    const point_t &pt,                 //
    const double &rad) {
    size_t level = 0;
    pointIndexArr nbh = neighborhood_(root, pt, rad, level);
    pointVec nbhp;
    nbhp.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                   [](pointIndex x) { return x.first; });
    return nbhp;
}

indexArr KDTree::neighborhood_indices(  //
    const point_t &pt,                  //
    const double &rad) {
    size_t level = 0;
    pointIndexArr nbh = neighborhood_(root, pt, rad, level);
    indexArr nbhi;
    nbhi.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhi.begin(),
                   [](pointIndex x) { return x.second; });
    return nbhi;
}

namespace py = pybind11;

// for python module
// you can add new function here to expose to python
// atttion: some complex types are not supported
PYBIND11_MODULE(kdtree, m) {
    py::class_<KDTree>(m, "KDTree")
        .def(py::init<pointVec &>())
        .def("nearest_topk_point_array", &KDTree::nearest_topk_point_array)
        .def("nearest_topk_index", &KDTree::nearest_topk_index);
}