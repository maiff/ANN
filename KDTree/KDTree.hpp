#pragma once

/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
#include <queue>


using point_t = std::vector< double >;
using indexArr = std::vector< size_t >;
using pointIndex = typename std::pair< std::vector< double >, size_t >;

class KDNode {
   public:
    using KDNodePtr = std::shared_ptr< KDNode >;
    size_t index;
    point_t x;
    KDNodePtr left;
    KDNodePtr right;
    KDNodePtr parent;
    size_t level;

    // initializer
    KDNode();
    KDNode(const point_t &, const size_t &, const KDNodePtr &,
           const KDNodePtr &);
    KDNode(const pointIndex &, const KDNodePtr &, const KDNodePtr &);
    ~KDNode();

    // getter
    double coord(const size_t &);

    // conversions
    explicit operator bool();
    explicit operator point_t();
    explicit operator size_t();
    explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr< KDNode >;

KDNodePtr NewKDNodePtr();

using pointDist = typename std::pair< KDNodePtr, double >;
struct cmp_by_second {
  bool operator()(pointDist const &a, pointDist const &b) {
    return a.second < b.second;
  }
};

using pointQueue = std::priority_queue<pointDist, std::vector<pointDist>, cmp_by_second>;
void queue_insert(pointQueue &queue, pointDist &item, int k);

// square euclidean distance
inline double dist2(const point_t &, const point_t &);
inline double dist2(const KDNodePtr &, const KDNodePtr &);

// euclidean distance
inline double dist(const point_t &, const point_t &);
inline double dist(const KDNodePtr &, const KDNodePtr &);

// Need for sorting
class comparer {
   public:
    size_t idx;
    explicit comparer(size_t idx_);
    inline bool compare_idx(
        const std::pair< std::vector< double >, size_t > &,  //
        const std::pair< std::vector< double >, size_t > &   //
    );
};

using pointIndexArr = typename std::vector< pointIndex >;

inline void sort_on_idx(const pointIndexArr::iterator &,  //
                        const pointIndexArr::iterator &,  //
                        size_t idx);

using pointVec = std::vector< point_t >;

class KDTree {
    KDNodePtr leaf;

    KDNodePtr make_tree(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        const size_t &length,                  //
                        const size_t &level                    //
    );

   public:
    KDTree() = default;
    KDNodePtr root;
    explicit KDTree(pointVec point_array);
    void view(KDNodePtr node);

   private:
    KDNodePtr nearest_(           //
        const KDNodePtr &branch,  //
        const point_t &pt,        //
        const size_t &level,      //
        const KDNodePtr &best,    //
        const double &best_dist,   //
        pointQueue &nearest_set,
        const int topk
    );

    void nearest_(const point_t &pt, pointQueue &nearest_set, const KDNodePtr &node, const int topk,  
                        const int level, std::set<size_t> visited);

    void nearest_(const point_t &pt, pointQueue &nearest_set, 
                        const KDNodePtr &node, const int topk);
    // default caller
    KDNodePtr nearest_(const point_t &pt, const int topk);

   public:
    point_t nearest_point(const point_t &pt);
    size_t nearest_index(const point_t &pt);
    pointIndex nearest_pointIndex(const point_t &pt);
    void nearest_topk_point(const point_t &pt, pointQueue &nearest_set, const int topk) ;
    pointVec nearest_topk_point_array(const point_t &pt, const int topk);
    std::vector<int> nearest_topk_index(const point_t &pt, const int topk);

   private:
    pointIndexArr neighborhood_(  //
        const KDNodePtr &branch,  //
        const point_t &pt,        //
        const double &rad,        //
        const size_t &level       //
    );

   public:
    pointIndexArr neighborhood(  //
        const point_t &pt,       //
        const double &rad);

    pointVec neighborhood_points(  //
        const point_t &pt,         //
        const double &rad);

    indexArr neighborhood_indices(  //
        const point_t &pt,          //
        const double &rad);
};
