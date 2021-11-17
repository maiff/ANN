#include <iostream>
#include <vector>

#include "KDTree.hpp"

// using point_t = std::vector< double >;
// using pointVec = std::vector< point_t >;

// point_t pt(2);

int main() {
    pointVec points;
    point_t pt;

    // pt = {0.0, 0.0}; // 0.68
    // points.push_back(pt);
    // pt = {1.0, 0.0}; // 0.08
    // points.push_back(pt);
    // pt = {0.0, 1.0}; // 1.28
    // points.push_back(pt);
    // pt = {1.0, 1.0}; // 0.68
    // points.push_back(pt);
    // pt = {0.5, 0.5}; // 0.18
    // points.push_back(pt);
    // pt = {0.8, 0.2}; // 0.0
    // points.push_back(pt);

    // KDTree tree(points);

    // std::cout << "nearest test\n";
    // pt = {0.8, 0.2};
    // auto res = tree.nearest_point(pt);
    // auto res_index = tree.nearest_index(pt);
    // pointQueue res_queue;
    // tree.nearest_topk_point(pt, res_queue, 2);
    // while (!res_queue.empty()) 
    // {
    //     std::cout << (res_queue.top().first->index) << "//////////////////" << std::endl;
    //     res_queue.pop();
    // } 
    // for (double b : res) {
    //     std::cout << b << " ";
    // }
    // std::cout << '\n';
    // std::cout << "res_index: " << res_index << std::endl;

    pt = {1.0, 1.0}; // 0.68
    points.push_back(pt);
    pt = {2.0, 2.0}; // 0.08
    points.push_back(pt);
    pt = {3.0, 3.0}; // 1.28
    points.push_back(pt);
    // pt = {2.0, 2.0}; // 0.68
    // points.push_back(pt);
    // pt = {3.0, 1.0}; // 0.18
    // points.push_back(pt);
    // pt = {3.0, 2.0}; // 0.0
    // points.push_back(pt);
    // pt = {3.0, 3.0}; // 0.0
    // points.push_back(pt);

    KDTree tree(points);
    tree.view(tree.root);
    pt = {2.0, 2.0};
    auto res = tree.nearest_point(pt);
    auto res_index = tree.nearest_index(pt);
    pointQueue res_queue;
    tree.nearest_topk_point(pt, res_queue, 3);
    while (!res_queue.empty()) 
    {
        std::cout << (res_queue.top().first->index) << "//////////////////" << std::endl;
        res_queue.pop();
    } 
    for (double b : res) {
        std::cout << b << " ";
    }
    std::cout << '\n';
    std::cout << "res_index: " << res_index << std::endl;
    

    /*
    std::cout << "going down the tree\n";

    for (auto b : point_t(*tree)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->right)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->left->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->right->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    std::cout << "printing nbh\n";

    pt = {.0, .5};

    */
    auto res2 = tree.neighborhood_points(pt, .55);

    for (point_t a : res2) {
        for (double b : a) {
            std::cout << b << " ";
        }
        std::cout << '\n';
    }
    return 0;
}
