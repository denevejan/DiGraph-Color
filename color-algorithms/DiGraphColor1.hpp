#ifndef DIGRAPH_COLOR_1
#define DIGRAPH_COLOR_1

#include "macros.hpp"

#include "help_structures.hpp"

#include <vector>
#include <queue>

class DiGraphColor1 {
    public:
    DiGraphColor1(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub);

    std::vector<std::vector<int>> color();

    void rotateFanCase1(int i, std::vector<std::pair<int,int>>& fan, FanGraph2& G, FreeGraph& H);

    const std::vector<std::vector<int>>& demand_matrix;
    const int T;
    int colors;

    int num_failures = 0;

    #ifdef ANALYZE_VIZING
    std::vector<int> cases = std::vector<int>(3);
    std::vector<int> fan_lengths;
    #endif

    std::vector<std::vector<int>> edges_per_color;
};

#endif