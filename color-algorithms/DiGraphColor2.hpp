#ifndef DiGRAPH_COLOR_2
#define DiGRAPH_COLOR_2

#include "macros.hpp"

#include "help_structures.hpp"

#include <vector>
#include <queue>
#include <tuple>

class DiGraphColor2 {
    public:
    DiGraphColor2(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub);

    std::vector<std::vector<int>> color();

    void rotateFanCase1(int i, std::vector<std::pair<int,int>>& fan, FanGraph2& G, FreeGraph& H);

    int findFanCase2(int u, int v, FanGraph2& G);
    void swapAugmentingPath(int alpha, int beta, int i, FanGraph2& G, FreeGraph& H);

    const std::vector<std::vector<int>>& demand_matrix;
    const int T;
    int colors;

    #ifdef ANALYZE_VIZING
    std::vector<int> cases = std::vector<int>(3);
    std::vector<int> fan_lengths;
    #endif

    std::vector<std::vector<int>> edges_per_color;
};

#endif