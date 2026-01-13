#ifndef HYBRID_VIZING
#define HYBRID_VIZING

#include <vector>

#include "DiGraphColor1.hpp"
#include "DiGraphColor2.hpp"

enum AlgorithmType {
    CLASSIC,
    GRAPH_1,
    GRAPH_2
};

class HybridVizing {
    public:
    HybridVizing(std::vector<std::vector<int>>& demand_matrix, const int T, const int ub, AlgorithmType alg_type);

    std::vector<std::vector<int>> color();

    void sequentialColor();
    bool colorGraph(bool only_case1);
    void colorClassic();

    void classicRotateFanCase1(int i, std::vector<std::pair<int,int>>& fan);
    void classicRotateFanCase2(int i, std::vector<std::pair<int,int>>& fan, int x, int y);

    void graphRotateFanCase1(int i, std::vector<std::pair<int,int>>& fan, FanGraph2* G, FreeGraph* H);
    int graphFindFanCase2(int u, int v, FanGraph2* G);
    void swapAugmentingPath(int alpha, int beta, int i, FanGraph2* G, FreeGraph* H);

    std::vector<std::vector<int>>& demand_matrix;
    const int T;
    int colors;

    AlgorithmType alg_type;

    std::vector<std::vector<int>> edges_per_color;

    std::vector<std::pair<int,int>> case1_pairs;
    std::vector<std::pair<int,int>> case2_pairs;

    int num_failures = 0;

    #ifdef ANALYZE_VIZING
    std::vector<int> cases = std::vector<int>(3);
    std::vector<int> fan_lengths;
    #endif
};

#endif