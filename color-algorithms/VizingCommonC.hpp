#ifndef VIZING_COMMON_C
#define VIZING_COMMON_C

#include "macros.hpp"

#include <vector>

// the majority of code in this class is identical to GeneralVizing, as VizingCommonC is the same algorithm extended with the common color heuristic
class VizingCommonC {
    public:
    VizingCommonC(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub);

    std::vector<std::vector<int>> color();

    void rotateFan_case1(int i, std::vector<std::pair<int,int>>& fan);
    void rotateFan_case2(int i, std::vector<std::pair<int,int>>& fan, int x, int y);
    
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
