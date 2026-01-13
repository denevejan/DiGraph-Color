#ifndef CHECK_SOLUTIONS_HEADER_FILE
#define CHECK_SOLUTIONS_HEADER_FILE

#define CHECK_SOLUTIONS // comment out this line if you do not want the overhead of checking the validity of solutions

#ifdef CHECK_SOLUTIONS

#include <cassert>
#include <iostream>

#include <vector>

void checkSolution(const std::vector<std::vector<int>>& demand_matrix, const int& T, std::vector<std::vector<int>> solution, int num_failures=0, bool with_excess=false);

int channelsDesired(const std::vector<std::vector<int>>& demand_matrix, int T);

void checkTotalCases(const std::vector<std::vector<int>>& demand_matrix, const int& T, std::vector<std::vector<int>> cases, int num_alg, std::vector<int> solution_sizes, std::vector<int> num_failures, bool with_excess=false);

#endif

#endif