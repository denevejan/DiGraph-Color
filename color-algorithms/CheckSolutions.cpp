#include "CheckSolutions.hpp"
#include "InitializeTestInput.hpp"

using namespace std;

#ifdef CHECK_SOLUTIONS

void checkSolution(const std::vector<std::vector<int>>& demand_matrix, const int& T, std::vector<std::vector<int>> solution, int num_failures, bool with_excess) {
    std::vector<std::vector<int>> color_checker(T, std::vector<int>(T));
    for (int c = 0; c < solution.size(); c++) {
        bool color_used = false;
        for (int i = 0; i < T; i++) {
            int j = solution[c][i];
            if (j == -1) {
                continue;
            }
            color_used = true;
            color_checker[i][j]++;

            // CHECK THAT NO TWO EDGES OF THE SAME COLOR ARRIVE IN THE SAME NODE (RIGHT-HAND SIDE OF THE BIPARTITE GRAPH)
            for (int ii = 0; ii < T; ii++) {
                if (ii == i) {
                    continue;
                }
                assert(j != solution[c][ii]);
            }
        }
    }

    // CHECK THAT EACH MULTI-EDGE IS ASSIGNED THE CORRECT NUMBER OF COLORS
    for (int i = 0; i < T; i++) {
        for (int j = i+1; j < T; j++) {
            if (with_excess) {
                assert(demand_matrix[i][j] <= color_checker[i][j]);
            } else {
                assert(demand_matrix[i][j] >= color_checker[i][j]);
                // if (demand_matrix[i][j] > color_checker[i][j]) {
                //     cout << "shortage of " << demand_matrix[i][j] - color_checker[i][j] << " in " << i << "," << j << endl;
                // }
                num_failures -= demand_matrix[i][j] - color_checker[i][j];
            }
        }
    }
    assert(num_failures == 0);
}

int channelsDesired(const std::vector<std::vector<int>>& demand_matrix, int T) {
    int desired = 0;
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++) {
            int demand = demand_matrix[i][j];
            desired += demand;
        }
    }
    return desired;
}

void checkTotalCases(const vector<vector<int>>& demand_matrix, const int& T, vector<vector<int>> cases, int num_alg, vector<int> solution_sizes, vector<int> num_failures, bool with_excess) {
    int matrix_sum = 0;
    for (int i = 0; i < T; i++) {
        for (int j = i+1; j < T; j++) {
            matrix_sum += demand_matrix[i][j];
        }
    }
    for (int i = 0; i < num_alg; i++) {
        int sum_cases = 0;
        for (int j = 0; j < 3; j++) {
            sum_cases += cases[i][j];
        }
        if (with_excess) {
            assert(matrix_sum <= sum_cases+num_failures[i]);
        } else {
            assert(matrix_sum == sum_cases+num_failures[i]);
        }
    }
}

#endif
