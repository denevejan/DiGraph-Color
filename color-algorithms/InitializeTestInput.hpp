#ifndef INITIALIZE_TEST_INPUT
#define INITIALIZE_TEST_INPUT

#include <cstdlib>
#include <cmath>

#include <vector>
#include <map>

// class for initializing the demand matrix for testing the coloring algorithms

class InitializeTestInput {
    public:
    InitializeTestInput();

    int initializeApplicationInput(std::vector<std::vector<int>>& demand_matrix, int T, int S, int num_peaks=0);
    void phase_1(std::vector<std::vector<int>>& d, unsigned int N, unsigned int T);
    void phase_2(std::vector<std::vector<int>>& orig_d, std::vector<std::vector<int>>& a, float alpha,unsigned int N,unsigned int T);

    int randomERGraph(std::vector<std::vector<int>>& demand_matrix, int T, int S);

    int maxDegree(std::vector<std::vector<int>>& demand_matrix, int T);
    int maxDegree(std::vector<std::vector<int>>& demand_matrix, int T, int& max_i);
    int upperBound(std::vector<std::vector<int>>& demand_matrix, int T);

    static void showMatrix(const std::vector<std::vector<int>>& demand_matrix, int T);

    int initializations_done = 0;
};


#endif