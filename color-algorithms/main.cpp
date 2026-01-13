#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdio>

#include <algorithm>

#include "InitializeTestInput.hpp"
#include "CheckSolutions.hpp"

#include "Vizing.hpp"
#include "DiGraphColor1.hpp"
#include "DiGraphColor2.hpp"
#include "HybridVizing.hpp"
#include "VizingCommonC.hpp"
#include "Karloff.hpp"

void testEverythingRandom(const std::vector<std::vector<int>>& demand_matrix, const int& T, const int& max_degree, const int ub, std::vector<int>& runtimes, std::vector<int>& num_colors, std::vector<std::vector<double>>& cases_distribution, std::vector<std::vector<int>>& fan_lengths);

int getPairIndex(int T, int S);
void writeHeadersLine(std::ofstream& output);
void writeLineToFile(std::ofstream& output, int T, int S, std::vector<std::vector<int>>& all_runtimes, std::vector<std::vector<int>>& all_num_colors, double upper_bound, std::vector<std::vector<std::vector<double>>>& all_cases_distribution, std::vector<std::vector<int>>& all_fan_lengths, int current_exp_ctr);


using namespace std;

#define NUM_EXPERIMENTS 10
#define NUM_ALGORITHMS 8

// note that we use "T" as the number of nodes in our code, while the more common "n" is used in the paper
// this is because in the research about the multi-GPU system, a node corresponds to a gpu Tile
#define T_START 4
#define T_END 32
// p is a helper variable to easily obtain the logarithmic-type (instead of linear/uniform) distribution of the input variable S
#define P_START 6
#define P_END 12
#define SAMPLE_RATE_LG 3 // lg of the number of different input sizes chosen between subsequent powers of two
#define INTERMEDIATE_OUTPUT_RATE 1

int main(int argc, char** argv) {
    vector<pair<int,int>> T_S_pairs;
    for (int T = T_START; T <= T_END; T*=2) {
        for (int p = P_START; p < P_END; p++) {
            for (int S = 1 << p; S < (1 << (p+1)); S += 1 << (p-SAMPLE_RATE_LG)) {
                T_S_pairs.push_back(pair<int,int>(T,S));
            }
        }
    }
    int num_pairs = T_S_pairs.size();

    InitializeTestInput iti;

    vector<vector<vector<int>>> all_runtimes(num_pairs, vector<vector<int>>(NUM_EXPERIMENTS, vector<int>(NUM_ALGORITHMS)));
    vector<vector<vector<int>>> all_num_colors(num_pairs, vector<vector<int>>(NUM_EXPERIMENTS, vector<int>(NUM_ALGORITHMS)));
    vector<vector<vector<vector<double>>>> all_cases_distribution(num_pairs, vector<vector<vector<double>>>(NUM_EXPERIMENTS, vector<vector<double>>(NUM_ALGORITHMS)));
    vector<vector<vector<int>>> all_fan_lengths(num_pairs, vector<vector<int>>(NUM_ALGORITHMS, vector<int>()));

    vector<double> sum_upper_bound(num_pairs);
    
    for (int experiment_counter = 0; experiment_counter < NUM_EXPERIMENTS; ) {
        // we include randomness in the order that we test each combination of T and S, to spread out as much as possible any runtime irregularities stemming from the testing environment (and are not inherent to the alorithms themselves)
        random_shuffle(T_S_pairs.begin(), T_S_pairs.end());

        for (int pair_ctr = 0; pair_ctr < num_pairs; pair_ctr++) {
            pair<int,int> T_S_pair = T_S_pairs[pair_ctr];
            int T = T_S_pair.first;
            int S = T_S_pair.second;
            int pair_index = getPairIndex(T,S);

            vector<vector<int>> demand_matrix(T, vector<int>(T));
            int ub = iti.randomERGraph(demand_matrix, T, S*T/2);
            sum_upper_bound[pair_index] += ub;
            int max_degree = iti.maxDegree(demand_matrix, T);
            testEverythingRandom(demand_matrix, T, max_degree, ub, all_runtimes[pair_index][experiment_counter], all_num_colors[pair_index][experiment_counter], all_cases_distribution[pair_index][experiment_counter], all_fan_lengths[pair_index]);
        }

        experiment_counter++;

        if (experiment_counter % INTERMEDIATE_OUTPUT_RATE == 0) {
            ofstream output; output.open("X" + to_string(experiment_counter) + ".csv");
            writeHeadersLine(output);
            for (int T = T_START, pair_index = 0; T <= T_END; T *=2) {
                for (int p = P_START; p < P_END; p++) {
                    for (int S = 1 << p; S < (1 << (p+1)); S += 1 << (p-SAMPLE_RATE_LG), pair_index++) {
                        writeLineToFile(output, T, S, all_runtimes[pair_index], all_num_colors[pair_index], sum_upper_bound[pair_index] / experiment_counter, all_cases_distribution[pair_index], all_fan_lengths[pair_index], experiment_counter);
                    }
                }
            }
            output.close();
        }
        if (experiment_counter % 2 == 0) cout << experiment_counter << " experiments finished" << endl;
    }
    return 0;
}


void testEverythingRandom(const vector<vector<int>>& demand_matrix, const int& T, const int& max_degree, const int ub, vector<int>& runtimes, vector<int>& num_colors, vector<vector<double>>& cases_distribution, vector<vector<int>>& fan_lengths) {
    // for a given demand matrix, the different algorithms are also tested in a random order to avoid any runtime phenomena that are inherent to the testing order (e.g. warm vs cold caches) and not to the algorithms themselves
    vector<int> random_order(NUM_ALGORITHMS);
    for (int i = 0; i < NUM_ALGORITHMS; i++) random_order[i] = i;
    random_shuffle(random_order.begin(), random_order.end());

    chrono::_V2::system_clock::time_point start_time;
    chrono::_V2::system_clock::time_point end_time;

    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        vector<vector<int>> demand_matrix_copy(demand_matrix);
        vector<vector<int>> solution;
        int num_failures = 0;
        switch (random_order[i]) {
            case 0: {
                GeneralVizing vv(demand_matrix_copy, T, ub);
                start_time = chrono::high_resolution_clock::now();
                solution = vv.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : vv.fan_lengths) fan_lengths[0].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[0].push_back(((double)vv.cases[j]) / (vv.cases[0]+vv.cases[1]+vv.cases[2]));
                break;
            }
            case 1: {
                DiGraphColor1 dgc1(demand_matrix_copy, T, ub);
                start_time = chrono::high_resolution_clock::now();
                solution = dgc1.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : dgc1.fan_lengths) fan_lengths[1].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[1].push_back(((double)dgc1.cases[j]) / (dgc1.cases[0]+dgc1.cases[1]+dgc1.cases[2]));
                num_failures = dgc1.num_failures;
                break;
            }
            case 2: {
                DiGraphColor2 dgc2(demand_matrix_copy, T, ub);
                start_time = chrono::high_resolution_clock::now();
                solution = dgc2.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : dgc2.fan_lengths) fan_lengths[2].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[2].push_back(((double)dgc2.cases[j]) / (dgc2.cases[0]+dgc2.cases[1]+dgc2.cases[2]));
                break;
            }
            case 3: {
                HybridVizing hv(demand_matrix_copy, T, ub, AlgorithmType::CLASSIC);
                start_time = chrono::high_resolution_clock::now();
                solution = hv.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : hv.fan_lengths) fan_lengths[3].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[3].push_back(((double)hv.cases[j]) / (hv.cases[0]+hv.cases[1]+hv.cases[2]));
                break;
            }
            case 4: {
                HybridVizing hgv(demand_matrix_copy, T, ub, AlgorithmType::GRAPH_1);
                start_time = chrono::high_resolution_clock::now();
                solution = hgv.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : hgv.fan_lengths) fan_lengths[4].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[4].push_back(((double)hgv.cases[j]) / (hgv.cases[0]+hgv.cases[1]+hgv.cases[2]));
                num_failures = hgv.num_failures;
                break;
            }
            case 5: {
                HybridVizing hgv2(demand_matrix_copy, T, ub, AlgorithmType::GRAPH_2);
                start_time = chrono::high_resolution_clock::now();
                solution = hgv2.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : hgv2.fan_lengths) fan_lengths[5].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[5].push_back(((double)hgv2.cases[j]) / (hgv2.cases[0]+hgv2.cases[1]+hgv2.cases[2]));
                break;
            }
            case 6: {
                VizingCommonC vcc(demand_matrix_copy, T, ub);
                start_time = chrono::high_resolution_clock::now();
                solution = vcc.color();
                end_time = chrono::high_resolution_clock::now();
                for (int l : vcc.fan_lengths) fan_lengths[6].push_back(l);
                for (int j = 0; j < 3; j++) cases_distribution[6].push_back(((double)vcc.cases[j]) / (vcc.cases[0]+vcc.cases[1]+vcc.cases[2]));
                break;
            }
            case 7: {
                Karloff k(demand_matrix_copy, T, max_degree);
                start_time = chrono::high_resolution_clock::now();
                solution = k.color();
                end_time = chrono::high_resolution_clock::now();
                fan_lengths[7].push_back(1);
                for (int j = 0; j < 3; j++) cases_distribution[7].push_back(0);
                break;
            }
        }
        runtimes[random_order[i]] = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
        num_colors[random_order[i]] = solution.size();

        #ifdef CHECK_SOLUTIONS
        checkSolution(demand_matrix, T, solution, num_failures);
        #endif
    }
}


int getPairIndex(int T, int S) {
    int p_block = 8*sizeof(int) - __builtin_clz(S) - P_START - 1;
    return (__builtin_ctz(T / T_START) * ((P_END - P_START) << SAMPLE_RATE_LG)) + (p_block << SAMPLE_RATE_LG) + ((S ^ (1 << (p_block + P_START))) >> (p_block + (P_START-SAMPLE_RATE_LG)));
}


void writeHeadersLine(std::ofstream& output) {
    // HEADERS OF THE OUTPUT FILE
    output << "#exp," << "S," << "T,";
    output
    << "t_vizing," << "rmse_vizing,"
    << "t_digraphcolor1," << "rmse_digraphcolor1,"
    << "t_digraphcolor2," << "rmse_digraphcolor2,"
    << "t_hybridvizing," << "rmse_hybridvizing,"
    << "t_hybriddigraphcolor1," << "rmse_hybriddigraphcolor1,"
    << "t_hybriddigraphcolor2," << "rmse_hybriddigraphcolor2,"
    << "t_vizingcc," << "rmse_vizingcc,"
    << "t_shannon," << "rmse_shannon,"
    ;
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        output << "ave c_" << i << ",rmse c_" << i << ",";
    }
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        output << "case1_" << i << ",case2_" << i << ",fansize_" << i << ",";
    }
    output << "ave upper bound" << endl;
}

void writeLineToFile(ofstream& output, int T, int S, vector<vector<int>>& all_runtimes, vector<vector<int>>& all_num_colors, double upper_bound, vector<vector<vector<double>>>& all_cases_distribution, vector<vector<int>>& all_fan_lengths, int current_exp_ctr) {
    // CALCULATE THE AVERAGE RUNTIME AND NUMBER OF USED COLORS OF EACH ALGORITHM
    vector<float> ave_runtimes(NUM_ALGORITHMS);
    vector<float> ave_num_colors(NUM_ALGORITHMS);
    vector<vector<float>> ave_cases_distribution(NUM_ALGORITHMS, vector<float>(3));
    vector<double> ave_fan_lengths(NUM_ALGORITHMS);
    vector<int> fanless_experiments(NUM_ALGORITHMS);
    for (int j = 0; j < current_exp_ctr; j++) {
        for (int i = 0; i < NUM_ALGORITHMS; i++) {
            ave_runtimes[i] += all_runtimes[j][i];
            ave_num_colors[i] += all_num_colors[j][i];
            for (int k = 0; k < 3; k++) {
                ave_cases_distribution[i][k] += all_cases_distribution[j][i][k];
            }
        }
    }
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        ave_runtimes[i] /= current_exp_ctr;
        ave_num_colors[i] /= current_exp_ctr;
        for (int k = 0; k < 3; k++) {
            ave_cases_distribution[i][k] /= current_exp_ctr;
        }
        for (int l : all_fan_lengths[i]) ave_fan_lengths[i] += l;
        ave_fan_lengths[i] /= all_fan_lengths[i].size();
    }

    // CALCULATE THE RMSE OF THE RUNTIMES AND NUMBER OF USED COLORS FOR EACH ALGORITHM
    vector<float> rmse_runtimes(NUM_ALGORITHMS);
    vector<float> rmse_num_colors(NUM_ALGORITHMS);
    for(int j = 0; j < current_exp_ctr; j++){
        for (int i = 0; i < NUM_ALGORITHMS; i++) {
            if (all_num_colors[j][i] == 0) continue;
            rmse_runtimes[i] += pow((all_runtimes[j][i]-ave_runtimes[i]),2);
            rmse_num_colors[i] += pow((all_num_colors[j][i]-ave_num_colors[i]),2);
        }
    }
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        rmse_runtimes[i] = sqrt(rmse_runtimes[i] / current_exp_ctr);
        rmse_num_colors[i] = sqrt(rmse_num_colors[i] / current_exp_ctr);
    }

    output << current_exp_ctr << "," << S << "," << T << ",";
    // OUTPUT THE RUNTIMES
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        output << ave_runtimes[i] << "," << rmse_runtimes[i] << ",";
    }
    // OUTPUT THE NUMBER OF USED COLORS
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        output << ave_num_colors[i] << "," << rmse_num_colors[i] << ",";
    }
    // OUTPUT THE VIZING ANALYSIS OF CASES (will be zero if analysis is skipped)
    for (int i = 0; i < NUM_ALGORITHMS; i++) {
        output << ave_cases_distribution[i][1] << "," << ave_cases_distribution[i][2] << "," << ave_fan_lengths[i] << ",";
    }
    // OUTPUT THE AVERAGE CHROMATIC NUMBER
    output << upper_bound << endl;
}
