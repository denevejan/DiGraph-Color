#include "VizingCommonC.hpp"

#include <iostream>

using namespace std;

VizingCommonC::VizingCommonC(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub) : demand_matrix(demand_matrix), T(T), colors(ub) {
    // initialize all colors to be free at the start
    edges_per_color = vector<vector<int>>(colors, vector<int>(T, -1));
};

vector<vector<int>> VizingCommonC::color() {
    for (int i = 0; i < T-1; i++) {
        for (int j = i+1; j < T; j++) {
            int d = demand_matrix[i][j];

            // first find all (or as many as needed) common free colors between i and j before going to the fan building 
            for (int c_common = 0; c_common < colors && d > 0; c_common++) {
                if (edges_per_color[c_common][i] == -1 && edges_per_color[c_common][j] == -1) {
                    edges_per_color[c_common][i] = j;
                    edges_per_color[c_common][j] = i;
                    d--;
                    #ifdef ANALYZE_VIZING
                    cases[0]++;
                    #endif
                }
            }
            
            int c_i = -1; int c_j = -1;
            for (; d --> 0 ;) {
                int free_c_i = -1;
                for (c_i++; c_i < colors; c_i++) {
                    if (edges_per_color[c_i][i] == -1) {
                        free_c_i = c_i;
                        break;
                    }
                }
                int free_c_j = -1;
                for (c_j++; true; c_j++) {
                    if (c_j >= colors) c_j = 0; // search circular (more efficient than starting from zero each time)
                    if (edges_per_color[c_j][j] == -1) {
                        free_c_j = c_j;
                        break;
                    }
                }

                vector<pair<int,int>> fan;
                vector<bool> fan_colors_array(colors, false);
                vector<int> color_position_per_j(T, -1);
                fan.push_back(pair<int,int>(j, free_c_j));
                int next_j = j;
                while (true) {
                    fan_colors_array[free_c_j] = true;

                    color_position_per_j[next_j] = free_c_j;
                    next_j = edges_per_color[free_c_j][i];
                    // this is only for the first free_c_j, all others start from zero so their "circular" search is just a normal search from zero
                    for (free_c_j = color_position_per_j[next_j]+1; true; free_c_j++) {
                        if (free_c_j == colors) free_c_j = 0;
                        if (edges_per_color[free_c_j][next_j] < 0) break;
                    }

                    fan.push_back(pair<int,int>(next_j, free_c_j));

                    if (edges_per_color[free_c_j][i] == -1) {
                        rotateFan_case1(i, fan);
                        #ifdef ANALYZE_VIZING
                        cases[1]++;
                        fan_lengths.push_back(fan.size());
                        #endif
                        break;
                    }
                    if (fan_colors_array[free_c_j]) {
                        rotateFan_case2(i, fan, free_c_i, free_c_j);
                        #ifdef ANALYZE_VIZING
                        cases[2]++;
                        fan_lengths.push_back(fan.size());
                        #endif
                        break;
                    }
                }
            }
        }
    }
    return edges_per_color;
}


void VizingCommonC::rotateFan_case1(int i, vector<pair<int,int>>& fan) {
    int last_c = -1;
    for (pair<int,int> edge_hole_pair : fan) {
        int j = edge_hole_pair.first;
        int c = edge_hole_pair.second;

        if (last_c >= 0) {
            edges_per_color[last_c][j] = -1; // uncolor the edge incident to i that previously had this color 
        }

        edges_per_color[c][i] = j;
        edges_per_color[c][j] = i;

        last_c = c;
    }
}

void VizingCommonC::rotateFan_case2(int i, vector<pair<int,int>>& fan, int x, int y) {
    // alternate xy-path
    int ii = i;
    int jj = edges_per_color[y][ii];
    bool odd_edge = false;
    while (jj != -1) {
        swap(edges_per_color[x][ii], edges_per_color[y][ii]);
        ii = jj;
        jj = edges_per_color[odd_edge?y:x][jj];
        odd_edge = !odd_edge;
    }
    swap(edges_per_color[x][ii], edges_per_color[y][ii]);

    bool case_2b = false;
    int last_c = -1;
    for (pair<int,int> edge_hole_pair : fan) {
        if (case_2b) break;
        int j = edge_hole_pair.first;
        int c = edge_hole_pair.second;

        if (c == y) { // last blade of the fan before the blade that starts the augmenting path
            if (ii == j) { // case 2a
                c = x;
            } else { // case 2b
                case_2b = true;
            }
        }

        edges_per_color[c][i] = j;
        edges_per_color[c][j] = i;

        if (last_c >= 0) {
            edges_per_color[last_c][j] = -1; // uncolor the edge incident to i that previously had this color 
        }

        last_c = c;
    }
}
