#include "HybridVizing.hpp"

#include "macros.hpp"

#include <iostream>

using namespace std;

HybridVizing::HybridVizing(vector<vector<int>>& demand_matrix, const int T, const int ub, AlgorithmType alg_type) : demand_matrix(demand_matrix), T(T), colors(ub), alg_type(alg_type) {
    edges_per_color = vector<vector<int>>(colors, vector<int>(T, -1));
}

vector<vector<int>> HybridVizing::color() {
    sequentialColor();
    if (alg_type == CLASSIC) {
        colorClassic();
    } else {
        if(!colorGraph(true)) {
            return vector<vector<int>>();
        }
        if (alg_type == GRAPH_2) {
            colorGraph(false);
        }
    }
    return edges_per_color;
}

void HybridVizing::sequentialColor() {
    for (int i = 0; i < T-1; i++) {
        for (int j = i+1; j < T; j++) {
            int d = demand_matrix[i][j];
            #ifdef ANALYZE_VIZING
            cases[0] += d;
            #endif
            int c = -1;
            for (; d > 0; d--) {
                for (c++; c < colors; c++) {
                    if (edges_per_color[c][i] == -1 && edges_per_color[c][j] == -1) {
                        edges_per_color[c][i] = j;
                        edges_per_color[c][j] = i;
                        break;
                    }
                }
                if (c >= colors) break;
            }
            demand_matrix[i][j] = demand_matrix[j][i] = d;
            if (d > 0) {
                case1_pairs.push_back(make_pair(i,j));
            }
            #ifdef ANALYZE_VIZING
            cases[0] -= d;
            #endif
        }
    }
}

void HybridVizing::colorClassic() {
    for (pair<int,int> p : case1_pairs) {
        int i = p.first;
        int j = p.second;

        int d = demand_matrix[i][j];
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
                if (c_j >= colors) {
                    c_j = 0; // search circular (more efficient than starting from zero each time)
                }
                if (edges_per_color[c_j][j] == -1) {
                    free_c_j = c_j;
                    break;
                }
            }
            if (edges_per_color[free_c_i][j] == -1) { // common free color
                edges_per_color[free_c_i][j] = i;
                edges_per_color[free_c_i][i] = j;
                #ifdef ANALYZE_VIZING
                cases[0]++;
                #endif
                continue;
            }
            if (edges_per_color[free_c_j][i] == -1) { // common free color
                edges_per_color[free_c_j][j] = i;
                edges_per_color[free_c_j][i] = j;
                #ifdef ANALYZE_VIZING
                cases[0]++;
                #endif
                continue;
            }
            // directly assigning a common free color is strictly speaking not part of Vizing's Theorem. But we do it to save trouble with the code later on
            
            vector<pair<int,int>> fan;
            vector<bool> fan_colors_array(colors, false);
            vector<int> color_position_per_j(T, -1);
            fan.push_back(pair<int,int>(j, free_c_j));
            int next_j = j;
            while (true) {
                fan_colors_array[free_c_j] = true;

                color_position_per_j[next_j] = free_c_j;
                next_j = edges_per_color[free_c_j][i];
                for (free_c_j = color_position_per_j[next_j]+1; true; free_c_j++) {
                    if (free_c_j >= colors) free_c_j = 0;
                    if (edges_per_color[free_c_j][next_j] < 0) break;
                }

                fan.push_back(pair<int,int>(next_j, free_c_j));

                if (edges_per_color[free_c_j][i] == -1) {
                    classicRotateFanCase1(i, fan);
                    #ifdef ANALYZE_VIZING
                    cases[1]++;
                    fan_lengths.push_back(fan.size());
                    #endif
                    break;
                }
                if (fan_colors_array[free_c_j]) {
                    classicRotateFanCase2(i, fan, free_c_i, free_c_j);
                    #ifdef ANALYZE_VIZING
                    cases[2]++;
                    fan_lengths.push_back(fan.size());
                    #endif
                    break;
                }
            }
        }
        demand_matrix[i][j] = demand_matrix[j][i] = 0;
    }
}

void HybridVizing::classicRotateFanCase1(int i, vector<pair<int,int>>& fan) {
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

void HybridVizing::classicRotateFanCase2(int i, vector<pair<int,int>>& fan, int x, int y) {
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

bool HybridVizing::colorGraph(bool only_case1) {
    FreeGraph* H;
    FanGraph2* G;
    vector<pair<int,int>>& caseX_pairs = only_case1 ? case1_pairs : case2_pairs;
    int last_i = -1;
    for (pair<int,int> p : caseX_pairs) {
        int i = p.first;
        if (i != last_i) {
            if (last_i != -1) {
                delete H;
                delete G;
            }
            H = new FreeGraph(edges_per_color, T, colors, i);
            G = new FanGraph2(edges_per_color, T, colors, i, only_case1);
            last_i = i;
        }
        int j = p.second;
        int d = demand_matrix[i][j];
        for (; d --> 0 ;) {
            int c = H->colorForEdge(j);
            if (c >= 0) { // a common free color found
                H->removeEdge(j, c);
                edges_per_color[c][i] = j;
                edges_per_color[c][j] = i;
                for (int k = 0; k < T; k++) {
                    if (edges_per_color[c][k] == -1) {
                        G->addEdge(k, j, c);
                        H->removeEdge(k, c);
                    }
                }
                #ifdef ANALYZE_VIZING
                cases[0]++;
                #endif
                continue;
            }

            // create fan by finding a path in the FanGraph and FreeGraph
            vector<pair<int,int>> predecessors(T, make_pair(-1, -1));
            predecessors[j] = make_pair(j, -1);

            queue<int> q;
            q.push(j);
            bool path_found = false;
            while(!q.empty() && !path_found) {
                int u = q.front(); q.pop();
                for(int v = 0; v < T; v++) {
                    if (predecessors[v].first >= 0) continue;
                    int c_G = G->findEdge(u, v);
                    if (c_G == -1) continue;
                    predecessors[v] = make_pair(u, c_G);
                    q.push(v);
                    int c_H = H->colorForEdge(v);
                    if (c_H >= 0) {
                        predecessors[i] = make_pair(v, c_H);
                        path_found = true;
                        break;
                    }
                }
            }
            if (!path_found) {
                if (alg_type == GRAPH_1) {
                    num_failures += d+1;
                    break;
                }
                if (only_case1) {
                    case2_pairs.push_back(make_pair(i,j));
                    break;
                }
                int beta = graphFindFanCase2(i, j, G);
                int alpha = -1;
                for (int k = 0; k < T && alpha == -1; k++) {
                    alpha = H->colorForEdge(k); // this might not work for odd T, other way of finding free color of i is then required
                }
                if (alpha == -1) {
                    for (alpha=0; alpha < colors; alpha++) {
                        if (edges_per_color[alpha][i] == -1) {
                            break;
                        }
                    }
                }
                swapAugmentingPath(alpha, beta, i, G, H);
                #ifdef ANALYZE_VIZING
                cases[2]++;
                cases[1]--; // this "case 2" will be solved next iteration by case 1, but it stays case 2 for the statistics (though fan length may differ, I think..)
                #endif
                d++; continue; // finishing case 2 after swapping the augmenting path boils down to starting over with case 1
            }

            // extract fan (path) from predecessor data
            vector<pair<int,int>> fan;
            int v = i;
            while(v != j) {
                fan.push_back(predecessors[v]);
                v = predecessors[v].first;
            }

            // rotate fan and add update the edges in G accordingly (H was already updated)
            graphRotateFanCase1(i, fan, G, H);

            #ifdef ANALYZE_VIZING
            cases[1]++;
            fan_lengths.push_back(fan.size());
            #endif
        }
        demand_matrix[i][j] = demand_matrix[j][i] = d+1; // this +1 is what you get from weird syntax like goes to operator #noregrets
    }
    if (last_i != -1) {
        delete H;
        delete G;
    }
    return true;
}


// executes the rotation of the fan and updates FanGraph G accordingly
void HybridVizing::graphRotateFanCase1(int i, vector<pair<int,int>>& fan, FanGraph2* G, FreeGraph* H) {
    for (int jj = 0; jj < T; jj++) {
        if (i==jj) continue;
        if (edges_per_color[fan[0].second][jj] == -1) {
            H->removeEdge(jj, fan[0].second);
        }
    }
    
    int last_j = -1;
    for (pair<int,int> edge_hole_pair : fan) {
        int j = edge_hole_pair.first;
        int c = edge_hole_pair.second;

        if (last_j >= 0) {
            edges_per_color[c][last_j] = -1; // uncolor the edge incident to i that previously had this color 
            G->removeEdge(j, last_j, c); // edge has rotated, so G loses its edge [v,v_old]
        }

        for (int k = 0; k < T; k++) {
            if (k == i || k == j) continue;
            if (edges_per_color[c][k] == -1) {
                G->addEdge(k, j, c); // new edge {u,v} is added (or an existing one {u,v_old} has rotated to {u,v}), so G gets extra directed edges [k,v]
                if (last_j >= 0) {
                    G->removeEdge(k, last_j, c); // edge has rotated, so G loses its edges [k,v_old]
                }
            }
        }

        edges_per_color[c][i] = j;
        edges_per_color[c][j] = i;

        last_j = j;
    }
}


int HybridVizing::graphFindFanCase2(int u, int v, FanGraph2* G) {
    // BFS search of all paths between nodes (only paths reachable by v)
    vector<vector<int>> path_matrix(T, vector<int>(T, -1)); // path_matrix[i][j] = k --> node is is reachable by node j through node k
    queue<tuple<int,int,int>> q;
    vector<bool> has_been_reached(T, false);
    for (int ii = 0; ii < T; ii++) {
        if (ii==v) continue;
        if (G->findEdge(v, ii) == -1) continue;
        q.push(make_tuple(ii,v,v));
    }
    has_been_reached[v] = true;
    while(!q.empty()) {
        tuple<int,int,int> t = q.front(); q.pop();
        int ii = get<0>(t); int jj = get<1>(t); int kk = get<2>(t);
        if (path_matrix[ii][jj] >= 0) continue;
        path_matrix[ii][jj] = kk;
        for (int iii = 0; iii < T; iii++) {
            if (ii==iii) continue;
            if (G->findEdge(ii, iii) == -1) continue;
            if (!has_been_reached[ii]) {
                q.push(make_tuple(iii, ii, ii));
            }
            if (iii == jj) continue;
            if (path_matrix[iii][jj] >= 0) continue;
            q.push(make_tuple(iii, jj, ii));
        }
        has_been_reached[ii] = true;
    }

    // find a node w with a pair (i,j) of incident edges of the same color c
    int i = -1;
    int j;
    int w;
    int c;
    for (int ww = 0; ww < T && i == -1; ww++) {
        if (ww != v && path_matrix[ww][v] == -1) continue;
        for (int ii = 0; ii < T && i == -1; ii++) {
            if (G->findEdge(ii, ww) == -1) continue;
            if (ii != v && path_matrix[ii][v] == -1) continue;
            for (int jj = 0; jj < T; jj++) {
                if (jj == ii) continue;
                if (G->findEdge(jj, ww) == -1) continue;
                if (path_matrix[jj][ww] == -1) continue;
                if (G->Cx[ww][min(ii,jj)][max(ii,jj)].size() == 0) continue;
                c = *(G->Cx[ww][min(ii,jj)][max(ii,jj)].begin());
                i = ii; j = jj; w = ww;
                break;
            }
        }
    }
    return c;
}

void HybridVizing::swapAugmentingPath(int x, int y, int i, FanGraph2* G, FreeGraph* H) {
    // first edge
    int j = edges_per_color[y][i];
    for (int k = 0; k < T; k++) {
        if (k == i || k == j) continue;
        if (edges_per_color[y][k] == -1) {
            G->removeEdge(k, j, y);
            H->addEdge(k, y);
        }
        if (edges_per_color[x][k] == -1) {
            G->addEdge(k, j, x);
            H->removeEdge(k, x);
        }
    }

    // alternate xy-path
    int ii = i;
    int jj = j;
    bool odd_edge = false;
    while (jj != -1) {
        swap(edges_per_color[x][ii], edges_per_color[y][ii]);
        ii = jj;
        jj = edges_per_color[odd_edge?y:x][jj];
        odd_edge = !odd_edge;
    }
    swap(edges_per_color[x][ii], edges_per_color[y][ii]);

    // last edge
    if (!odd_edge) {
        G->addEdge(ii, j, x); // the last node of the path loses color x, so G must get an edge [ii,j] of color x
        H->removeEdge(ii, y); // color y is no longer free in node ii, so the edge in H added above should be removed again already
    } else {
        G->removeEdge(ii, j, x);
        H->addEdge(ii, y);
    }
}
