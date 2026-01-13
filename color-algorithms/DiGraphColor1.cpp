#include "DiGraphColor1.hpp"

#include <iostream>

using namespace std;


DiGraphColor1::DiGraphColor1(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub) : demand_matrix(demand_matrix), T(T), colors(ub) {
    edges_per_color = vector<vector<int>>(colors, vector<int>(T, -1));
}

vector<vector<int>> DiGraphColor1::color() {
    for (int i = 0; i < T-1; i++) {
        FreeGraph H(edges_per_color, T, colors, i);
        FanGraph2 G(edges_per_color, T, colors, i, true);
        for (int j = i+1; j < T; j++) {
            int d = demand_matrix[i][j];
            for (; d --> 0 ;) {
                int c = H.colorForEdge(j);
                if (c >= 0) { // a common free color found
                    H.removeEdge(j, c);
                    edges_per_color[c][i] = j;
                    edges_per_color[c][j] = i;
                    for (int k = 0; k < T; k++) {
                        if (edges_per_color[c][k] == -1) {
                            G.addEdge(k, j, c);
                            H.removeEdge(k, c);
                        }
                    }
                    #ifdef ANALYZE_VIZING
                    cases[0]++;
                    #endif
                    continue;
                }

                // create fan by finding a path in the FanGraph2 and FreeGraph
                vector<pair<int,int>> predecessors(T, make_pair(-1, -1));
                predecessors[j] = make_pair(j, -1);

                queue<int> q;
                q.push(j);
                bool path_found = false;
                while(!q.empty() && !path_found) {
                    int u = q.front(); q.pop();
                    for(int v = 0; v < T; v++) {
                        if (predecessors[v].first >= 0) continue;
                        int c_G = G.findEdge(u, v);
                        if (c_G == -1) continue;
                        predecessors[v] = make_pair(u, c_G);
                        q.push(v);
                        int c_H = H.colorForEdge(v);
                        if (c_H >= 0) {
                            predecessors[i] = make_pair(v, c_H);
                            path_found = true;
                            break;
                        }
                    }
                }
                if (!path_found) {
                    num_failures += d+1;
                    break;
                }

                // extract fan (path) from predecessor data
                vector<pair<int,int>> fan;
                int v = i;
                while(v != j) {
                    fan.push_back(predecessors[v]);
                    v = predecessors[v].first;
                }

                // rotate fan and add update the edges in G accordingly (H was already updated)
                rotateFanCase1(i, fan, G, H);

                #ifdef ANALYZE_VIZING
                cases[1]++;
                fan_lengths.push_back(fan.size());
                #endif

            }
        }
    }
    return edges_per_color;
}

// executes the rotation of the fan and updates FanGraph2 G accordingly
void DiGraphColor1::rotateFanCase1(int i, vector<pair<int,int>>& fan, FanGraph2& G, FreeGraph& H) {
    for (int jj = 0; jj < T; jj++) {
        if (i==jj) continue;
        if (edges_per_color[fan[0].second][jj] == -1) {
            H.removeEdge(jj, fan[0].second);
        }
    }
    
    int last_j = -1;
    for (pair<int,int> edge_hole_pair : fan) {
        int j = edge_hole_pair.first;
        int c = edge_hole_pair.second;

        if (last_j >= 0) {
            edges_per_color[c][last_j] = -1; // uncolor the edge incident to i that previously had this color 
            G.removeEdge(j, last_j, c); // edge has rotated, so G loses its edge [v,v_old]
        }

        for (int k = 0; k < T; k++) {
            if (k == i || k == j) continue;
            if (edges_per_color[c][k] == -1) {
                G.addEdge(k, j, c); // new edge {u,v} is added (or an existing one {u,v_old} has rotated to {u,v}), so G gets extra directed edges [k,v]
                if (last_j >= 0) {
                    G.removeEdge(k, last_j, c); // edge has rotated, so G loses its edges [k,v_old]
                }
            }
        }

        edges_per_color[c][i] = j;
        edges_per_color[c][j] = i;

        last_j = j;
    }
}
