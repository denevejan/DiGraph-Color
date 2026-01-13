#include "DiGraphColor2.hpp"

#include <iostream>

using namespace std;

DiGraphColor2::DiGraphColor2(const std::vector<std::vector<int>>& demand_matrix, const int T, const int ub) : demand_matrix(demand_matrix), T(T), colors(ub) {
    edges_per_color = vector<vector<int>>(colors, vector<int>(T, -1));
};

vector<vector<int>> DiGraphColor2::color() {
    for (int i = 0; i < T-1; i++) {
        FreeGraph H(edges_per_color, T, colors, i);
        FanGraph2 G(edges_per_color, T, colors, i, false);
        for (int j = i+1; j < T; j++) {
            int d = demand_matrix[i][j];
            for (; d --> 0 ;) {

                int c = H.colorForEdge(j);
                if (c >= 0) { // a common free color found
                    H.removeEdge(j, c);
                    edges_per_color[c][i] = j;
                    edges_per_color[c][j] = i;
                    for (int k = 0; k < T; k++) {
                        if (k==i && k==j) continue;
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
                    int beta = findFanCase2(i, j, G);
                    int alpha = -1;
                    for (int k = 0; k < T && alpha == -1; k++) {
                        alpha = H.colorForEdge(k);
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
                    cases[2]++; cases[1]--; // this "case 2" will be solved next iteration by case 1, but it stays case 2 for the statistics
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

                // rotate fan and add update the edges in G and H accordingly
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

// executes the rotation of the fan and updates G and H accordingly
void DiGraphColor2::rotateFanCase1(int i, vector<pair<int,int>>& fan, FanGraph2& G, FreeGraph& H) {
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

int DiGraphColor2::findFanCase2(int u, int v, FanGraph2& G) {

    // BFS search of all paths between nodes (only paths reachable by v)
    vector<vector<int>> path_matrix(T, vector<int>(T, -1)); // path_matrix[i][j] = k --> node is is reachable by node j through node k
    queue<tuple<int,int,int>> q;
    vector<bool> has_been_reached(T, false);
    for (int ii = 0; ii < T; ii++) {
        if (ii==v) continue;
        if (G.findEdge(v, ii) == -1) continue;
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
            if (G.findEdge(ii, iii) == -1) continue;
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
            if (G.findEdge(ii, ww) == -1) continue;
            if (ii != v && path_matrix[ii][v] == -1) continue;
            for (int jj = 0; jj < T; jj++) {
                if (jj == ii) continue;
                if (G.findEdge(jj, ww) == -1) continue;
                if (path_matrix[jj][ww] == -1) continue;
                if (G.Cx[ww][min(ii,jj)][max(ii,jj)].size() == 0) continue;
                c = *(G.Cx[ww][min(ii,jj)][max(ii,jj)].begin());
                i = ii; j = jj; w = ww;
                break;
            }
        }
    }

    return c;
}

void DiGraphColor2::swapAugmentingPath(int x, int y, int i, FanGraph2& G, FreeGraph& H) {
    // first edge
    int j = edges_per_color[y][i];
    for (int k = 0; k < T; k++) {
        if (k == i || k == j) continue;
        if (edges_per_color[y][k] == -1) {
            G.removeEdge(k, j, y);
            H.addEdge(k, y);
        }
        if (edges_per_color[x][k] == -1) {
            G.addEdge(k, j, x);
            H.removeEdge(k, x);
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
        G.addEdge(ii, j, x); // the last node of the path loses color x, so G must get an edge [ii,j] of color x
        H.removeEdge(ii, y); // color y is no longer free in node ii, so the edge in H added above should be removed again already
    } else {
        G.removeEdge(ii, j, x);
        H.addEdge(ii, y);
    }
}
