#include "help_structures.hpp"

using namespace std;

FreeGraph::FreeGraph(vector<vector<int>>& edges_per_color, int T, int colors, int i) : T(T), colors(colors), free_edges(T) {
    for (int c = 0; c < colors; c++) {
        if (edges_per_color[c][i] >= 0) continue;
        for (int j = 0; j < T; j++) {
            if (j == i) continue;
            if (edges_per_color[c][j] == -1) {
                addEdge(j, c);
            }
        }
    }
}

// returns the color of edge {i,v} if there is one, else -1
int FreeGraph::colorForEdge(int v) {
    return free_edges[v].size() > 0 ? *(free_edges[v].begin()) : -1;
}

void FreeGraph::addEdge(int v, int c) {
    free_edges[v].insert(c);
}
void FreeGraph::removeEdge(int v, int c) {
    free_edges[v].erase(c);
}



FanGraph2::FanGraph2(vector<vector<int>>& edges_per_color, int T, int colors, int i, bool only_case1) : T(T), colors(colors), edges(T, vector<set<int>>(T)), only_case1(only_case1) {
    if (!only_case1) {
        Cx = vector<vector<vector<set<int>>>>(T, vector<vector<set<int>>>(T, vector<set<int>>(T)));
    }
    
    for (int c = 0; c < colors; c++) {
        int j = edges_per_color[c][i];
        if (j == -1) continue;
        for (int k = 0; k < T; k++) {
            if (edges_per_color[c][k] >= 0) continue;
            addEdge(k, j, c);
            if (only_case1) continue;
            for (int l = k+1; l < T; l++) {
                if (edges_per_color[c][l] == -1) {
                    Cx[j][k][l].insert(c);
                }
            }
        }
    }
}

void FanGraph2::addEdge(int u, int v, int c) {
    edges[u][v].insert(c);
    if (only_case1) return;
    for (int k = 0; k < T; k++) {
        if (k == u) continue;
        if (!hasEdge(k, v, c)) continue;
        Cx[v][min(u,k)][max(u,k)].insert(c);
    }
}
void FanGraph2::removeEdge(int u, int v, int c) {
    edges[u][v].erase(c);
    if (only_case1) return;
    for (int k = 0; k < T; k++) {
        if (k == u) continue;
        Cx[v][min(u,k)][max(u,k)].erase(c);
    }
}

int FanGraph2::findEdge(int u, int v) {
    return edges[u][v].size() > 0 ? *(edges[u][v].begin()) : -1;
}
int FanGraph2::findEdgeRestricted(int u, int v, set<int>& avoid_doubles) {
    for (int c : edges[u][v]) {
        if (avoid_doubles.count(c) > 0) continue;
        return c;
    }
    return -1;
}

bool FanGraph2::hasEdge(int u, int v, int c) {
    return edges[u][v].count(c);
}
