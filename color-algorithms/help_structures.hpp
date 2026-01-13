#ifndef HELP_STRUCTURES
#define HELP_STRUCTURES

#include <vector>
#include <set>

#include <cassert>

// this data structure corresponds to all the type 2 arcs of $D^u$ in the paper
class FreeGraph {
    public:
    FreeGraph(std::vector<std::vector<int>>& edges_per_color, int T, int colors, int i);

    int colorForEdge(int v);

    void addEdge(int v, int c);
    void removeEdge(int v, int c);

    std::vector<std::set<int>> free_edges;

    int T;
    int colors;
};

// this data structure corresponds to $D^u$ in the paper
// the data structure $C^x$ is a member of this class, but it is only initialized and used in DiGraphColor2
class FanGraph2 {
    public:
    FanGraph2(std::vector<std::vector<int>>& edges_per_color, int T, int colors, int i, bool only_case1);

    void addEdge(int u, int v, int c);
    void removeEdge(int u, int v, int c);
    int findEdge(int u, int v);
    int findEdgeRestricted(int u, int v, std::set<int>& avoid_doubles);
    bool hasEdge(int u, int v, int c);

    std::vector<std::vector<std::set<int>>> edges;

    std::vector<std::vector<std::vector<std::set<int>>>> Cx;

    int T;
    int colors;
    
    bool only_case1;
};

#endif