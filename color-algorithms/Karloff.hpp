#ifndef KARLOFF
#define KARLOFF

#include "InitializeTestInput.hpp"

#include <vector>
#include <list>
#include <queue>

enum EulerAlg {
    BEST, BEST_Z, HIERHOLZER, HIERHOLZER_Z
};

struct ListNode {
    ListNode(int destination) : destination(destination) {};

    int destination;
    std::list<ListNode>::iterator destination_list_it;
};

class Karloff {
    public:
    Karloff(std::vector<std::vector<int>>& demand_matrix, const int T, const int max_degree);

    std::vector<std::vector<int>> color();

    void bipartiteEC();
    void finalColoring();

    void eulerTour_Hierholzer();
    std::vector<std::list<ListNode>> getEdgeLists();
    
    std::vector<std::vector<int>>& demand_matrix;
    const int T;
    int lb;

    std::vector<std::vector<int>> bipartite_demand_matrix;
    std::vector<std::vector<int>> bipartite_solution;
    int bipartite_colors;

    std::vector<std::vector<int>> edges_per_color;

};

class BipartiteEdgeColor {
    public:
    BipartiteEdgeColor(std::vector<std::vector<int>>& demand_matrix, int T, int S);

    std::vector<std::vector<int>> color();
    void augment(int v, int w, int alpha, int beta);

    std::vector<std::vector<int>> convertICJtoCIJ();

    int T;
    int S;

    std::vector<std::vector<int>> edges_per_color;
    std::vector<std::vector<int>> edges_per_color_reverse;

    std::vector<std::vector<int>> num_excess_colors;
};

#endif
