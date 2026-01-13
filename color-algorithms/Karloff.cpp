#include "Karloff.hpp"
#include "InitializeTestInput.hpp"

using namespace std;

Karloff::Karloff(vector<vector<int>>& demand_matrix, const int T, const int S) : demand_matrix(demand_matrix), T(T), lb(S),
bipartite_demand_matrix(T, vector<int>(T)) {}


vector<vector<int>> Karloff::color() {
    eulerTour_Hierholzer();
    bipartiteEC();
    finalColoring();

    return edges_per_color;
}

void Karloff::eulerTour_Hierholzer() {
    vector<list<ListNode>> nodes = getEdgeLists();
    list<pair<int,int>> Euler_tour;
    
    queue<pair<int,list<pair<int,int>>::iterator>> insert_positions; // for each vertex u, remember some (for example the first) position it appears in the Euler tour to add extra Euler tour later that starts and ends with u
    insert_positions.push({0, Euler_tour.begin()});
    vector<bool> position_initialized(T, false); position_initialized[0] = true;
    
    list<pair<int,int>> partial_Euler_tour;
    while (!insert_positions.empty()) {
        int i = insert_positions.front().first;
        list<pair<int,int>>::iterator current_insert_position = insert_positions.front().second;
        insert_positions.pop();
        int u = i;
        while (nodes[u].size() > 0) {
            int v = nodes[u].front().destination;
            list<ListNode>::iterator reverse_edge_it = nodes[u].front().destination_list_it;
            nodes[u].pop_front();
            nodes[v].erase(reverse_edge_it);
            partial_Euler_tour.push_back({u,v});
            if (u < T && !position_initialized[u]) {
                insert_positions.push({u, prev(partial_Euler_tour.end())});
                position_initialized[u] = true;
            }
            u = v;
        }
        Euler_tour.splice(current_insert_position, partial_Euler_tour); // the iterators of partial_Euler_tour remain valid in Euler_tour
    }
    
    for (pair<int,int> se : Euler_tour) {
        int u = se.first;
        int v = se.second;
        if (u == T || v == T) continue;
        bipartite_demand_matrix[u][v]++;
    }
}

vector<list<ListNode>> Karloff::getEdgeLists() {
    int degrees[T];
    for (int i = 0; i < T; i++) {
        degrees[i] = 0;
        for (int j = 0; j < T; j++) {
            degrees[i] += demand_matrix[i][j];
        }
    }
    
    vector<list<ListNode>> nodes(T+1);
    for (int i = 0; i < T ; i++) {
        for (int j = i+1; j < T; j++) {
            for (int d = 0; d < demand_matrix[i][j]; d++) {
                nodes[i].push_back(ListNode(j));
                nodes[j].push_back(ListNode(i));
                nodes[i].back().destination_list_it = prev(nodes[j].end());
                nodes[j].back().destination_list_it = prev(nodes[i].end());
            }
        }
        if (degrees[i] % 2 == 1) { // add dummy edge to dummy vertex
            nodes[i].push_back(ListNode(T));
            nodes[T].push_back(ListNode(i));
            nodes[i].back().destination_list_it = prev(nodes[T].end());
            nodes[T].back().destination_list_it = prev(nodes[i].end());
        }
    }
    return nodes;
}

void Karloff::bipartiteEC() {
    bipartite_colors = (lb+1)/2; // max degree of bipartite matrix is ceil(max degree / 2)
    BipartiteEdgeColor b_ec(bipartite_demand_matrix, T, bipartite_colors);
    bipartite_solution = b_ec.color();
    bipartite_solution = b_ec.convertICJtoCIJ();
}

void Karloff::finalColoring() {
    for (int c = 0; c < bipartite_colors; c++) {
        vector<int> matching = bipartite_solution[c];
        vector<int> reverse_matching(T, -1);
        for (int i = 0; i < T; i++) {
            int j = matching[i];
            if (j >= 0) {
                reverse_matching[j] = i;
            }
        }
        vector<vector<int>> double_c(2, vector<int>(T, -1));
        vector<int> extra_c(T, -1);
        int num_extra_c = 0;
        for (int i = 0; i < T; i++) {
            int current_c = 0;
            int current_rc = 1;
            int j = matching[i];
            int last_j = i;
            int rj = reverse_matching[i];
            int last_rj = i;
            if (rj >= 0) matching[rj] = -1;
            if (j >= 0) reverse_matching[j] = -1;
            while (j >= 0 || rj >= 0) {
                if (j == rj) { // end of even cycle: no third color needed
                    double_c[current_c][last_j] = j;
                    double_c[current_c][j] = last_j;
                    matching[last_j] = -1;
                    matching[j] = -1;
                    double_c[current_rc][last_rj] = rj;
                    double_c[current_rc][rj] = last_rj;
                    reverse_matching[last_rj] = -1;
                    reverse_matching[rj] = -1;
                    break;
                } else if (j == last_rj) { // also (rj == last_j) // end of odd cycle: third color needed
                    extra_c[last_j] = j;
                    extra_c[j] = last_j;
                    matching[last_j] = -1;
                    matching[j] = -1;
                    num_extra_c++;
                    break;
                }
                if (j >= 0) {
                    double_c[current_c][last_j] = j;
                    double_c[current_c][j] = last_j;
                    matching[last_j] = -1;
                    last_j = j;
                    j = matching[last_j];
                    current_c = 1 - current_c;
                    if (j >= 0) reverse_matching[j] = -1;
                }
                if (rj >= 0) {
                    double_c[current_rc][last_rj] = rj;
                    double_c[current_rc][rj] = last_rj;
                    reverse_matching[last_rj] = -1;
                    last_rj = rj;
                    rj = reverse_matching[last_rj];
                    current_rc = 1 - current_rc;
                    if (rj >= 0) matching[rj] = -1;
                }
            }
            matching[i] = -1;
            reverse_matching[i] = -1;
        }
        edges_per_color.push_back(double_c[0]);
        edges_per_color.push_back(double_c[1]);
        if (num_extra_c > 0) {
            edges_per_color.push_back(extra_c);
        }
    }
}




BipartiteEdgeColor::BipartiteEdgeColor(std::vector<std::vector<int>>& demand_matrix, int T, int S) : num_excess_colors(demand_matrix), T(T), S(S) {
    edges_per_color = vector<vector<int>>(T, vector<int>(S,-1));
    edges_per_color_reverse = vector<vector<int>>(T, vector<int>(S,-1));
}



vector<vector<int>> BipartiteEdgeColor::color() {
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++) {
            int excess = num_excess_colors[i][j];
            int alpha = -1; int beta = -1;
            int c_i = 0; int c_j = 0;
            while (excess --> 0) {
                for (; c_i < S; c_i++) {
                    if (edges_per_color[i][c_i] == -1) {
                        alpha = c_i++;
                        break;
                    }
                }
                for (; c_j < S; c_j++) {
                    if (edges_per_color_reverse[j][c_j] == -1) {
                        // c_j needs no increment in case that beta is not "used" in augment, meaning that alpha was not incident to j
                        beta = c_j;
                        break;
                    }
                }
                augment(i, j, alpha, beta);
            }
        }
    }
    return edges_per_color;
}

// alpha is absent color from v, beta is absent color from w
// an edge (v,w) with color alpha will be added
// v is in the left part of the bipartite graph, w in the right part
void BipartiteEdgeColor::augment(int v, int w, int alpha, int beta) {
    int vertex = w;
    int next_vertex = edges_per_color_reverse[w][alpha];
    int next_next_vertex;
    bool left = false;

    while (next_vertex >= 0) {
        if(!left){
            next_next_vertex = edges_per_color[next_vertex][beta];

            edges_per_color[next_vertex][beta] = vertex;
            edges_per_color_reverse[vertex][beta] = next_vertex;

            edges_per_color[next_vertex][alpha] = -1;
        } else {
            next_next_vertex = edges_per_color_reverse[next_vertex][alpha];

            edges_per_color_reverse[next_vertex][alpha] = vertex;
            edges_per_color[vertex][alpha] = next_vertex;

            edges_per_color_reverse[next_vertex][beta] = -1;
        }
        vertex = next_vertex;
        next_vertex = next_next_vertex;
        left = !left;
    }
    edges_per_color[v][alpha] = w;
    edges_per_color_reverse[w][alpha] = v;
}

// edges_per_color has the indexing order [i][c][j], rather than the desired [c][i][j] (which we fix in the function)
// this is done because the bipartite edge coloring algorithm above executes faster with the [i][c][j] indexing order due to cache locality
vector<vector<int>> BipartiteEdgeColor::convertICJtoCIJ() {
    vector<vector<int>> c_first(S, vector<int>(T, -1));
        for (int c = 0; c < S; c++) {
            for (int i = 0; i < T; i++) {
                int j = edges_per_color[i][c];
                if (j == -1) {
                    continue;
                }
                c_first[c][i] = j;
            }
        }
    return c_first;
}
