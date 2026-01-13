#include "InitializeTestInput.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

InitializeTestInput::InitializeTestInput() {
    srand(10);
}

// code of this function is based on the source code by Ziyue Zhang in the work https://doi.org/10.1364/JOCN.474187
int InitializeTestInput::initializeApplicationInput(vector<vector<int>>& demand_matrix, int T, int S, int num_peaks) {
    int maxD = ceil(0.5*(float)S/T);
    for (int i = 0; i < T; i++) {
        for (int j = i+1; j < T; j++) {
            demand_matrix[i][j] = demand_matrix[j][i] = rand() % (maxD +1) ;
        }
    }  
    //==============set some peak demand_matrixs===========
    for (int i = 0; i < num_peaks; i++) {
        int peak_x=rand() % T;
        int peak_y=rand() % (T-1);
        if (peak_y >= peak_x) {
            peak_y++;
        }
        demand_matrix[peak_x][peak_y] += maxD*10;
        demand_matrix[peak_y][peak_x] += maxD*10;
    }

    vector<vector<int>> orig_demand(demand_matrix);
    phase_1(demand_matrix, S, T);
    phase_2(orig_demand, demand_matrix, 0.0, S, T);

    initializations_done++;
    return upperBound(demand_matrix, T);
}

void InitializeTestInput::phase_1(vector<vector<int>>& d, unsigned int N, unsigned int T) {
    int sum_max = maxDegree(d, T);

    for (size_t i = 0; i < T; i++) {
        for (size_t j = 0; j < T; j++) {
            if (i==j) { 
                continue;
            }
            //(N-(T-1)) because we take T-1 wavelength channels for a full mesh connection
            d[i][j] = (d[i][j]*(N-T+1)/sum_max) + 1;
        }
    }
}

void InitializeTestInput::phase_2(vector<vector<int>>& orig_d, vector<vector<int>>& a, float alpha,unsigned int N,unsigned int T){

    vector<size_t> sum_row(T,0);

    multimap<float, pair<int,int>> Q_m;
    for (size_t i = 0; i < T; i++)
    {
        for (size_t j = 0; j < T; j++)
        {
            sum_row[i]+=a[i][j];
            #ifdef DEBUG
            assert(sum_col[i]<=N && sum_row[j]<=N);
            #endif
        }
    }
    for (size_t i = 0; i < T; i++)
    {
        for (size_t j = 0; j < T; j++)
        {
            if (sum_row[j]<N && sum_row[i]<N && i!=j && orig_d[i][j]>0)
            {
                Q_m.insert(make_pair((a[i][j]-orig_d[i][j])/pow(orig_d[i][j],alpha) , make_pair(i,j)));
            }
        }
    }
    while (!Q_m.empty())
    {
        auto it=Q_m.begin();//the minimal margin
        if (sum_row[it->second.second]==N || sum_row[it->second.first]==N)
        {
            it=Q_m.erase(it);
            continue;
        }
        
        int i = it->second.first;
        int j = it->second.second;
        a[i][j]++;
        a[j][i]++;
        sum_row[j]++;
        sum_row[i]++;
        #ifdef DEBUG
        assert(sum_col[i]<=N && sum_row[j]<=N);
        #endif
        if (sum_row[j]==N || sum_row[i]==N)
        {
            it=Q_m.erase(it);
        }else{
            float new_m=it->first+1/pow(orig_d[i][j],alpha);
            it=Q_m.erase(it);
            Q_m.insert(make_pair(new_m, make_pair(i,j)));
        }
    }
}


int InitializeTestInput::randomERGraph(vector<vector<int>>& demand_matrix, int T, int m) {
    for (int i = 0; i < T; i++) {
        for (int j = i+1; j < T; j++) {
            demand_matrix[i][j] = demand_matrix[j][i] = 0;
        }
    }
    for (; m --> 0 ;) {
        int i = rand() % T;
        int j = rand() % (T-1); if (j >= i) j++;
        demand_matrix[i][j] = demand_matrix[j][i] += 1;
    }
    
    return upperBound(demand_matrix, T);
}


// method to check for a demand matrix what the max degree of the corresponding graph is
int InitializeTestInput::maxDegree(vector<vector<int>>& demand_matrix, int T, int& max_i) {
    int max_degree = 0;

    // check degrees of left nodes of bipartite graph (rows in the matrix)
    for (int i = 0; i < T; i++) {
        int outgoing_degree = 0;
        for (int j = 0; j < T; j++) {
            outgoing_degree += demand_matrix[i][j];
        }
        if (outgoing_degree > max_degree) {
            max_degree = outgoing_degree;
            max_i = i;
        }
    }
    return max_degree;
}
int InitializeTestInput::maxDegree(vector<vector<int>>& demand_matrix, int T) {
    int dummy = 0;
    return maxDegree(demand_matrix, T, dummy);
}

// method to check for a demand matrix the Delta+µ upper bound
int InitializeTestInput::upperBound(vector<vector<int>>& demand_matrix, int T) {
    int max_degree = 0;
    int max_multiplicity = 0;

    // check degrees of all vertices (rows in the matrix, columns would be the same because symmetric)
    for (int i = 0; i < T; i++) {
        int degree = 0;
        for (int j = 0; j < T; j++) {
            int d = demand_matrix[i][j];
            max_multiplicity = max(max_multiplicity, d);
            degree += d;
        }
        max_degree = max(max_degree, degree);
    }

    return max_degree+max_multiplicity;
}


void InitializeTestInput::showMatrix(const vector<vector<int>>& demand_matrix, int T) {
    int w = 1;
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++) {
            while (demand_matrix[i][j] >= (int)pow(10, w)) {
                w++;
            }
        }
    }

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++) {
            cout << setw(w) << demand_matrix[i][j] << " ";
        }
        cout << endl;
    }
}