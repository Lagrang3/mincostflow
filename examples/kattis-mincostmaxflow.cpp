// https://open.kattis.com/problems/mincostmaxflow
#include <mincostflow/mincostflow.hpp>
#include <iostream>

typedef ln::maxflow_augmenting_path<ln::pathSearch_BFS> maxflow_t;
// typedef ln::maxflow_augmenting_path<ln::pathSearch_labeling> maxflow_t;
// typedef ln::maxflow_scaling<ln::pathSearch_BFS> maxflow_t;
// typedef ln::maxflow_scaling<ln::pathSearch_labeling> maxflow_t;
// typedef ln::maxflow_preflow maxflow_t;

// typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_FIFO> mincostflow_t; // 0.94s
// typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_BellmanFord> mincostflow_t; // 1.98s
// typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_Dijkstra> mincostflow_t; // expected failure

// typedef ln::mincostflow_PrimalDual<ln::shortestPath_FIFO,maxflow_t> mincostflow_t; // 0.95s, 1.11s, 1.18s, 1.38s, TLE
// typedef ln::mincostflow_PrimalDual<ln::shortestPath_BellmanFord,maxflow_t> mincostflow_t; // 2.13s, 2.23s, 2.20s, 2.47s, TLE
typedef ln::mincostflow_PrimalDual<ln::shortestPath_Dijkstra,maxflow_t> mincostflow_t; // 0.63s, 0.74s, 0.87s, 1.07s, TLE

int main()
{
    int N,M,S,T;
    std::cin >> N >> M >> S >> T;
    
    ln::digraph G(N);
    std::vector<int> capacity(M,0);
    std::vector<int> weight(M,0);
    
    for(int e=0;e<M;++e)
    {
        int a,b,c,w;
        std::cin>>a>>b>>c>>w;
        G.add_edge(a,b);
        capacity[e] = c;
        weight[e] = w;
    }
    
    mincostflow_t f(G);
    f.set_capacity(capacity);
    
    auto Flow = f.solve(S,T,weight);
    long long Cost = 0;
    
    for(int e=0;e<M;++e)
    {
        Cost+= static_cast<long long>(weight[e])*f.flow_at(e);
    }
    
    std::cout << Flow << " " << Cost << '\n';
    
    return 0;
}

