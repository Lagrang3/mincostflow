// https://open.kattis.com/problems/mincostmaxflow
#include <mincostflow/mincostflow.hpp>
#include <iostream>

typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_FIFO> mincostflow_t; // 0.94s
// typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_BellmanFord> mincostflow_t; // 1.98s
// typedef ln::mincostflow_EdmondsKarp<ln::shortestPath_Dijkstra> mincostflow_t; // expected failure

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

