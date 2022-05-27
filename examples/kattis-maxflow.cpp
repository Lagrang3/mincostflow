// https://open.kattis.com/problems/maxflow

#include <mincostflow/graph.hpp>
#include <mincostflow/shortestpath.hpp>
#include <mincostflow/maxflow.hpp>
#include <iostream>

// typedef ln::maxflow_augmenting_path<ln::pathSearch_BFS> maxflow_t; // 0.43s
// typedef ln::maxflow_augmenting_path<ln::pathSearch_labeling> maxflow_t; // 0.14s
// typedef ln::maxflow_scaling<ln::pathSearch_BFS> maxflow_t; // 0.15s
typedef ln::maxflow_scaling<ln::pathSearch_labeling> maxflow_t; // 0.02s
// typedef ln::maxflow_preflow maxflow_t;

int main()
{
    int N,M,S,T;
    std::cin >> N >> M >> S >> T;
    
    ln::digraph G(N);
    std::vector<int> capacity(M,0);
    
    for(int e=0;e<M;++e)
    {
        int a,b,c;
        std::cin>>a>>b>>c;
        G.add_edge(a,b);
        capacity[e] = c;
    }
    
    maxflow_t f(G);
    f.set_capacity(capacity);
    
    const int max_flow = f.solve(S,T);
    int M_count  =0 ;
    for(int e=0;e<G.n_edges();++e)
    {
        int my_f = f.flow_at(e);
        M_count += (my_f > 0 ? 1 : 0);
    }
    
    std::cout << N << ' ' << max_flow  << ' ' << M_count << '\n';
    
    for(int e=0;e<G.n_edges();++e)
    {
        int my_f = f.flow_at(e);
        if(my_f==0)
            continue;
            
        auto [a,b] = G.get_edge(e);
        std::cout << a << ' ' << b << ' ' << my_f << '\n';
    }
    
    return 0;
}
