// https://open.kattis.com/problems/shortestpath1

#include "mincostflow/graph.hpp"
#include "mincostflow/paths.hpp"
#include <iostream>
// runtime 0.23s

int main()
{
    while(1)
    {
        int N_vertex,N_edges,Q,S;
        std::cin>>N_vertex>>N_edges>>Q>>S;
        if(N_vertex==0)break;
        
        ln::digraph Graph(N_vertex);
        std::vector<int> weights(N_edges);
        
        
        for(int e=0;e<N_edges;++e)
        {
            int a,b,w;
            std::cin>> a>>b>>w;
            Graph.add_edge(a,b);
            weights.at(e) = w;
        }
        
        ln::shortest_path_dijkstra solver(Graph);
        solver(S,weights,[](int e){return true;});
        
        for(int v;Q--;)
        {
            std::cin>>v; 
            if(!solver.is_connected(v))
                std::cout << "Impossible\n";
            else
            {
                auto d = solver.get_distance(v); 
                std::cout << d << '\n';
            }
        }
        
        std::cout<<'\n';
    }
    return 0;
}
