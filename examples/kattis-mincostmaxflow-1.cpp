// https://open.kattis.com/problems/mincostmaxflow
#include "network.hpp"
#include <iostream>
// runtime: 0.86s

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
    
    ln::network_flow_cost_EdmondsKarp f(G);
    f.set_capacity(capacity);
    
    auto Flow = f.send(S,T,weight);
    long long Cost = 0;
    
    for(int e=0;e<M;++e)
    {
        Cost+= static_cast<long long>(weight[e])*f.flow_at(e);
    }
    
    std::cout << Flow << " " << Cost << '\n';
    
    return 0;
}

