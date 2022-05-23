// https://open.kattis.com/problems/maxflow

#include <iostream>
// runtime 1.75s

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
    
    ln::network_flow_shortest_path_scaling f(G);
    f.set_capacity(capacity);
    
    const int max_flow = f.send(S,T);
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
