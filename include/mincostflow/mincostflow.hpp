#pragma once

namespace ln
{
    template<typename base_network_flow>
    class network_flow_cost_solver : public base_network_flow
    {
        public:
        
        network_flow_cost_solver(const digraph& in_graph):
            base_network_flow(in_graph)
        {}
    };
    
    class network_flow_cost_EdmondsKarp : 
        public network_flow_cost_solver<network_flow_EdmondsKarp>
    {
        using base_network = network_flow_cost_solver<network_flow_EdmondsKarp>;
        
        int greedy_augmenting_path(
            const int Source, const int Dest,
            const std::vector<int> weight,
            int Flow_demand = INF)
        // augmenting path
        // kattis mincostflow 0.86s
        {   
            std::vector<int> weight_ex(nedges*2);
            for(int e=0;e<nedges;++e)
            {
                weight_ex.at(e) = weight.at(e);
                weight_ex.at(dual(e)) = -weight.at(e);
            } 
            int sent =0 ;
            shortest_path_bfs bfs(Graph);
            
            while(Flow_demand>0)
            {
                bfs(Source,
                    weight_ex,
                    // edge is valid if
                    [this](int e){
                        return residual_cap.at(e)>0;
                    });
                
                auto path = bfs.find_path(Dest);
                
                if(path.empty())
                    break;
                    
                int k = Flow_demand;
                for(auto e : path)
                {
                    k = std::min(k,residual_cap.at(e));
                }
                
                for(auto e: path)
                {
                    residual_cap[e] -= k;
                    residual_cap[dual(e)] += k;
                } 
                
                sent += k;
                Flow_demand -= k;
            }
            return sent;
        }
        
        public:
        network_flow_cost_EdmondsKarp(const digraph& in_graph):
            base_network{in_graph}    
        {}
    
        int send(
            const int Source, const int Dest,
            const std::vector<int> weight,
            int Flow_demand = INF)
        {
            if(weight.size()!=nedges)
                throw std::runtime_error(
                    "send: weight.size() != "
                    "number of edges");
            return greedy_augmenting_path(Source,Dest,weight,Flow_demand); 
        }
    };
    
    
    // TODO: 
    class network_flow_cost_PrimalDual : 
        public network_flow_cost_solver<network_flow_PushRelabel>
    {
        using base_network = network_flow_cost_solver<network_flow_PushRelabel>;
        
        int primal_dual(
            const int Source, const int Dest,
            const std::vector<int> weight,
            int Flow_demand = INF)
        // Primal-Dual
        {   
            std::vector<int> weight_ex(nedges*2);
            std::vector<int> potential(Graph.n_vertex(),0);
            
            for(int e=0;e<nedges;++e)
            {
                weight_ex.at(e) = weight.at(e);
                weight_ex.at(dual(e)) = -weight.at(e);
            } 
            
            int sent =0 ;
            shortest_path_dijkstra dijkstra(Graph);
            while(Flow_demand>0)
            {
                dijkstra(
                    Source, 
                    weight_ex,
                    // edge is valid if
                    [this](int e) -> bool
                    {
                        return residual_cap.at(e)>0;
                    });
                const auto& distance{dijkstra.distance};
                
                for(int e = 0 ;e<nedges;++e)
                {
                    auto [a,b] = Graph.get_edge(e);
                    if(distance[a]<INF && distance[b]<INF)
                    {
                        weight_ex[e]       -= distance[b]-distance[a];
                        weight_ex[dual(e)] += distance[b]-distance[a];
                    }
                }
                
                int k = base_network::send(
                    Source,Dest,Flow_demand,
                    // admissibility
                    [weight_ex](int e)
                    {
                        return weight_ex[e]==0;
                    });
                
                if(k==0)
                    break;
                
                sent += k;
                Flow_demand -= k;
            }
            return sent;
        }
        
        public:
        network_flow_cost_PrimalDual(const digraph& in_graph):
            base_network{in_graph}    
        {}
        
        int send(
            const int Source, const int Dest,
            const std::vector<int> weight,
            int Flow_demand = INF)
        {
            if(weight.size()!=nedges)
                throw std::runtime_error(
                    "send: weight.size() != "
                    "number of edges");
            return primal_dual(Source,Dest,weight,Flow_demand); 
        }
    };

}
