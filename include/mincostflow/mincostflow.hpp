#pragma once

#include <mincostflow/maxflow.hpp>

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
    
    template<typename path_optimizer_type>
    class mincostflow_EdmondsKarp : public network_flow_solver
    {
        int execute(
            const int Source, const int Dest,
            const std::vector<int> weight)
        // augmenting path
        {   
            std::vector<int> weight_ex(nedges*2);
            for(int e=0;e<nedges;++e)
            {
                weight_ex.at(e) = weight.at(e);
                weight_ex.at(dual(e)) = -weight.at(e);
            } 
            int sent =0 ;
            path_optimizer_type path_opt(Graph);
            
            while(true)
            {
                bool found = path_opt.solve(
                    Source,Dest,
                    weight_ex,
                    // edge is valid if
                    [this](int e){
                        return residual_cap.at(e)>0;
                    });
                
                if(!found)
                    break;
                
                auto path = path_opt.get_path(Dest);
                
                int k = INF;
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
            }
            return sent;
        }
        
        public:
        mincostflow_EdmondsKarp(const digraph& in_graph):
            network_flow_solver{in_graph}    
        {}
    
        int solve(
            const int Source, const int Dest,
            const std::vector<int> weight)
        {
            if(weight.size()!=nedges)
                throw std::runtime_error(
                    "send: weight.size() != "
                    "number of edges");
            return execute(Source,Dest,weight); 
        }
    };
    
    
    // // TODO: 
    // class network_flow_cost_PrimalDual : 
    //     public network_flow_cost_solver<network_flow_PushRelabel>
    // {
    //     using base_network = network_flow_cost_solver<network_flow_PushRelabel>;
    //     
    //     int primal_dual(
    //         const int Source, const int Dest,
    //         const std::vector<int> weight,
    //         int Flow_demand = INF)
    //     // Primal-Dual
    //     {   
    //         std::vector<int> weight_ex(nedges*2);
    //         std::vector<int> potential(Graph.n_vertex(),0);
    //         
    //         for(int e=0;e<nedges;++e)
    //         {
    //             weight_ex.at(e) = weight.at(e);
    //             weight_ex.at(dual(e)) = -weight.at(e);
    //         } 
    //         
    //         int sent =0 ;
    //         shortest_path_dijkstra dijkstra(Graph);
    //         while(Flow_demand>0)
    //         {
    //             dijkstra(
    //                 Source, 
    //                 weight_ex,
    //                 // edge is valid if
    //                 [this](int e) -> bool
    //                 {
    //                     return residual_cap.at(e)>0;
    //                 });
    //             const auto& distance{dijkstra.distance};
    //             
    //             for(int e = 0 ;e<nedges;++e)
    //             {
    //                 auto [a,b] = Graph.get_edge(e);
    //                 if(distance[a]<INF && distance[b]<INF)
    //                 {
    //                     weight_ex[e]       -= distance[b]-distance[a];
    //                     weight_ex[dual(e)] += distance[b]-distance[a];
    //                 }
    //             }
    //             
    //             int k = base_network::send(
    //                 Source,Dest,Flow_demand,
    //                 // admissibility
    //                 [weight_ex](int e)
    //                 {
    //                     return weight_ex[e]==0;
    //                 });
    //             
    //             if(k==0)
    //                 break;
    //             
    //             sent += k;
    //             Flow_demand -= k;
    //         }
    //         return sent;
    //     }
    //     
    //     public:
    //     network_flow_cost_PrimalDual(const digraph& in_graph):
    //         base_network{in_graph}    
    //     {}
    //     
    //     int send(
    //         const int Source, const int Dest,
    //         const std::vector<int> weight,
    //         int Flow_demand = INF)
    //     {
    //         if(weight.size()!=nedges)
    //             throw std::runtime_error(
    //                 "send: weight.size() != "
    //                 "number of edges");
    //         return primal_dual(Source,Dest,weight,Flow_demand); 
    //     }
    // };

}
