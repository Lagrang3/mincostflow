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
    
    
    // TODO: 
    template<typename path_optimizer_type, typename maxflow_type>
    class mincostflow_PrimalDual : 
        public maxflow_type
    {
        using maxflow_type::nedges;
        using maxflow_type::dual;
        using maxflow_type::Graph;
        using maxflow_type::residual_cap;
    
        int execute(
            const int Source, const int Dest,
            const std::vector<int> weight)
        {   
            std::vector<int> weight_ex(nedges*2);
            
            for(int e=0;e<nedges;++e)
            {
                weight_ex.at(e) = weight.at(e);
                weight_ex.at(dual(e)) = -weight.at(e);
            } 
            
            int sent =0 ;
            path_optimizer_type path_opt(Graph);
            
            // int cycle=0;
            while(true)
            {
                bool found = path_opt.solve(
                    Source,Dest,
                    weight_ex,
                    // edge is valid if
                    [this](int e) -> bool
                    {
                        return residual_cap.at(e)>0;
                    });
                    
                if(!found)
                    break;
                    
                const auto& distance{path_opt.distance};
                
                for(int e = 0 ;e<nedges;++e)
                {
                    auto [a,b] = Graph.get_edge(e);
                    if(distance[a]<INF && distance[b]<INF)
                    {
                        weight_ex[e]       += distance[a]-distance[b];
                        weight_ex[dual(e)] -= distance[a]-distance[b];
                    }
                }
                
                
                int F = maxflow_type::solve(
                    Source,Dest,
                    // admissibility
                    [weight_ex](int e)
                    {
                        return weight_ex[e]==0;
                    });
                
                sent += F;
            }
            return sent;
        }
        
        public:
        mincostflow_PrimalDual(const digraph& in_graph):
            maxflow_type{in_graph}    
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

}
