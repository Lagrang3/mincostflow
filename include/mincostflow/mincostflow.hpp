#pragma once

#include <mincostflow/maxflow.hpp>

namespace ln
{
    template<typename T, typename path_optimizer_type>
    class mincostflow_EdmondsKarp : public maxflow_type<T>
    {
        public:
        using base_type = maxflow_type<T>;
        using value_type = typename base_type::value_type;    
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        template<typename graph_t>
        value_type solve(
            const graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& residual_cap
            )
        // augmenting path
        {   
            value_type sent =0 ;
            path_optimizer_type path_opt;
            
            while(true)
            {
                path_opt.solve(
                    g,
                    Source,
                    weight,
                    // edge is valid if
                    [&residual_cap](arc_pos_t e){
                        return residual_cap.at(e)>0;
                    });
                
                if(! path_opt.is_reacheable(Dest))
                    break;
                
                auto path = path_opt.get_path(g,Dest);
                
                value_type k = INFINITY;
                for(auto e : path)
                {
                    k = std::min(k,residual_cap.at(e));
                }
                
                for(auto e: path)
                {
                    residual_cap[e] -= k;
                    residual_cap[g.arc_dual(e)] += k;
                } 
                
                sent += k;
            }
            return sent;
        }
        
        mincostflow_EdmondsKarp()
        {}
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
