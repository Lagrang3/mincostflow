#pragma once

#include <mincostflow/maxflow.hpp>

namespace ln
{
    template<typename T, typename path_optimizer_type>
    class mincostflow_EdmondsKarp : public maxflow_base<T>
    {
        public:
        using base_type = maxflow_base<T>;
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
    
    
    template<typename path_optimizer_type, typename maxflow_type>
    class mincostflow_PrimalDual : public maxflow_type
    {
        public:
        using base_type = maxflow_type;
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
        {   
            std::vector<value_type> reduced_weight = weight;
            
            value_type sent =0 ;
            path_optimizer_type path_opt;
            
            while(true)
            {
                path_opt.solve(
                    g,
                    Source,
                    reduced_weight,
                    // edge is valid if
                    [&residual_cap](arc_pos_t e) -> bool
                    {
                        return residual_cap.at(e)>0;
                    });
                    
                if(! path_opt.is_reacheable(Dest))
                    break;
                    
                const auto& distance{path_opt.distance};
                
                for(auto e : g.arcs())
                {
                
                    auto [a,b] = g.arc_ends(e);
                    if(distance[a]<INFINITY && distance[b]<INFINITY)
                    {
                        reduced_weight[e]       += distance[a]-distance[b];
                    }
                }
                
                
                auto F = base_type::solve(
                    g,
                    Source,Dest,
                    residual_cap,
                    // admissibility
                    [&reduced_weight](arc_pos_t e)->bool
                    {
                        return reduced_weight[e]==0;
                    });
                
                sent += F;
            }
            return sent;
        }
        
        mincostflow_PrimalDual()
        {}
    };

}
