#pragma once

#include <mincostflow/maxflow.hpp>
#include <mincostflow/scope_guard.hpp>

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
    
    template<typename path_optimizer_type, typename maxflow_type>
    class mincostflow_capacityScaling : public maxflow_type
    {
        public:
        using base_type = maxflow_type;
        using value_type = typename base_type::value_type;
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t  = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        template<typename graph_t>
        value_type solve(
            graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& residual_cap)
        {
        
            std::vector<value_type> reduced_weight = weight;
            value_type maxflow{0};
            
            // find the max-flow-anycost
            {
                std::vector<value_type> copy_residual_cap = residual_cap;
                maxflow = maxflow_type::solve(
                    g,Source,Dest,
                    copy_residual_cap,
                    [](arc_pos_t)->bool{return true;});
            }
            value_type cap_flow = lower_bound_power2(maxflow);
            
            std::vector<value_type> excess(g.max_num_nodes(),0);
            excess.at(Source) = maxflow;
            excess.at(Dest) = -maxflow;
            
            
            std::vector<value_type> weight_ex = weight;
            
            auto update_reduced_costs = 
                [&](const std::vector<value_type>& potential)
            {
                for(auto e : g.arcs())
                {
                    auto [src,dst] = g.arc_ends(e);
                    auto p_src = potential.at(src), p_dst = potential.at(dst);
                    
                    if(p_src<INFINITY && p_dst<INFINITY)
                        weight_ex.at(e) +=  p_src - p_dst;
                }
            };
            
            auto push_flow = 
                [&](arc_pos_t e,value_type delta)
            {
                    // std::cerr << " push flow at " << e << " delta = " << delta << "\n";
                    auto [src,dst] = g.arc_ends(e);
                    
                    residual_cap[e]-=delta;
                    residual_cap[g.arc_dual(e)]+=delta;
                    
                    excess.at(src) -= delta;
                    excess.at(dst) += delta;
            };
            
            // auto report = 
            // [&]()
            // {
            //     std::cerr << "residual cap + mod. costs\n";
            //     for(auto e : g.arcs())
            //     {
            //         std::cerr << " " << e << " -> " << residual_cap[e] << " " << weight_ex[e] << "\n";
            //     }
            //     std::cerr << "potential + excess\n";
            //     for(auto v : g.nodes())
            //     {
            //         std::cerr << " " << v << " -> " << excess[v] << "\n";
            //     }
            // };
            // 
            // std::cerr << " maxflow = " << maxflow << "\n";
            
            // int cycle=0;
            for(;cap_flow>0;cap_flow/=2)
            {
                // cycle++;
                // std::cerr << "cycle " << cycle << " cap_flow = " << cap_flow << '\n';
                // report();
                
                // saturate edges with negative cost
                for(auto e : g.arcs()) 
                while(residual_cap.at(e)>=cap_flow && weight_ex.at(e)<0)
                {
                    push_flow(e,cap_flow);
                }
                
                path_optimizer_type path_opt;
                
                // build S and T
                std::set<node_pos_t> Sset,Tset;
                for(auto v : g.nodes())
                {
                    if(excess.at(v)>=cap_flow)
                        Sset.insert(v);
                    if(excess.at(v)<=-cap_flow)
                        Tset.insert(v);
                }
                
                const auto multi_source_node = g.new_node();
                excess.resize(g.max_num_nodes());
                excess.at(multi_source_node) = 0;
                const Scope_guard rm_node = [&](){ g.erase(multi_source_node);};
                
                
                
                for(auto v : Sset)
                {
                    auto arc1 = g.new_arc(multi_source_node,v);
                    auto arc2 = g.new_arc(v,multi_source_node);
                    
                    g.set_dual(arc1,arc2);
                    
                    weight_ex.resize(g.max_num_arcs());
                    residual_cap.resize(g.max_num_arcs());
                    
                    weight_ex.at(arc1) = 0;
                    residual_cap.at(arc1) = excess.at(v);
                    
                    weight_ex.at(arc2) = 0;
                    residual_cap.at(arc2) = 0;
                    
                    excess.at(multi_source_node) += excess.at(v);
                    excess.at(v) = 0;
                }
                
                const Scope_guard restore_excess = [&]()
                {
                    for(auto e : g.out_arcs(multi_source_node))
                    {
                        auto [src,dst] = g.arc_ends(e);
                        excess.at(dst) = residual_cap.at(e);
                    }
                };
                
                while(!Sset.empty() && !Tset.empty())
                { 
                    path_opt.solve(
                        g, multi_source_node,
                        weight_ex,
                        [cap_flow,&residual_cap](arc_pos_t e)->bool
                        {
                            return residual_cap.at(e)>=cap_flow;
                        }
                    );
                    
                    const auto& distance{path_opt.distance};
                    
                    auto it = std::find_if(Tset.begin(),Tset.end(),
                        [&](node_pos_t v)->bool {
                            return distance.at(v)<INFINITY;
                        });
                    
                    if(it==Tset.end())
                        break;
                    
                    auto dst = *it;
                    
                    
                    // std::cerr << " vertex distance to pivot\n";
                    // for(int v=0;v<Graph.n_vertex();++v)
                    // {
                    //     std::cerr << " " <<v <<" -> " << distance[v]<<"\n";
                    // }
                    
                    update_reduced_costs(distance);
                    
                    auto path = path_opt.get_path(g,dst);
                    for(auto e: path)
                    {
                        // auto [src,dst] = g.arc_ends(e);
                        push_flow(e,cap_flow);
                    }
                    
                    if(excess.at(dst)>-cap_flow)
                        Tset.erase(dst);
                }
            }
            
            return maxflow;
        }
        
        public:
        mincostflow_capacityScaling()
        {}
    };

}
