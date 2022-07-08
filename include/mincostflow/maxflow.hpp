#pragma once

#include <mincostflow/shortestpath.hpp>
    
namespace ln
{
    template<typename T, typename path_solver_type>
    class maxflow_augmenting_path : public digraph_types
    {
        public:
        
        using value_type = T;    
        static constexpr value_type INFINITY = std::numeric_limits<value_type>::max();
        
        
        template<typename graph_t>
        value_type flow_at(
            const graph_t& g,
            const arc_pos_t e,
            const std::vector<value_type>& capacity)
        {
            auto e2 = g.arc_dual(e);
            return capacity.at(e2.x);
        }
        
        template<typename graph_t, typename condition_t>
        value_type solve(
            const graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            std::vector<value_type>& capacity,
            condition_t valid_arc)
        {
            value_type sent=0;
            path_solver_type path_solver;
            
            while(1)
            {
                bool found = path_solver.solve(
                    g,
                    Source,Dest,
                    [valid_arc,&capacity](arc_pos_t e)
                    {
                        return capacity.at(e)>0 && valid_arc(e);
                    });
                
                
                if(!found)
                    break;
                
                auto path = path_solver.get_path(g,Dest);
                
                value_type k = INFINITY;
                for(auto e : path)
                {
                    k = std::min(k,capacity.at(e));
                }
                
                for(auto e: path)
                {
                    capacity.at(e) -= k;
                    capacity.at(g.arc_dual(e)) += k;
                } 
                
                sent += k;
            }
            return sent;
        }
        
        maxflow_augmenting_path()
        {}
    };
   
    template<typename path_solver_type>
    class maxflow_scaling : public network_flow_solver
    {
        
        template<class condition_t>
        int execute(
            const int Source, const int Dest,
            condition_t admissible)
        // augmenting path
        {
            int sent=0;
            path_solver_type search_algo(Graph);
            
            int cap_flow = 1;
            for(int e : Graph.out_edges(Source))
                cap_flow = std::max(cap_flow,residual_cap.at(e));
            
            cap_flow = lower_bound_power2(cap_flow);
            
            // int cycle=0;
            for(;cap_flow>0;)
            {
                // cycle++;
                // std::cerr << "augmenting path cycle: " << cycle << '\n';
                // std::cerr << "flow sent: " << sent << '\n';
                // std::cerr << "cap flow: " << cap_flow << '\n';
            
                bool found = search_algo.solve(
                    Source,Dest,
                    // edge is valid if
                    [this,admissible,cap_flow](int e)
                    {
                        return residual_cap.at(e)>=cap_flow && admissible(e);
                    });
                
                if(! found)
                {
                    cap_flow/=2;
                    // std::cerr << "path not found!\n";
                    search_algo.reset();
                    continue;
                }
                
                auto path = search_algo.get_path(Dest);
                
                // std::cerr << "path found!\n";
                
                for(auto e: path)
                {
                    residual_cap[e] -= cap_flow;
                    residual_cap[dual(e)] += cap_flow;
                } 
                
                sent += cap_flow;
            }
            return sent;
        }
        public:
        maxflow_scaling(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int solve(
            const int Source, const int Dest,
            condition_t admissible)
        {
            return execute(Source,Dest,admissible);
        }
        int solve(
            const int Source, const int Dest)
        {
            return execute(Source,Dest,[](int){return true;});
        }
        
    };
    
    class maxflow_preflow : public network_flow_solver
    {
        
        std::vector<int> distance;
        std::vector<int> excess;
        
        template<class condition_t>
        void initialize_distance(
            const int Dest,
            condition_t valid_edge)
        {
            std::fill(distance.begin(),distance.end(),INF);
            distance.at(Dest)=0;
            
            std::queue<int> q;
            q.push(Dest);
            
            while(!q.empty())
            {
                auto n = q.front();
                q.pop();
                
                for(int e: Graph.in_edges(n))
                if( valid_edge(e) ) 
                {
                    // assert b==n
                    auto [a,b] = Graph.get_edge(e);
                    int dnew = distance[b] + 1;
                    
                    if(distance[a]==INF)
                    {
                        distance[a] = dnew;
                        q.push(a);
                    }
                }
            }
        }
        
        template<class condition_t>
        int execute(
            const int Source, const int Dest,
            condition_t valid_edge)
        {
            std::fill(excess.begin(),excess.end(),0);
            
            initialize_distance(Dest,valid_edge);
            std::queue<int> q;
            
            auto push = [&](int e)
            {
                auto [a,b] = Graph.get_edge(e);
                const int delta = std::min(excess[a],residual_cap.at(e));
                residual_cap.at(e) -= delta;
                residual_cap.at(dual(e)) += delta;
                
                assert(delta>=0);
                
                excess.at(a) -= delta;
                excess.at(b) += delta;
                
                if(delta>0 && excess.at(b)==delta)
                    q.push(b);
            };
            
            auto relabel = [&](int v)
            {
                int hmin = INF;
                for(int e : Graph.out_edges(v))
                    if(valid_edge(e) && residual_cap.at(e)>0)
                        hmin = std::min(hmin,distance.at(Graph.to_node(e)));
                if(hmin<INF)    
                    distance.at(v) = hmin+1;
            };
            
            auto discharge = [&](int a)
            {
                while(true)
                {
                    for(int e : Graph.out_edges(a))
                        if(valid_edge(e) && residual_cap.at(e)>0)
                        {
                            int b = Graph.to_node(e);
                            if(distance[a]== distance[b]+1)
                                push(e);
                        }
                    
                    if(excess.at(a)==0)
                        break;
                    
                    relabel(a);
                }
            };
            
            excess.at(Source) = INF;
            distance.at(Source) = Graph.n_vertex();
            
            for(int e : Graph.out_edges(Source))
                if(valid_edge(e))
                    push(e);
            
            while(!q.empty())
            {
                int a = q.front();
                q.pop();
                
                if(a!=Dest && a!=Source)
                    discharge(a);
            }
            return excess.at(Dest);
        }
        public:
        maxflow_preflow(const digraph& in_graph):
            network_flow_solver{in_graph},
            distance(Graph.n_vertex()),
            excess(Graph.n_vertex())
        {}
        
        template<class condition_t>
        int solve(
            const int Source, const int Dest,
            condition_t admissible)
        {
            return execute(Source,Dest,admissible);
        }
        int solve(
            const int Source, const int Dest)
        {
            return execute(Source,Dest,[](int){return true;});
        }
    };
}
