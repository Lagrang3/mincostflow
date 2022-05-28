#pragma once

#include <mincostflow/shortestpath.hpp>
    
namespace ln
{
    class network_flow_solver
    {
        protected:
        
        const std::size_t nedges;
        digraph Graph; // extended graph
        
        std::vector<int> residual_cap;
        
        int dual(int e)const
        {
            return e<nedges ? e+nedges : e-nedges;
        }
        int prime_edge(int e)const
        {
            return e<nedges ? e : e-nedges;
        }
        
        
        public:
        network_flow_solver(const digraph& in_graph):
            nedges{in_graph.n_edges()},
            
            Graph(in_graph.n_vertex()),
            residual_cap(nedges*2,0)
        {
            for(auto [a,b]: in_graph.edges())
                Graph.add_edge(a,b);
            for(auto [a,b]: in_graph.edges())
                Graph.add_edge(b,a);
        }
        
        void set_capacity(const std::vector<int>& capacity)
        {
            if(capacity.size()!=nedges)
                throw std::runtime_error(
                    "set_capacity: capacity.size() != "
                    "number of edges");
            
            for(int e=0;e<nedges;++e)
            {
                residual_cap.at(e) = capacity.at(e);
                residual_cap.at(dual(e)) = 0;
            }
        }
        
        int capacity_at(int e)const
        {
            if(e<0 || e>=nedges)
                throw std::runtime_error("capacity_at: edge id is not valid");
                
            return residual_cap.at(e) + residual_cap.at(dual(e));
        }
        int flow_at(int e)const
        {
            if(e<0 || e>=nedges)
                throw std::runtime_error("capacity_at: edge id is not valid");
            
            return residual_cap.at(dual(e));
        }
    };
    
    template<typename path_solver_type>
    class maxflow_augmenting_path : public network_flow_solver
    {
        
        template<class condition_t>
        int execute(
            const int Source, const int Dest,
            condition_t valid_edge)
        {
            int sent=0;
            path_solver_type path_solver(Graph);
            
            //int cycle=0;
            while(1)
            {
                // cycle++;
                // std::cerr << "augmenting path cycle: " << cycle << '\n';
                // std::cerr << "flow sent: " << sent << '\n';
            
                bool found = path_solver.solve(
                    Source, Dest,
                    // edge is valid if
                    [this,valid_edge](int e)
                    {
                        return residual_cap.at(e)>0 && valid_edge(e);
                    });
                    
                if(!found)
                    break;
                
                auto path = path_solver.get_path(Dest);
                
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
        maxflow_augmenting_path(const digraph& in_graph):
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
