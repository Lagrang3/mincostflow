// https://open.kattis.com/problems/maxflow

//#include "network.hpp"

#include <iostream>
#include <set>
#include <algorithm>
#include <queue>
#include <vector>
#include <cassert>

namespace ln
{
    static constexpr int INF = std::numeric_limits<int>::max();
    
    inline long long int lower_bound_power2(long long int n)
    {
        if(n<=2) return n;
        while(n != (n & -n))
            n -= (n & -n);
        return n;
    }
    
    class digraph
    {
        std::vector< std::pair<int,int> > Edges;
        const int N_vertex;
        
        std::vector< std::vector<int> > In;
        std::vector< std::vector<int> > Out;
        
        public:
        digraph(int N):
            N_vertex{N},
            In(N),Out(N)
        {}
        
        int add_edge(int /* from */ a, int /* to */ b)
        {
            int e = Edges.size();
            Edges.push_back({a,b});
            In.at(b).push_back(e);
            Out.at(a).push_back(e);
            
            return e;
        }
        
        int n_vertex() const 
        {
            return N_vertex;
        }
        int n_edges() const
        {
            return Edges.size();
        }
        
        const auto& edges() const
        {
            return Edges;
        }
        
        auto get_edge(int e)const
        {
            return Edges.at(e);
        }
        
        auto from_node(int e)const
        {
            return Edges.at(e).first;
        }
        auto to_node(int e)const
        {
            return Edges.at(e).second;
        }
        
        const auto& out_edges(int v) const
        {
            return Out.at(v);
        }
        const auto& in_edges(int v) const
        {
            return In.at(v);
        }
    };
    
    class shortest_path_tree
    {
        public:
        
        const digraph& Graph;
        int root{-1};
        std::vector<int> distance;
        std::vector<int> parent_edge;
        
        
        shortest_path_tree(const digraph& graph):
            Graph{graph},
            distance(Graph.n_vertex(),INF),
            parent_edge(Graph.n_vertex(),-1)
        {
        }
        
        bool has_root()const
        {
            return root>=0;
        }
       
        auto find_path(int v) const
        {
            if(! has_root())
                throw std::runtime_error("find_path: root is not set");
            
            std::vector<int> path;
            
            if(parent_edge.at(v)<0) 
                // I hope this was c++20, I could use std::format
                // throw std::runtime_error("find_path: ");
                return path;
                
            
            if(v==root)
            // we are already there
                return path;
            
            while(v!=root)
            {
                int e = parent_edge[v];
                path.push_back(e);
                v = Graph.from_node(e);
            }
            return path;
        }
    };
    
    // class shortest_path_solver
    // {
    //     protected:
    //     shortest_path_tree sp_tree;
    //     
    //     public:
    //     shortest_path_solver(const digraph& Graph):
    //         sp_tree{Graph}
    //     {}
    // };
    
    class shortest_path_bfs : public shortest_path_tree
    {
        public:
        
        shortest_path_bfs(const digraph& graph):
            shortest_path_tree{graph}
        {}
        
        template<class condition_t>
        auto operator() (
            const int Source,
            const std::vector<int>& weight,
            condition_t valid_edge)
        // BFS shortest path
        {
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_bfs: operator() : weight.size() different"
                    " from the number of edges");
            
            std::fill(distance.begin(),distance.end(),INF);
            std::fill(parent_edge.begin(),parent_edge.end(),-1);
            root = Source;
            distance.at(root) = 0;
            
            std::queue<int> q;
            q.push(root);
            
            while(!q.empty())
            {
                auto a = q.front();
                q.pop();
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    auto [a,b] = Graph.get_edge(e);
                    int dnew = distance[a] + weight[e];
                    
                    if(distance[b]==INF || distance[b]>dnew)
                    {
                        distance[b] = dnew;
                        parent_edge[b] = e;
                        q.push(b);
                    }
                }
            }
        }
    };
    
    class shortest_path_dijkstra : public shortest_path_tree
    {
        public:
        
        shortest_path_dijkstra(const digraph& graph):
            shortest_path_tree{graph}
        {}
        
        template<class condition_t>
        auto operator() (
            const int Source,
            const std::vector<int>& weight,
            condition_t valid_edge)
        {
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_bfs: operator() : weight.size() different"
                    " from the number of edges");
            
            std::fill(distance.begin(),distance.end(),INF);
            std::fill(parent_edge.begin(),parent_edge.end(),-1);
            root = Source;
            distance.at(root) = 0;
            
            std::set< std::pair<int,int> > q;
            q.insert( {0,root} );
            
            while(!q.empty())
            {
                auto a = q.begin()->second;
                q.erase(q.begin());
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    auto [a,b] = Graph.get_edge(e);
                    
                    int dnew = distance.at(a) + weight.at(e);
                    if(distance.at(b)>dnew)
                    {
                        distance[b] = dnew;
                        parent_edge[b] = e;
                        q.insert({dnew,b});
                    }
                }
            }
        }
    };
    
    class network_flow_solver
    {
        protected:
        
        const int nedges;
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
   
    class shortest_path_label : public shortest_path_tree
    {
        const int mySource,myDest;
        public:
        
        shortest_path_label(const digraph& graph,int Source, int Dest):
            shortest_path_tree{graph},
            mySource{Source},
            myDest{Dest}
        {
            std::fill(distance.begin(),distance.end(),INF);
            std::fill(parent_edge.begin(),parent_edge.end(),-1);
            root = mySource;
        }
        template<class condition_t>
        void initialize (
            condition_t valid_edge)
        {
            std::queue<int> q;
            
            distance.at(myDest)=0;
            q.push(myDest);
            
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
                    
                    if(distance[a]==INF || distance[a]>dnew)
                    {
                        distance[a] = dnew;
                        q.push(a);
                    }
                }
            }
            
            // std::cerr << "distance vector: " ;
            // for(auto d : distance)
            // {
            //     std::cerr << d << ' ';
            // }
            // std::cerr << '\n';
        }
        
        void reset()
        {
            std::fill(distance.begin(),distance.end(),INF);
        }
        template<class condition_t>
        auto operator() (
            condition_t valid_edge)
        {
            if(distance.at(mySource)==INF)
                initialize(valid_edge);
            
            // std::fill(parent_edge.begin(),parent_edge.end(),-1);
            parent_edge.at(myDest)=-1;
            
            // int cycle = 0;
            for(int current = mySource;
                distance.at(mySource)<Graph.n_vertex() && current!=myDest;)
            {
                // cycle++;
                // std::cerr << "advance-relabel: " << cycle << '\n';
                // std::cerr << "current node: " << current << '\n';
                // std::cerr << "distance(current): " << distance.at(current) << '\n';
                // if(cycle>100) break;
                
               // advance
               bool found_next=false;
               for(int e : Graph.out_edges(current))
               {
                    int next = Graph.to_node(e);
                    if(valid_edge(e) && distance.at(current)==distance.at(next)+1)
                    {
                        found_next = true;
                        parent_edge.at(next) = e;
                        current = next;
                        break;
                    }
               }
               if(found_next) continue; // advance success
               
               // relabel
               int min_dist = Graph.n_vertex()+10;
               for(int e : Graph.out_edges(current))
               {
                    int next = Graph.to_node(e);
                    if(valid_edge(e))
                    {
                        min_dist= std::min(min_dist,distance.at(next));
                    }
               }
               distance.at(current) = min_dist+1;
               
               // retreat
               if(parent_edge.at(current)>=0)
               {
                    int e = parent_edge.at(current);
                    current = Graph.from_node(e);
               }
            }
        }
    };
    class network_flow_shortest_path : public network_flow_solver
    {
        public:
        network_flow_shortest_path(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int augmenting_path(
            const int Source, const int Dest,
            condition_t admissible)
        // augmenting path
        {
            int sent=0;
            shortest_path_label search_algo(Graph,Source,Dest);
            
            //int cycle=0;
            while(1)
            {
                // cycle++;
                // std::cerr << "augmenting path cycle: " << cycle << '\n';
                // std::cerr << "flow sent: " << sent << '\n';
            
                search_algo(
                    // edge is valid if
                    [this,admissible](int e)
                    {
                        return residual_cap.at(e)>0 && admissible(e);
                    });
                
                auto path = search_algo.find_path(Dest);
                
                if(path.empty())
                    break;
                    
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
        
        template<class condition_t>
        int send(
            const int Source, const int Dest,
            condition_t admissible)
        {
            return augmenting_path(Source,Dest,admissible);
        }
        int send(
            const int Source, const int Dest)
        {
            return augmenting_path(Source,Dest,[](int e){return true;});
        }
        
    };
   
    class network_flow_shortest_path_scaling : public network_flow_solver
    {
        public:
        network_flow_shortest_path_scaling(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int augmenting_path(
            const int Source, const int Dest,
            condition_t admissible)
        // augmenting path
        {
            int sent=0;
            shortest_path_label search_algo(Graph,Source,Dest);
            
            long long int cap_flow = 1;
            for(int e : Graph.out_edges(Source))
                cap_flow += residual_cap.at(e);
            
            cap_flow = lower_bound_power2(cap_flow);
            
            // int cycle=0;
            for(;cap_flow>0;)
            {
                // cycle++;
                // std::cerr << "augmenting path cycle: " << cycle << '\n';
                // std::cerr << "flow sent: " << sent << '\n';
                // std::cerr << "cap flow: " << cap_flow << '\n';
            
                search_algo(
                    // edge is valid if
                    [this,admissible,cap_flow](int e)
                    {
                        return residual_cap.at(e)>=cap_flow && admissible(e);
                    });
                
                auto path = search_algo.find_path(Dest);
                
                if(path.empty())
                {
                    cap_flow/=2;
                    // std::cerr << "path not found!\n";
                    search_algo.reset();
                    continue;
                }
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
        
        template<class condition_t>
        int send(
            const int Source, const int Dest,
            condition_t admissible)
        {
            return augmenting_path(Source,Dest,admissible);
        }
        int send(
            const int Source, const int Dest)
        {
            return augmenting_path(Source,Dest,[](int e){return true;});
        }
        
    };
    
    class network_flow_EdmondsKarp : public network_flow_solver
    {
        public:
        network_flow_EdmondsKarp(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int augmenting_path(
            const int Source, const int Dest,
            int Flow_demand,
            condition_t admissible)
        // augmenting path
        {
            int sent=0;
            shortest_path_bfs bfs(Graph);
            std::vector<int> weight(Graph.n_edges(),1);
            
            while(Flow_demand>0)
            {
                bfs(Source,weight,
                    // edge is valid if
                    [this,admissible](int e)
                    {
                        return residual_cap.at(e)>0 && admissible(e);
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
        
        template<class condition_t>
        int send(
            const int Source, const int Dest,
            int Flow_demand, 
            condition_t admissible)
        {
            return augmenting_path(Source,Dest,Flow_demand,admissible);
        }
        int send(
            const int Source, const int Dest,
            int Flow_demand=INF)
        {
            return augmenting_path(Source,Dest,Flow_demand,[](int e){return true;});
        }
    };
    
    class network_flow_capScaling : public network_flow_solver
    {
        public:
        network_flow_capScaling(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int augmenting_path(
            const int Source, const int Dest,
            condition_t admissible)
        // augmenting path with capacity scaling
        {
            int sent=0;
            shortest_path_bfs bfs(Graph);
            std::vector<int> weight(Graph.n_edges(),1);
            
            long long int cap_flow = 1;
            for(int e : Graph.out_edges(Source))
                cap_flow += residual_cap.at(e);
            
            
            cap_flow = lower_bound_power2(cap_flow);
            
            for(;cap_flow>0;)
            {   
                bfs(Source,weight,
                    // edge is valid if
                    [this,admissible,cap_flow](int e)
                    {
                        return residual_cap.at(e)>=cap_flow && admissible(e);
                    });
                
                auto path = bfs.find_path(Dest);
                
                if(path.empty())
                {
                    cap_flow /= 2;
                    continue;
                }
                
                for(auto e: path)
                {
                    residual_cap[e] -= cap_flow;
                    residual_cap[dual(e)] += cap_flow;
                } 
                sent += cap_flow;
            }
            return sent;
        }
        
        template<class condition_t>
        int send(
            const int Source, const int Dest,
            condition_t admissible)
        {
            return augmenting_path(Source,Dest,admissible);
        }
        int send(
            const int Source, const int Dest)
        {
            return augmenting_path(Source,Dest,[](int e){return true;});
        }
    };
    
    class network_flow_PushRelabel : public network_flow_solver
    {
        public:
        network_flow_PushRelabel(const digraph& in_graph):
            network_flow_solver{in_graph}
        {}
        
        template<class condition_t>
        int push_relabel(
            const int Source, const int Dest,
            int Flow_demand,
            condition_t admissible)
        {
            std::vector<int> height(Graph.n_vertex(),0);
            std::vector<int> excess(Graph.n_vertex(),0);
            std::queue<int> q;
            
            {
                shortest_path_bfs bfs(Graph); 
                std::vector<int> w(Graph.n_edges(),1);
                bfs(Dest,w,[admissible](int e){ return admissible(e); });
                height = std::move(bfs.distance);    
            }
            
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
                    if(admissible(e) && residual_cap.at(e)>0)
                        hmin = std::min(hmin,height.at(Graph.to_node(e)));
                if(hmin<INF)    
                    height.at(v) = hmin+1;
            };
            
            auto discharge = [&](int a)
            {
                while(true)
                {
                    for(int e : Graph.out_edges(a))
                        if(admissible(e) && residual_cap.at(e)>0)
                        {
                            int b = Graph.to_node(e);
                            if(height[a]== height[b]+1)
                                push(e);
                        }
                    
                    if(excess.at(a)==0)
                        break;
                    
                    relabel(a);
                }
            };
            
            excess.at(Source) = Flow_demand;
            height.at(Source) = Graph.n_vertex();
            
            for(int e : Graph.out_edges(Source))
                if(admissible(e))
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
        
        template<class condition_t>
        int send(
            const int Source, const int Dest,
            int Flow_demand, 
            condition_t admissible)
        {
            return push_relabel(Source,Dest,Flow_demand,admissible);
        }
        int send(
            const int Source, const int Dest,
            int Flow_demand=INF)
        {
            return push_relabel(Source,Dest,Flow_demand,[](int e){return true;});
        }
    };
    
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
};
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
