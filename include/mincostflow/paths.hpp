#pragma once

#include <iostream>
#include <set>
#include <algorithm>
#include <queue>
#include <vector>
#include <cassert>

namespace ln
{
    inline long long int lower_bound_power2(long long int n)
    {
        if(n<=2) return n;
        while(n != (n & -n))
            n -= (n & -n);
        return n;
    }
    
    
    class shortest_path_tree
    /*
        Represents: ?
        Invariant: 
        -> the root is the root of the tree,
        
        -> the root has no parent edge
        
        -> a node v is in the tree if and only if parent_edge[v]>=0.
        
        -> parent_edge.size() == Graph.n_vertex()
        
        User interface: ?
    */
    {
        
        const digraph& Graph;
        int root{-1};
        std::vector<int> parent_edge;
        
        public:
        
        shortest_path_tree(const digraph& graph):
            Graph{graph},
            parent_edge(Graph.n_vertex(),-1)
        {
        }
        
        void set_root(int v)
        // O(|V|)
        {
            root = v;
            if(root>=0)
            {
                std::fill(parent_edge.begin(),parent_edge.end(),-1);
            }
        }
        
        auto get_root() const
        // O(1)
        {
            return root;
        }
       
        bool has_root()const
        // O(1)
        {
            return root>=0;
        }
        
        auto get_parent(int v)const
        // O(1)
        {
            return parent_edge.at(v);
        }
        
        bool has_parent(int v)const
        // O(1)
        {
            return parent_edge.at(v)>=0;
        }
        
        void set_parent(int b, int e)
        // O(1)
        {
            if(b==root)
                throw std::runtime_error("set_parent: root cannot have a parent edge");
            
            
            int a = Graph.from_node(e);
            if( (!has_parent(a)) && a!=root)
                throw std::runtime_error("set_parent: parent node is not connected");
            
            if(has_parent(b)) 
            // weaker but less expensive than
            // if (a is in the subtree rooted at b) throw ...;
                throw std::runtime_error("set_parent: node is already connected");
            
            parent_edge.at(b)=e;
        }
       
        auto find_path(int v) const
        // O(|V|)
        {
            if(! has_root())
                throw std::runtime_error("find_path: root is not set");
            
            std::vector<int> path;
            
            if(parent_edge.at(v)<0) 
            // node v is disconnected
                return path;
            
            if(v==root)
            // we are already there
                return path;
            
            while(v!=root)
            {
                int e = parent_edge.at(v);
                path.push_back(e);
                v = Graph.from_node(e);
            }
            return path;
        }
    };
        
    constexpr int INF = std::numeric_limits<int>::max();
    
    class shortest_path_bfs : public shortest_path_tree
    /*
        Represents: shortest path using BFS
        Invariant:
        
        User interface: 
        Complexity:
    */
    {
        std::vector<int> distance;
        public:
        
        shortest_path_bfs(const digraph& graph):
            shortest_path_tree{graph},
            distance(Graph.n_vertex());
        {}
        
        template<class condition_t>
        auto operator() (
            const int Source,
            condition_t valid_edge)
        // BFS shortest path
        // each call resets the state of the tree and distances
        // O(|E|+|V|)
        {
            std::fill(distance.begin(),distance.end(),INF);
            
            set_root(Source);
            std::queue<int> q;
            distance.at(get_root()) = 0;
            q.push(get_root());
            
            while(!q.empty())
            {
                auto a = q.front();
                q.pop();
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    // assert(distance[a]!=INF);
                    auto b = Graph.to_node(e);
                    
                    if(distance[b]==INF)
                    {
                        distance[b] = distance[a]+1;
                        set_parent(b,e);
                        q.push(b);
                    }
                }
            }
        }
    };
    
    class shortest_path_dijkstra : public shortest_path_tree
    /*
        Represents: shortest path with weights using Dijkstra
        Invariant:
        
        User interface: 
        Complexity:
    */
    {
        std::vector<int> distance;
        
        public:
        
        shortest_path_dijkstra(const digraph& graph):
            shortest_path_tree{graph},
            distance(Graph.n_vertex())
        {}
        
        template<class condition_t>
        auto operator() (
            const int Source,
            const std::vector<int>& weight,
            condition_t valid_edge)
        // Dijkstra algorithm 
        // precondition: doesnt work with negative weights!
        // O( |E|+|V| log |V| )
        {
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_dijkstra: operator() : weight.size() different"
                    " from the number of edges");
            
            std::fill(distance.begin(),distance.end(),INF);
            set_root(Source);
            
            const int r = get_root();
            distance.at(r) = 0;
            std::priority_queue< 
                std::pair<int,int>, 
                std::vector< std::pair<int,int> >, 
                std::greater<std::pair<int,int> > 
                > q;
            q.push( {0,r} );
            
            while(!q.empty())
            {
                auto [dist,a] = q.top();
                q.pop();
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    auto b = Graph.to_node(e);
                    
                    int dnew = dist + weight.at(e);
                    if(distance.at(b)>dnew)
                    {
                        distance[b] = dnew;
                        set_parent(b,e);
                        q.push({dnew,b});
                    }
                }
            }
        }
    };
    
   
    // class shortest_path_label : public shortest_path_tree
    // /*
    //     Represents: shortest path with labeling
    //     Invariant:
    //     
    //     User interface: 
    //     Complexity:
    // */
    // {
    //     const int mySource,myDest;
    //     public:
    //     
    //     shortest_path_label(const digraph& graph,int Source, int Dest):
    //         shortest_path_tree{graph},
    //         mySource{Source},
    //         myDest{Dest}
    //     {
    //         std::fill(distance.begin(),distance.end(),INF);
    //         std::fill(parent_edge.begin(),parent_edge.end(),-1);
    //         root = mySource;
    //     }
    //     template<class condition_t>
    //     void initialize (
    //         condition_t valid_edge)
    //     {
    //         std::queue<int> q;
    //         
    //         distance.at(myDest)=0;
    //         q.push(myDest);
    //         
    //         while(!q.empty())
    //         {
    //             auto n = q.front();
    //             q.pop();
    //             
    //             for(int e: Graph.in_edges(n))
    //             if( valid_edge(e) ) 
    //             {
    //                 // assert b==n
    //                 auto [a,b] = Graph.get_edge(e);
    //                 int dnew = distance[b] + 1;
    //                 
    //                 if(distance[a]==INF || distance[a]>dnew)
    //                 {
    //                     distance[a] = dnew;
    //                     q.push(a);
    //                 }
    //             }
    //         }
    //         
    //         // std::cerr << "distance vector: " ;
    //         // for(auto d : distance)
    //         // {
    //         //     std::cerr << d << ' ';
    //         // }
    //         // std::cerr << '\n';
    //     }
    //     
    //     void reset()
    //     {
    //         std::fill(distance.begin(),distance.end(),INF);
    //     }
    //     template<class condition_t>
    //     auto operator() (
    //         condition_t valid_edge)
    //     {
    //         if(distance.at(mySource)==INF)
    //             initialize(valid_edge);
    //         
    //         // std::fill(parent_edge.begin(),parent_edge.end(),-1);
    //         parent_edge.at(myDest)=-1;
    //         
    //         // int cycle = 0;
    //         for(int current = mySource;
    //             distance.at(mySource)<Graph.n_vertex() && current!=myDest;)
    //         {
    //             // cycle++;
    //             // std::cerr << "advance-relabel: " << cycle << '\n';
    //             // std::cerr << "current node: " << current << '\n';
    //             // std::cerr << "distance(current): " << distance.at(current) << '\n';
    //             // if(cycle>100) break;
    //             
    //            // advance
    //            bool found_next=false;
    //            for(int e : Graph.out_edges(current))
    //            {
    //                 int next = Graph.to_node(e);
    //                 if(valid_edge(e) && distance.at(current)==distance.at(next)+1)
    //                 {
    //                     found_next = true;
    //                     parent_edge.at(next) = e;
    //                     current = next;
    //                     break;
    //                 }
    //            }
    //            if(found_next) continue; // advance success
    //            
    //            // relabel
    //            int min_dist = Graph.n_vertex()+10;
    //            for(int e : Graph.out_edges(current))
    //            {
    //                 int next = Graph.to_node(e);
    //                 if(valid_edge(e))
    //                 {
    //                     min_dist= std::min(min_dist,distance.at(next));
    //                 }
    //            }
    //            distance.at(current) = min_dist+1;
    //            
    //            // retreat
    //            if(parent_edge.at(current)>=0)
    //            {
    //                 int e = parent_edge.at(current);
    //                 current = Graph.from_node(e);
    //            }
    //         }
    //     }
    // };
};
