#pragma once

#include <mincostflow/graph.hpp>

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
    
    
    class shortestPath_tree
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
        
        int root{-1};
        std::vector<int> parent_edge;
        
        protected:
        
        const digraph& Graph;
        
        public:
        
        shortestPath_tree(const digraph& graph):
            parent_edge(graph.n_vertex(),-1),
            Graph{graph}
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
            // assert(v>=0 && v<parent_edge.size());
            return parent_edge.at(v)>=0;
        }
        
        bool is_connected(int v)const
        // O(1)
        {
            return v==root || parent_edge.at(v)>=0;
        }
        
        void set_parent(int b, int e)
        // O(1)
        {
            if(b==root)
                throw std::runtime_error("set_parent: root cannot have a parent edge");
            
            
            int a = Graph.from_node(e);
            //assert(a>=0 && a<Graph.n_vertex());
            //assert(b>=0 && b<Graph.n_vertex());
            //assert(e>=0 && e<Graph.n_edges());
            
            if( (!has_parent(a)) && a!=root)
                throw std::runtime_error("set_parent: parent node is not connected");
            
            // if (a is in the subtree rooted at b) throw ...;
            //    throw std::runtime_error("set_parent: node is already connected");
            
            parent_edge.at(b)=e;
        }
       
        auto get_path(int v) const
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
    
    class pathSearch_BFS : public shortestPath_tree
    /*
        Represents: weigthless path using BFS
        Invariant:
        
        User interface: 
        Complexity: |E|+|V|
    */
    {
        std::vector<int> distance;
        public:
        
        pathSearch_BFS(const digraph& graph):
            shortestPath_tree{graph},
            distance(Graph.n_vertex())
        {
        }
        
        void reset()
        {
        }
        
        template<class condition_t>
        bool solve (
            int Source, int Dest,
            condition_t valid_edge)
        // each call resets the state of the tree and distances
        // O(|E|+|V|)
        {
            bool found = false;
            
            set_root(Source);
            std::fill(distance.begin(),distance.end(),INF);
            distance.at(get_root()) = 0;
            
            std::queue<int> q;
            q.push(get_root());
            
            while(!q.empty())
            {
                auto a = q.front();
                q.pop();
                
                if(a==Dest)
                {
                    found = true;
                    break;
                }
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
            return found;
        }
    };
    
    class pathSearch_labeling : public shortestPath_tree
    /*
        Represents: shortest path with labeling
        Invariant:
        
        User interface: 
        Complexity:
    */
    {
        int last_source{-1},last_dest{-1};
        std::vector<int> distance,dist_freq;
        public:
        
        pathSearch_labeling(const digraph& graph):
            shortestPath_tree{graph},
            distance(Graph.n_vertex()),
            dist_freq(Graph.n_vertex()+1)
        {
        }
        template<class condition_t>
        void initialize (
            condition_t valid_edge)
        {
            set_root(last_source);
            std::fill(distance.begin(),distance.end(),INF);
            std::fill(dist_freq.begin(),dist_freq.end(),0);
            std::queue<int> q;
            
            distance.at(last_dest)=0;
            
            // TODO: write a general purpose BFS label solver
            q.push(last_dest);
            
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
                        dist_freq.at(dnew)++;
                        q.push(a);
                    }
                }
            }
        }
        
        void reset()
        {
            last_source = last_dest = -1;
        }
        template<class condition_t>
        bool solve (
            const int Source, const int Dest,
            condition_t valid_edge)
        {
            if(last_source<0 || last_dest<0 || last_source!=Source || last_dest!=Dest)
            {
                last_source = Source;
                last_dest = Dest;
                initialize(valid_edge);
            }
            set_root(Source);
            
            for(int current = Source;
                distance.at(Source)<Graph.n_vertex() && current!=Dest;)
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
                        set_parent(next,e);
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
               {
                    const int new_dist = min_dist+1;
                    const int old_dist = distance.at(current);
                    distance.at(current) = new_dist;
                    if(new_dist<dist_freq.size())
                        dist_freq.at(new_dist)++;
                    dist_freq.at(old_dist)--;
                    if(dist_freq.at(old_dist)==0)
                        break;
               }
               
               // retreat
               if(get_parent(current)>=0)
               {
                    int e = get_parent(current);
                    current = Graph.from_node(e);
               }
            }
            return has_parent(Dest);
        }
        // TODO: optimize page 219 Ahuja
    };
    
    class shortestPath_FIFO : public shortestPath_tree
    /*
        Represents: shortest path label-correcting FIFO
        Invariant:
        
        User interface: 
        Complexity: pseudo-polynomial
    */
    {
        public:
        std::vector<int> distance;
        
        shortestPath_FIFO(const digraph& graph):
            shortestPath_tree{graph},
            distance(Graph.n_vertex())
        {}
        
        auto get_distance(int v)const
        // O(1)
        {
            if(! is_connected(v))
                throw std::runtime_error("shortest_path_bfs_weigthed::get_distance "
                "node is disconnected");
            
            return distance.at(v);
        }
        
        template<class condition_t>
        bool solve(
            const int Source, const int Dest,
            const std::vector<int>& weight,
            condition_t valid_edge)
        // shortest path FIFO
        // each call resets the state of the tree and distances
        // O( pseudo-polynomial )
        {
            bool found=false;
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_bfs_weigthed: operator() : weight.size() different"
                    " from the number of edges");
            std::fill(distance.begin(),distance.end(),INF);
            
            set_root(Source);
            std::queue<int> q;
            distance.at(get_root()) = 0;
            q.push(get_root());
            
            while(!q.empty())
            {
                auto a = q.front();
                q.pop();
                
                if(a==Dest)
                    found=true;
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    // assert(distance[a]!=INF);
                    auto b = Graph.to_node(e);
                    const int dnew = distance.at(a)+weight.at(e);
                    
                    if(distance[b]==INF || distance[b]>dnew)
                    {
                        distance[b] = dnew;
                        set_parent(b,e);
                        q.push(b);
                    }
                }
            }
            return found;
        }
    };
    
    class shortestPath_BellmanFord : public shortestPath_tree
    /*
        Represents: shortest path using Bellman-Ford
        Invariant:
        
        User interface: 
        Complexity: |V| |E|
    */
    {
        
        public:
        std::vector<int> distance;
        
        shortestPath_BellmanFord(const digraph& graph):
            shortestPath_tree{graph},
            distance(Graph.n_vertex())
        {}
        
        auto get_distance(int v)const
        // O(1)
        {
            if(! is_connected(v))
                throw std::runtime_error("shortest_path_BellmanFord::get_distance "
                "node is disconnected");
            
            return distance.at(v);
        }
        
        template<class condition_t>
        bool solve (
            const int Source, const int Dest,
            const std::vector<int>& weight,
            condition_t valid_edge)
        // shortest path Bellman-Ford
        // each call resets the state of the tree and distances
        // O(|V||E|)
        {
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_BellmanFord: operator() : weight.size() different"
                    " from the number of edges");
            std::fill(distance.begin(),distance.end(),INF);
            
            set_root(Source);
            distance.at(get_root()) = 0;
            
            for(int i=0;i<Graph.n_vertex();++i)
            {
                bool updates = false;
                for(int e=0;e<Graph.n_edges();++e)
                if(valid_edge(e))
                {
                    const auto [a,b] = Graph.get_edge(e);
                    if(distance.at(a)==INF)
                        continue;
                    
                    const auto dnew = distance[a]+weight.at(e);
                    if(distance.at(b)==INF || distance.at(b)>dnew)
                    {
                        distance[b]=dnew;
                        set_parent(b,e);
                        updates = true;
                    }
                }
                if(! updates)
                    break;
            }
            // TODO: check for negative cycles
            return distance.at(Dest)<INF;
        }
    };
    
    class shortestPath_Dijkstra : public shortestPath_tree
    /*
        Represents: shortest path with weights using Dijkstra
        Invariant:
        
        User interface: 
        Complexity: |E| + |V| log |V|
    */
    {
        
        public:
        std::vector<int> distance;
        
        shortestPath_Dijkstra(const digraph& graph):
            shortestPath_tree{graph},
            distance(Graph.n_vertex())
        {}
        
        auto get_distance(int v)const
        // O(1)
        {
            if(! is_connected(v))
                throw std::runtime_error("shortest_path_dijkstra::get_distance "
                "node is disconnected");
            
            return distance.at(v);
        }
        template<class condition_t>
        bool solve (
            const int Source, const int Dest,
            const std::vector<int>& weight,
            condition_t valid_edge,
            bool prune=true)
        // Dijkstra algorithm 
        // precondition: doesnt work with negative weights!
        // O( |E|+|V| log |V| )
        {
            if(weight.size()!=Graph.n_edges())
                throw std::runtime_error(
                    "shortest_path_dijkstra: operator() : weight.size() different"
                    " from the number of edges");
            
            std::vector<bool> visited(Graph.n_vertex(),false);
            std::fill(distance.begin(),distance.end(),INF);
            set_root(Source);
            
            const int r = get_root();
            // assert(r>=0 && r<distance.size());
            distance.at(r) = 0;
            std::priority_queue< 
                std::pair<int,int>, 
                std::vector< std::pair<int,int> >, 
                std::greater<std::pair<int,int> > 
                > q;
            q.push( {0,r} );
            
            while(!q.empty())
            {
                const auto [dist,a] = q.top();
                q.pop();
                
                if(visited.at(a))
                    continue;
                
                visited[a]=true;
                
                if(a==Dest && prune)
                    break;
                
                for(int e: Graph.out_edges(a))
                if( valid_edge(e) ) 
                {
                    // assert(e>=0 && e<Graph.n_edges());
                    auto b = Graph.to_node(e);
                    // assert(b>=0 && b<Graph.n_vertex());
                    
                    if(weight.at(e)<0)
                        throw std::runtime_error("Dijkstra found a negative edge");
                    
                    int dnew = dist + weight.at(e);
                    if(distance.at(b)>dnew)
                    {
                        distance[b] = dnew;
                        set_parent(b,e);
                        q.push({dnew,b});
                    }
                }
            }
            return visited.at(Dest);
        }
    };
    
   
}
