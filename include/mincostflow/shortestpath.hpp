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
    
    template<typename T>
    class shortestPath_FIFO : public digraph_types
    /*
        Represents: shortest path label-correcting FIFO
        Invariant:
        
        User interface: 
        Complexity: pseudo-polynomial
    */
    {
        public:
        using value_type = T;
        static constexpr value_type INFINITY = std::numeric_limits<value_type>::max();
        std::vector<value_type> distance;
        std::vector<arc_pos_t>  parent;
        
        shortestPath_FIFO()
        {}
        
        template<typename graph_t, typename condition_t>
        void solve(
            const graph_t& g,
            const node_pos_t Source,
            const std::vector<value_type>& weight,
            condition_t valid_arc)
        // shortest path FIFO
        // each call resets the state of the tree and distances
        // O( pseudo-polynomial )
        {
            if(!g.is_valid(Source))
                throw std::runtime_error(
                    "shortestPath_FIFO::solve source node is not valid");
            
            const auto num_arcs = g.max_num_arcs();
            const auto num_nodes = g.max_num_nodes();
        
            if(weight.size()<num_arcs)
                throw std::runtime_error(
                    "shortestPath_FIFO::solve weight does not map arc property");
            
            distance.resize(num_nodes);
            parent.resize(num_nodes);
            
            std::fill(distance.begin(),distance.end(),INFINITY);
            std::fill(parent.begin(),parent.end(),arc_pos_t{NONE});
            
            std::queue<node_pos_t> q;
            q.push(Source);
            distance.at(Source.x)=0;
            
            while(!q.empty())
            {
                auto node = q.front();
                q.pop();
                
                for(auto e: g.out_arcs(node))
                if( valid_arc(e) ) 
                {
                    auto [a,b] = g.arc_ends(e);
                    const value_type dnew = distance.at(a.x)+weight.at(e.x);
                    
                    if(distance.at(b.x)>dnew)
                    {
                        distance.at(b.x) = dnew;
                        parent.at(b.x) = e;
                        q.push(b);
                    }
                }
            }
        }
    };
    
    template<typename T>
    class shortestPath_BellmanFord : public digraph_types
    /*
        Represents: shortest path using Bellman-Ford
        Invariant:
        
        User interface: 
        Complexity: |V| |E|
    */
    {
        public:
        using value_type = T;
        static constexpr value_type INFINITY = std::numeric_limits<value_type>::max();
        
        std::vector<value_type> distance;
        std::vector<arc_pos_t>  parent;
        
        shortestPath_BellmanFord()
        {}
        
        template<typename graph_t, typename condition_t>
        void solve (
            const graph_t& g,
            const node_pos_t Source,
            const std::vector<value_type>& weight,
            condition_t valid_arc)
        // shortest path Bellman-Ford
        // each call resets the state of the tree and distances
        // O(|V||E|)
        {
            if(!g.is_valid(Source))
                throw std::runtime_error(
                    "shortestPath_BellmanFord::solve source node is not valid");
            
            const auto num_arcs = g.max_num_arcs();
            const auto num_nodes = g.max_num_nodes();
        
            if(weight.size()<num_arcs)
                throw std::runtime_error(
                    "shortestPath_BellmanFord::solve weight does not map arc property");
            
            distance.resize(num_nodes);
            parent.resize(num_nodes);
            
            std::fill(distance.begin(),distance.end(),INFINITY);
            std::fill(parent.begin(),parent.end(),arc_pos_t{NONE});
            
            distance.at(Source.x) = 0;
            
            // TODO: use here the right number of nodes
            for(int i=0;i<num_nodes;++i)
            {
                bool updates = false;
                for(auto e: g.arcs())
                if(valid_arc(e))
                {
                    const auto [a,b] = g.arc_ends(e);
                    if(distance.at(a.x)==INFINITY)
                        continue;
                    
                    const value_type dnew = distance[a.x]+weight.at(e.x);
                    if(distance.at(b.x)>dnew)
                    {
                        distance[b.x]=dnew;
                        parent.at(b.x) = e;
                        updates = true;
                    }
                }
                if(! updates)
                    break;
            }
            // TODO: check for negative cycles
        }
    };
    
    template<typename T>
    class shortestPath_Dijkstra : public digraph_types
    /*
        Represents: shortest path with weights using Dijkstra
        Invariant:
        
        User interface: 
        Complexity: |E| + |V| log |V|
    */
    {
        public:
        
        using value_type = T;
        static constexpr value_type INFINITY = std::numeric_limits<value_type>::max();
        std::vector<value_type> distance;
        std::vector<arc_pos_t>  parent;
        
        
        shortestPath_Dijkstra()
        {}
        
        template<typename graph_t, typename condition_t>
        void solve (
            const graph_t& g,
            const node_pos_t Source,
            const std::vector<value_type>& weight,
            condition_t valid_arc)
        // Dijkstra algorithm 
        // precondition: doesnt work with negative weights!
        // O( |E|+|V| log |V| )
        {
            if(!g.is_valid(Source))
                throw std::runtime_error(
                    "shortestPath_Dijkstra::solve source node is not valid");
            
            const auto num_arcs = g.max_num_arcs();
            const auto num_nodes = g.max_num_nodes();
        
            if(weight.size()<num_arcs)
                throw std::runtime_error(
                    "shortestPath_Dijkstra::solve weight does not map arc property");
            
            distance.resize(num_nodes);
            parent.resize(num_nodes);
            std::vector<bool> visited(num_nodes,false);
            
            std::fill(distance.begin(),distance.end(),INFINITY);
            std::fill(parent.begin(),parent.end(),arc_pos_t{NONE});
            
            distance.at(Source.x) = 0;
            std::priority_queue< 
                std::pair<value_type,node_pos_t>, 
                std::vector< std::pair<value_type,node_pos_t> >, 
                std::greater<std::pair<value_type,node_pos_t> > 
                > q;
            q.push( {0,Source} );
            
            while(!q.empty())
            {
                const auto [dist,node] = q.top();
                q.pop();
                
                if(visited.at(node.x))
                    continue;
                
                visited[node.x]=true;
                
                for(auto e: g.out_arcs(node))
                if( valid_arc(e) ) 
                {
                    auto [a,b] = g.arc_ends(e);
                    
                    if(weight.at(e.x)<0)
                        throw std::runtime_error(
                            "shortestPath_Dijkstra::solve found a negative edge");
                    
                    value_type dnew = dist + weight.at(e.x);
                    if(distance.at(b.x)>dnew)
                    {
                        distance[b.x] = dnew;
                        parent[b.x] = e;
                        q.push({dnew,b});
                    }
                }
            }
        }
    };
}
