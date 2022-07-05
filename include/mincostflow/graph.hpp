#pragma once

#include <stdexcept>
#include <limits>
#include <cassert>
#include <cstdint>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

namespace ln
{   

    class digraph_types
    {
        public:
        
        // a lot more efficient if we can internally identify arcs and nodes by their unique
        // position in the buffer array
        typedef std::size_t pos_type;
        static constexpr pos_type NONE = std::numeric_limits<pos_type>::max();
        
        struct node_pos_t
        {
            pos_type x{NONE};
            bool operator<(const node_pos_t& that)const
            {
                return x < that.x;
            }
            bool operator==(const node_pos_t& that)const
            {
                return x==that.x;
            }
        };
        struct arc_pos_t
        {
            pos_type x{NONE};
            bool operator<(const arc_pos_t& that)const
            {
                return x < that.x;
            }
            bool operator==(const arc_pos_t& that)const
            {
                return x==that.x;
            }
        };
        
        struct arc_data_t
        {
            node_pos_t a{NONE},b{NONE};
            arc_pos_t dual{NONE};
            
            void clear()
            {
                a=b=node_pos_t{};
                dual=arc_pos_t{};
            }
            
            bool is_valid()const
            {
                return a.x!=NONE && b.x!=NONE;
            }
            void init(node_pos_t src, node_pos_t dst)
            {
                a = src;
                b = dst;
                dual = arc_pos_t{};
            }
        };
        
        struct node_data_t
        {
            bool is_valid_val{false};
            std::vector<arc_pos_t> out_arcs,in_arcs;
            
            void init()
            {
                is_valid_val = true;
                in_arcs.clear();
                out_arcs.clear();
            }
            
            void clear()
            {
                is_valid_val=false;
                out_arcs.clear();
                in_arcs.clear();
            }
            bool is_valid()const
            {
                return is_valid_val;
            }
            
            void rm_arc(arc_pos_t arc)
            {
                auto rm_arc_vec=[arc](std::vector<arc_pos_t>& V)
                {
                    if (auto ptr = std::find(V.begin(),V.end(),arc);
                             ptr!=V.end())
                    {
                        *ptr = V.back();
                        V.pop_back();
                    }
                };
                rm_arc_vec(in_arcs);
                rm_arc_vec(out_arcs);
            }
            void add_in_arc(arc_pos_t arc)
            {
                in_arcs.push_back(arc);
            }
            void add_out_arc(arc_pos_t arc)
            {
                out_arcs.push_back(arc);
            }
        };
    };
    
    // TODO: template on custom allocator
    template<typename node_id_t, typename arc_id_t>
    class digraph : public digraph_types
    /*
        Represents: a directed graph with dual arcs to simulate the residual network.
        This data structure represents only the topological information.
        Nodes and arcs have fixed positions.
    */
    /*
        ideally I would like to:
        
        digraph<node_id_t,arc_id_t> g;
        
        g.add_arc(node_a,node_b,arc_ab); // adds also the dual
        
        g::node_id_t n_1 = g.source_node(arc_x); 
        g::node_id_t n_2 = g.dest_node(arc_x); 
        
        g::node_id_t n = g.add_node(node_factory); // generates a new node
        
        
        
        // an optimized interface allows to pass additional structure as arrays
        g.max_nodes(); // size of the node array
        g.max_arcs() ; // size of arc array
        
        for()
    */
    {
        std::vector< arc_data_t > arcs_vec;
        std::set< arc_pos_t > free_arcs;
        std::unordered_map<arc_id_t,arc_pos_t> arcs_htable;
        std::vector<arc_id_t> arcs_ids;
        
        
        std::vector< node_data_t > nodes_vec;
        std::set< node_pos_t > free_nodes;
        std::unordered_map<node_id_t,node_pos_t> nodes_htable;
        std::vector<node_id_t> nodes_ids;
            
        void free_arc_space()
        {
            while(!arcs_vec.empty() && !arcs_vec.back().is_valid())
            // eliminate unused arcs from the back of the buffer
            {
                auto arc = arc_pos_t{arcs_vec.size()-1};
                assert(free_arcs.find(arc)!=free_arcs.end());
                
                free_arcs.erase(arc);
                arcs_vec.pop_back();
                arcs_ids.pop_back();
            }
        }
        void free_node_space()
        {
            while(!nodes_vec.empty() && !nodes_vec.back().is_valid())
            // eliminate unused nodes from the back of the buffer
            {
                auto node = node_pos_t{nodes_vec.size()-1};
                assert(free_nodes.find(node)!=free_nodes.end());
                
                free_nodes.erase(node);
                nodes_vec.pop_back();
                nodes_ids.pop_back();
            }
        }
        
        
        public:
        bool is_valid(arc_pos_t arc)const
        {
            return arc.x!=NONE && arc.x<arcs_vec.size() && arcs_vec.at(arc.x).is_valid();
        }
        bool is_valid(node_pos_t node)const
        {
            return node.x!=NONE && node.x<nodes_vec.size() && nodes_vec.at(node.x).is_valid();
        }
        
        auto arc_ends(arc_pos_t arc)const
        {
            return std::pair<node_pos_t,node_pos_t>{
                arcs_vec.at(arc.x).a,
                arcs_vec.at(arc.x).b
                };
        }
        arc_pos_t arc_dual(arc_pos_t arc)const
        {
            return arcs_vec.at(arc.x).dual;
        }
        
        void rm(arc_pos_t arc)
        {
            if(! is_valid(arc))
                return;
            
            auto [a,b] = arc_ends(arc);
            
            nodes_vec.at(a.x).rm_arc(arc);
            nodes_vec.at(b.x).rm_arc(arc);
            
            auto id = arcs_ids.at(arc.x);
            arcs_htable.erase(id);
            arcs_vec.at(arc.x).clear();
            free_arcs.insert(arc);
            
            // done
            free_arc_space();
        }
        void rm(node_pos_t node)
        {
            if(! is_valid(node))
                return;
            std::vector<arc_pos_t> ls_arcs;
            std::copy(nodes_vec.at(node.x).in_arcs.begin(),
                      nodes_vec.at(node.x).in_arcs.end(),
                      std::back_inserter(ls_arcs));
            std::copy(nodes_vec.at(node.x).out_arcs.begin(),
                      nodes_vec.at(node.x).out_arcs.end(),
                      std::back_inserter(ls_arcs));
                      
            // first remove all incoming and outgoin arcs
            for(auto arc: ls_arcs)
                rm(arc);
            
            auto id = nodes_ids.at(node.x);
            nodes_htable.erase(id);
            nodes_vec.at(node.x).clear();
            free_nodes.insert(node);
            
            // done
            free_node_space();
        }
        node_pos_t new_node()
        {
            node_pos_t node{NONE};
            if(!free_nodes.empty())
            // take a node from a free slot
            {
                node = *free_nodes.begin();
                free_nodes.erase(node);
            }else
            // no free slots,
            // take a node from the end of the buffer
            {
                node = node_pos_t{nodes_vec.size()};
                nodes_vec.emplace_back();
                nodes_ids.emplace_back();
            }
            nodes_vec.at(node.x).init();
            return node;
        }
        
        const auto& out_arcs(node_pos_t node)const
        {
            if(!is_valid(node))
                throw std::runtime_error(
                    "digraph::out_arcs invalid node");
            return nodes_vec.at(node.x).out_arcs;
        }
        
        arc_pos_t new_arc(node_pos_t a, node_pos_t b)
        {
            if(!is_valid(a) || !is_valid(b))
                throw std::runtime_error("digraph::new_arc add a new arc with invalid end nodes");
            
            arc_pos_t arc{NONE};
            if(!free_arcs.empty())
            // take an arc from a free slot
            {
                arc = *free_arcs.begin();
                free_arcs.erase(arc);
            }else
            // no free slots,
            // take an arc from the end of the buffer
            {
                arc = arc_pos_t{arcs_vec.size()};
                arcs_vec.emplace_back();
                arcs_ids.emplace_back();
            }
            arcs_vec.at(arc.x).init(a,b);
            nodes_vec.at(a.x).add_out_arc(arc); 
            nodes_vec.at(b.x).add_in_arc(arc); 
            return arc;
        }
        void set_dual(arc_pos_t arc1, arc_pos_t arc2)
        {
            if(!is_valid(arc1) || !is_valid(arc2))
                throw std::runtime_error("digraph::set_dual invalid arcs");
                
            arcs_vec.at(arc1.x).dual = arc2;
            arcs_vec.at(arc2.x).dual = arc1;
        }
        
        auto max_num_arcs()const
        {
            assert(arcs_vec.size()==arcs_ids.size());
            return arcs_vec.size();
        }
        auto max_num_nodes()const
        {
            assert(nodes_vec.size()==nodes_ids.size());
            return nodes_vec.size();
        }
        
        // translation 
        node_pos_t get_node(node_id_t id)const
        {
            node_pos_t node{NONE};
            if(auto ptr = nodes_htable.find(id);
               ptr != nodes_htable.end())
            {
                node = ptr->second;
            }
            return node;
        }
        arc_pos_t get_arc(arc_id_t id)const
        {
            arc_pos_t arc{NONE};
            if(auto ptr = arcs_htable.find(id);
               ptr != arcs_htable.end())
            {
                arc = ptr->second;
            }
            return arc;
        }
        node_pos_t add_node(node_id_t id)
        {
            auto node = get_node(id);
            if(!is_valid(node))
            {
                node = new_node();
                nodes_ids.at(node.x) = id;
                nodes_htable[id] = node;
            }
            return node;
        }
        std::pair<arc_pos_t,arc_pos_t> add_arc(node_id_t a, node_id_t b, arc_id_t id)
        {
            auto n_a = add_node(a);
            auto n_b = add_node(b);
            
            if(auto arc = get_arc(id); is_valid(arc))
                throw std::runtime_error("digraph::add_arc arc id already exists");
            
            auto arc1 = new_arc(n_a,n_b);
            auto arc2 = new_arc(n_b,n_a);
            set_dual(arc1,arc2);
            
            arcs_ids.at(arc1.x) = id;
            arcs_ids.at(arc2.x) = id;
            
            arcs_htable[id] = arc1;
            
            return {arc1,arc2};
        }
        void remove_node(node_id_t id)
        {
            auto node = get_node(id);
            if(!is_valid(node))
                return;
            nodes_htable.erase(id);
            rm(node);
        }
        void remove_arc(arc_id_t id)
        {
            auto arc = get_arc(id);
            if(!is_valid(arc))
                return;
            
            auto arc2 = arc_dual(arc);
            
            arcs_htable.erase(id);
            rm(arc);
            rm(arc2);
        }
        
        digraph()
        {}
    };
}
