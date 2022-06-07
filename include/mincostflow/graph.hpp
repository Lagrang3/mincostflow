#pragma once

#include <vector>
#include <algorithm>

namespace ln
{
    class digraph
    /*
        Represents: a directed graph
        Interface: be able to construct a directed graph and retrieve
        first hand information about it,
        
        Invariant: 
    */
    {
        std::vector< std::pair<int,int> > Edges;
        std::size_t N_vertex{0};
        
        std::vector< std::vector<int> > In;
        std::vector< std::vector<int> > Out;
        
        public:
        digraph(std::size_t N):
            N_vertex{N},
            In(N),Out(N)
        {}
        digraph(){}
        
        void add_vertex()
        {
            ++N_vertex;
            In.emplace_back();
            Out.emplace_back();
        }
        
        int add_edge(int /* from */ a, int /* to */ b)
        {
            int e = Edges.size();
            Edges.push_back({a,b});
            In.at(b).push_back(e);
            Out.at(a).push_back(e);
            
            return e;
        }
        
        auto n_vertex() const 
        {
            return N_vertex;
        }
        auto n_edges() const
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

}
