#include <boost/python.hpp>
#include <mincostflow/mincostflow.hpp>
#include <vector>
#include <map>

namespace {
typedef ln::maxflow_augmenting_path<ln::pathSearch_BFS> maxflow_t;
typedef ln::mincostflow_PrimalDual<ln::shortestPath_Dijkstra,maxflow_t> mincostflow_t;

class simple_mcf
{
    ln::digraph Graph;
    std::vector<int> edge_cap,edge_cost;
    const int SRC{-1},DEST{-2};
    std::map<int,int> Id;
    std::vector<int> vertex_list;
    mincostflow_t* f{nullptr};
    // 
    int get_vertex(const int n)
    {
        // std::cerr << "calling get_vertex: " << n << '\n';
        int pos=vertex_list.size();
        // 
        if(auto it=Id.find(n); it==Id.end())
        {
            // std::cerr << "new vertex\n";
            // new vertex
            Id[n] = pos = vertex_list.size();
            Graph.add_vertex();
            vertex_list.push_back(n);
        }else
        {
            pos = it->second;
        }
        return pos;
    }
    
    public:
    simple_mcf()
    {
        // std::cerr << "simple_mcf ctor\n";
        get_vertex(SRC);
        get_vertex(DEST);
    }
    void add_arc(int a,int b, int cap, int cost)
    {
        // std::cerr << "adding arc (" << a << ", " << b << ") "
        //            << "cap: " << cap << " cost: " << cost << '\n';
        a = get_vertex(a);
        b = get_vertex(b);
        // std::cerr << "arc is (" << a << ", " << b << ")\n";
        edge_cap.push_back(cap);
        edge_cost.push_back(cost);
        // 
        Graph.add_edge(a,b);
    }
    void SetNodeSupply(int n, int supply)
    {
        if(supply==0)
            return;
            
        n = get_vertex(n);
        int src = get_vertex(SRC);
        int dest = get_vertex(DEST);
        
        if(supply>0)
        {
            Graph.add_edge(src,n);
            edge_cap.push_back(supply);
            edge_cost.push_back(0);
        }else
        {
            Graph.add_edge(n,dest);
            edge_cap.push_back(0-supply);
            edge_cost.push_back(0);
        }
    }
    ~simple_mcf()
    {
        if(f!=nullptr)
            delete f;
    }
    int Solve()
    {
        if(f==nullptr)
            f = new mincostflow_t(Graph);
        f->set_capacity(edge_cap);
        int src = get_vertex(SRC),
            dest = get_vertex(DEST);
            
        return f->solve(src,dest,edge_cost);
    }
};
}

BOOST_PYTHON_MODULE(pymcf)
{
    using namespace boost::python;
    
    class_<simple_mcf>("simple_mcf",init<>())
        .def("add_arc",&simple_mcf::add_arc)
        .def("SetNodeSupply",&simple_mcf::SetNodeSupply)
        .def("Solve",&simple_mcf::Solve)
    ;   
}
