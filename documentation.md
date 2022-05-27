Graph representation
===
```
    class digraph;
```

Path solvers
===

Finds a path in a directed graph with un-weighted edges.
```
    class pathSearch_BFS;
```

Finds a path in a directed graph with un-weighted edges.
This label-relabel search algorithm has meaninful state.
Figure 7.6 of Ahuja 93.
```
    class pathSearch_labeling;
```

Pseudo-polynomial generic path optimization.
```
    class shortestPath_FIFO;
```
    
Bellman-Ford path optimization
```
    class shortestPath_BellmanFord;
```

Dijkstra path optimization.
```
    class shortestPath_Dijkstra;
```

Max-flow
===

Generic augmenting path algorithm, template on the path finder algorithm.
Ahuja figure 6.12.
```
    template<typename path_solver_type>
    class maxflow_augmenting_path;
```

Capacity scaling algorithm template on the path finder.
Ahuja figure 7.3.
```
    template<typename path_solver_type>
    class maxflow_scaling;
```

Preflow-push algorithm.
Ahuja figure 7.12.
```
    class maxflow_preflow;
```

Min-Cost-Flow
===
