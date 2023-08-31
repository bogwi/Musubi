# All purpose Graph in ZIG.

## Description

An implementation of an adjacency map graph where all edges incident to a vertex are gathered into a map, with the adjacent vertex serving as a key. \
The graph can be initiated as one of the four variants, subjects of graph's type: `.directed,.undirected` and mode: `.weighted, .unweighted`.
```zig
const Graph = Musubi(VertexId, EdgeId, EdgeWt, .undirected, .unweighted);
var graph: Graph = .{};
graph.init(allocator);
defer graph.deinit();
``` 
`VertexId`, a vertex type, can be anything hash-able, like numeric types, structs, arrays, pieces of code, or anything save floats and untagged enums. \
`EdgeId`, an edge type, can be anything. \
`EdgeWt`, type of edge's weight if the graph was initiated as weighted; any numeric type. If the graph is unweighted, the type is void. 

Here is the full API at the moment:

*GENERAL API* \
A set of standard procedures typically found in graph ADTs.
```
Vertex:
    .id:                         VertexId
    init:                            void

Edge:
    .origin:                       Vertex
    .destination:                  Vertex
    .id:                           EdgeId
    .weight:                       EdgeWt
    init:                            void
    endpoints: PairV: 
                origin
                destination
    opposite:                      Vertex



init:                                void
deinit:                              void

clearAndFree:                        void
clearRetainingCapacity:              void
ensureTotalCapacity:                !void

cloneIntoSelf:                      !void
cloneIntoSelfWithAllocator:         !void
mergeIntoSelf:                      !void

makeVertex:                        Vertex
insertVertex:                     !Vertex
insertVertexIfVertex:               !void
removeVertex:                        bool
gotVertex:                           bool
vertexCount:                          u64
vertices:                       ?[]Vertex
adjacentVertices:               ?[]Vertex
verticesIntoSet:             !AllVertices
    AllVertices:
        .vertices:           ArrayHashMap
        deinit:                      void
        list:                    []Vertex
        count:                      usize
        gotVertex:                   bool
        deleteVertex:                bool

makeEdge:                            Edge
insertEdge:                         !Edge
insertEdgeIfEdge:                   !void
removeEdge:                          bool
gotEdge:                             bool
gotEdgeIfEdge:                       bool
getEdge:                            ?Edge
edgeCount:                          usize
degree:                             usize
incidentEdges:                    ?[]Edge
edgesIntoSet:                    AllEdges
    AllEdges:
        .edges:              ArrayHashMap
        deinit:                      void
        list:                      []Edge
        count:                      usize
        gotEdge:                     bool
        deleteEdge:                  bool

```
*SPECIAL API* 

*tree Traversing*
```
traverseTree:                  !ArrayList
traverseTreeIfTarget:          !ArrayList
```
The tree traversing procedure supports four algorithms by passing a corresponding enum to the above function:
>`TreeTraverseAlg:`

>`.bfs`, breadth-first, iterative \
`.pre`, preorder, recursive \
`.post`, postorder, recursive \
`.ino`, inorder, recursive 

*graph Traversing*
```
connectionTree:              !Connections
connectionTreeExcept:        !Connections
connectionTreeThrough:       !Connections

connectionTreeIfTarget:       !Connection
connectionTreeIfTargetExcept: !Connection
connectionTreeIfTargetThrough:!Connection
```
The graph traversing procedure supports four + 1 algorithm by passing a corresponding enum to the above functions: 
>`SearchAlg:`

>`.bfs`, breadth-first search, iterative \
`.dfsA`, depth-first search, iterative \
`.dfsB`, depth-first search, iterative, true recursion emulation \
`.dfsC`, depth-first search, pure recursive \
`.dij` , Dijkstra shortest path, iterative. 

Which algorithm is better? It depends. `.bfs` and `.dij` are both Shortest-Path algorithms. The only difference between them is that `.bfs` gives the shortest path based on how many edges it needs to travel to reach the goal, whereas `.dij` considers edges' weights, treating them akin to distances connecting vertices. If all edges have the same weight, then `.dij` will produce the same result as `.bfsd`. \
The depth-first search group of algorithms is different and graph-dependent. If the graph is undirected, where all the vertices are connected randomly, they will not necessarily produce the shortest paths from origin to destination. They explore the graph as a whole and are useful in finding the longest paths possible. If we have such an undirected tangly graph of 1M vertices connected at random, and if it is possible to travel from the first vertex to the last vertex and visit all the nodes, the recursive `.dfsC` algorithm will find that path of 1M - 1 vertex. \
`.dfsB` is the author's iterative algorithm, which emulates true recursion to a great extent. In some scenarios, the paths it produces are identical to true recursion with an identical stack trace, yet in branches, it might differ. It is designed only for undirected graphs as `dfsC` substitute. \
`.dfsA` is a lazy iterative algorithm often found in books and used worldwide. It is an inversion of `.bfs` where the queue is replaced for the stack. In the case of an undirected randomly connected graph, the paths it produces will be much shorter than those of the recursive `.dfsC`. 

Additional parameters are: \
`knockout`, a set of vertices that should be removed from traversing, or traversing exclusively through them \
`target`, the target of traversing, will stop once the target is reached \
`depth`, the depth of the traversing, which, depending on the algorithm, has slightly different objectives.

The traversing process computes a connection tree from the given origin vertex to all the other vertices in the map. \
The connection tree has its own documented API to work with the result:
```
Connection
    .found:                          bool
    .explore:                 Connections
    deinit:                          void

Connections:
    .origin:                       Vertex
    .path:                   ArrayHashMap
    .discovered:             ArrayHashMap
    .last_lookup:                  Vertex
    deinit:                          void
    connectedTo:                     bool
    getAllConnected:             []Vertex  
    getDistanceTo:                 EdgeWt          
    getPathTo:                  ![]Vertex
    walkPathTo:                 !WalkPath
        WalkPath:
            .cnt:            *Connections
            .idx:                     u64
            next:                 ?Vertex
            reset:                   void
    popPathTo:                   !PopPath
        PopPath:
            .cnt:            *Connections
            .dest:                 Vertex
            next:                 ?Vertex
            reset:                   void
```

*Common-problems algorithms and their APIs* \
*topological sort*
```
topologicalSort:                     TOPO
    TOPO:
        .topo:                  ArrayList
        .acyclic:                    bool
        getAll:                 ?[]Vertex
        getPositions:             ?Vertex
        getFirst:                 ?Vertex
        getLast:                  ?Vertex
        walk:                    WalkTopo   
```
*minimum spanning tree*
```  
primJarnikMST:                        MST
kruskalMST:                           MST
    MST
        .cost:                        u64
        .tree:               ArrayHashMap          
        .len:                       usize
        getEdges:                  []Edge
        getVertexPairs:           []PairV
        gotVertexPair:               bool  
```
## Performance
The backing ADT of Musubi is Zig's superior ArrayHashMap, which has an unmatched iteration speed over keys and values and can extract keys and values for granted. That speeds up the graph's routine considerably. Calling, say, `vertices()`, you get an array of all vertices in the graph without harvesting them all into a container and only then returning them to the user. Same goes with finding `incidentEdges()` of a vertex or its `adjacentVertices()`.

### Testing 
Apple M1 laptop with 32GB of RAM, \
ReleaseFast optimization

#### Complete Binary Tree
```
20M vertices: u64
20M-1 edges: void
creation:            time: 7.667

BFS                  time: 2.682
PRE                  time: 3.020
POST                 time: 3.033
INO                  time: 3.019
```
Although not as advertised, Musubi remembers the insertion order and thus can be easily used as a general or binary tree for your projects. The only implication is that broken links must be restored manually in case of vertex or edge removal. The graph is not a linked tree and cannot behave as such. Nerveless tree traversing is implemented for directed graphs, and it is pretty fast.

#### Undirected, weighted, randomly connected, cobweb-looking graph

```
25k vertices: u64 
500k edges:    u1
creation:            time: 0.105 sec

Tree - connection tree
Paths - origin -> others         25k

BFS Tree                 time: 0.019
BFS Paths                time: 0.001
DFS A Tree               time: 0.021
DFS A Paths              time: 0.342 a
DFS B Tree               time: 0.031
DFS B Paths              time: 5.056 a
DFS C Tree               time: 0.020
DFS C Paths              time: 6.181 a
DIJ Tree                 time: 0.027
DIJ Paths                time: 0.002

MST:
Prim-Jarnik: cost: 29751 time: 0.062, 
throughput: 8.089

Kruskal:     cost: 29751 time: 0.082, 
throughput: 6.108
```
(a) Constructing all 25k-1 paths computed by depth-first algorithms happens to be a costly task. As mentioned, dfs algorithms on undirected randomly built graphs tend to produce the longest paths possible, with dfsC as a true recursive algorithm producing the longest paths. Thus, in the further results, *Paths* test will be omitted. However, such graphs are not real-life scenarios but mere benchmarking vessels. It also does not mean that DFS traversing should not be used at all for finding a connection between two points of interest if you work with such a tangled graph.

```
50k vertices: u64 
1m edges:      u1
creation:            time: 0.313 sec

Tree - connection tree
Paths - origin -> others         50k 

BFS Tree                 time: 0.057
BFS Paths                time: 0.004
DFS A Tree               time: 0.057
DFS A Paths                
DFS B Tree               time: 0.086
DFS B Paths              
DFS C Tree               time: 0.056 a
DFS C Paths              
DIJ Tree                 time: 0.103
DIJ Paths                time: 0.005

MST:
Prim-Jarnik: cost: 59264 time: 0.146 
throughput: 6.830

Kruskal:     cost: 59264 time: 0.226 
throughput: 4.418
```
(a) By an experiment, it was found that recursive dfsC breaks to segmentation fault at around 1_200_000 edges given the graph described above; therefore no data for this algorithm implementation is for larger graphs. For small undirected random graphs < 1.2M edges, using a pure recursive dfsC algorithm should be fine. 

```
100k vertices: u64 
2M edges:       u1
creation:            time: 0.785 sec

Tree - connection tree
Paths - origin -> others        100k

BFS Tree                 time: 0.131
BFS Paths                time: 0.009
DFS A Tree               time: 0.129
DFS A Paths              
DFS B Tree               time: 0.199
DFS B Paths              
DFS C Tree               
DFS C Paths              
DIJ Tree                 time: 0.246
DIJ Paths                time: 0.015

MST:
Prim-Jarnik: 
cost: 118512             time: 0.354 
throughput: 5.650

Kruskal:
cost: 118512             time: 0.503 
throughput: 3.978
```

```
1M vertices: u64 
20M edges:    u1
creation:           time: 13.289 sec

Tree - connection tree
Paths - origin -> others          1M

BFS Tree                 time: 3.213
BFS Paths                time: 0.177
DFS A Tree               time: 3.221
DFS A Paths              
DFS B Tree               time: 3.550
DFS B Paths              
DFS C Tree               
DFS C Paths              
DIJ Tree                 time: 6.335
DIJ Paths                time: 0.347

MST:
Prim-Jarnik: 
cost: 1184658            time: 7.624
throughput: 2.623

Kruskal:
cost: 1184658            time: 9.918
throughput: 2.017
```

#### Directed, weighted, acyclic, randomly connected graph
```
1M vertices: u64 
20M+ edges:  u64
creation:            time: 7.083 sec

Tree - connection tree
Paths - origin -> others          1M

BFS Tree                 time: 0.483
BFS Paths                time: 0.105
DFS A Tree               time: 0.446
DFS A Paths              time: 0.138
DFS B Tree            not applicable
DFS B Paths           not applicable   
DFS C Tree               time: 0.469
DFS C Paths              time: 0.144
DIJ Tree                 time: 1.582
DIJ Paths                time: 0.239

Topological Sort         time: 1.684
```
```
5M vertices: u64 
102M+ edges: u64
creation:           time: 46.741 sec

Tree - connection tree
Paths - origin -> others          5M

BFS Tree                 time: 4.741
BFS Paths                time: 0.765
DFS A Tree               time: 3.743
DFS A Paths              time: 0.914
DFS B Tree            not applicable  
DFS B Paths           not applicable  
DFS C Tree               time: 4.026
DFS C Paths              time: 0.907
DIJ Tree                 time: 15.195
DIJ Paths                time: 2.739

Topological Sort         time: 20.585
```
In the case of a directed graph, the results are very different. The cost of obtaining every path from the origin to all other vertices is very modest. Since there are no cycles, the recursive .dfsC algorithm inspecting 102M edges works correctly and does not break.


