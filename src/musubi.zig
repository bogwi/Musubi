// Apache 2.0 (c) bogwi 2023 //

const std = @import("std");
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const eql = std.meta.eql;
const Tag = std.meta.activeTag;

const AutoContext = @import("autocontext.zig").AutoContext;
const pieQ = @import("pieq.zig").PieQ;
const Cache = @import("cache.zig").Cache;

pub const GraphType = enum { directed, undirected };
pub const GraphMode = enum { weighted, unweighted };

/// AN ADJACENCY MAP GRAPH VARIANT, where
/// all edges incident to a vertex are gathered into a map,
/// with the adjacent vertex serving as a key.
///
/// Parameters:
///
/// `vertexId`: anything hash-able can be a vertex, save float and untagged union types. \
/// `edgeId`: any type or void if edges need no ID, sign, or color attached to them. \
/// `edgeWt`: edge's weight; any numeric type, or void if the graph is `.unweighted`. \
/// `graphType`: `.directed`, or `.undirected`. \
/// `graphMode`: `.weighted`, or `.unweighted`.
pub fn Musubi(
    comptime vertexId: type,
    comptime edgeId: type,
    comptime edgeWt: type,
    comptime graphType: GraphType,
    comptime graphMode: GraphMode,
) type {
    if (graphMode == .weighted and (@typeInfo(edgeWt) != .Int and @typeInfo(edgeWt) != .Float))
        @compileError("The graphMode is .weighted but the edge's weight is not .Int nor .Float, yet: " ++ @typeName(edgeWt));
    if (graphMode == .weighted and @typeInfo(edgeWt) == .Void)
        @compileError("The graphMode is .weighted but edge's weight is " ++ @typeName(edgeWt));
    if (graphMode == .unweighted and @typeInfo(edgeWt) != .Void)
        @compileError("The graphMode is .unweighted but edge's weight is " ++ @typeName(edgeWt));
    return struct {
        pub const Vertex = struct {
            id: vertexId,

            /// Initiates a vertex with the given Id and returns it to the user.
            pub fn init(id: vertexId) Vertex {
                return Vertex{ .id = id };
            }
        };
        const PairV = struct { Vertex, Vertex };

        pub const Edge = struct {
            origin: Vertex,
            destination: Vertex,
            id: edgeId,
            weight: edgeWt,

            /// Initiates an Edge with given parameters and returns it to the user.
            pub fn init(
                origin: Vertex,
                destination: Vertex,
                id: edgeId,
                weight: edgeWt,
            ) Edge {
                return Edge{
                    .origin = origin,
                    .destination = destination,
                    .id = id,
                    .weight = weight,
                };
            }
            /// Returns a pair of vertices used to initiate this edge in the format of
            /// indexable struct `.{ origin, destination }`.
            pub fn endpoints(self: Edge) PairV {
                return .{ self.origin, self.destination };
            }
            /// Returns the opposite vertex to the given, both used to initiate this edge.
            pub fn opposite(self: Edge, v: Vertex) Vertex {
                return if (eql(v, self.origin)) self.destination else self.origin;
            }
        };

        /// Set of errors the client might encounter using the graph.
        pub const GraphError = error{
            // Graph's native errs
            NegativeWeight,
            MissingVertex,

            // Allocator's errs
            OutOfSize,
            OutOfMemory,

            // Traversing algorithm's errs
            UnknownAlgorithm,
            DijkstraAlg_With_Unweighted_Graph,
            dfsB_With_Directed_Graph,
            QueueIsEmpty,
            AttemptToModifyLockedRoot,
            InfiniteLoop,
        };

        pub const HashMap1 = std.ArrayHashMap(Vertex, HashMap2, AutoContext(Vertex), true);
        pub const HashMap2 = std.ArrayHashMap(Vertex, Edge, AutoContext(Vertex), true);

        const IteratorV = HashMap1.Iterator;
        const IteratorE = HashMap2.Iterator;

        const directed: bool = if (graphType == .directed) true else false;

        /// Map containing all outgoing vertices.
        outGoing: HashMap1 = undefined,

        /// Map containing all incoming vertices. If the graph is `.undirected`,
        /// inComing is a pointer to outGoing.
        inComing: if (directed) HashMap1 else *HashMap1 = undefined,

        /// Allocator used to initiate the graph.
        alloc: Allocator = undefined,

        /// Musubi
        pub const Self = @This();

        // MAIN API //

        /// Initiates the graph using the given allocator.
        pub fn init(self: *Self, alloc: Allocator) void {
            self.outGoing = HashMap1.init(alloc);
            self.inComing = if (directed) HashMap1.init(alloc) else &self.outGoing;
            self.alloc = alloc;
        }

        /// Clears the graph of all data and releases the backing allocation.
        pub fn deinit(self: *Self) void {
            var out = self.outGoing.iterator();
            while (out.next()) |entry| {
                entry.value_ptr.*.deinit();
            }
            self.outGoing.deinit();
            if (directed) {
                var in = self.inComing.iterator();
                while (in.next()) |entry| {
                    entry.value_ptr.*.deinit();
                }
                self.inComing.deinit();
            }
        }

        /// Clears the graph of all data.
        pub fn clearAndFree(self: *Self) void {
            for (self.outGoing.values()) |*dest_vtx| {
                dest_vtx.deinit();
            }
            self.outGoing.clearAndFree();
            if (directed) {
                for (self.inComing.values()) |*origin_vtx| {
                    origin_vtx.deinit();
                }
                self.inComing.clearAndFree();
            }
        }

        /// Clears the graph of all data but retains the backing allocation for future use.
        pub fn clearRetainingCapacity(self: *Self) void {
            for (self.outGoing.values()) |*dest_vtx| {
                dest_vtx.deinit();
            }
            self.outGoing.clearRetainingCapacity();
            if (directed) {
                for (self.inComing.values()) |*origin_vtx| {
                    origin_vtx.deinit();
                }
                self.inComing.clearRetainingCapacity();
            }
        }

        /// Adjusts the size of the graph to carry the given number of vertices
        /// without triggering new allocation. Call it whenever you know beforehand
        /// the number of vertices the graph needs to hold.
        pub fn ensureTotalCapacity(self: *Self, N: usize) GraphError!void {
            try self.outGoing.ensureTotalCapacity(@intCast(N));
            if (directed) try self.inComing.ensureTotalCapacity(@intCast(N));
        }

        /// Clones the *other* graph of the same type into self
        /// using the allocator of the *other*. The process is fast,
        /// but not entirely instant.
        /// Usage:
        ///
        /// `var new_graph: @TypeOf(existing_graph) = .{}`\
        /// `try new_graph.cloneIntoSelf(&existing_graph);`\
        /// `defer new_graph.deinit();`
        ///
        /// **@IMPORTANT:** To avoid allocation issues, `cloneIntoSelf()`
        /// **erases the content of *self*** (if was called on the existing graph),
        /// **changing it for the *other***. Therefore, the best practice for this
        /// method is as shown above.
        pub fn cloneIntoSelf(self: *Self, other: *Self) !void {
            self.deinit();

            self.outGoing = try other.outGoing.clone();
            self.inComing = if (directed) try other.inComing.clone() else &self.outGoing;
            self.alloc = other.alloc;

            var self_values = self.outGoing.values();
            var other_values = self.outGoing.values();
            for (self_values, other_values) |*sv, cv| {
                sv.* = try cv.clone();
            }

            if (directed) {
                self_values = self.inComing.values();
                other_values = self.inComing.values();
                for (self_values, other_values) |*sv, cv| {
                    sv.* = try cv.clone();
                }
            }
        }

        /// Clones the *other* graph of the same type into self
        /// using the allocator of the *other*. The process is fast,
        /// but not entirely instant.
        /// Usage:
        ///
        /// `var new_graph: @TypeOf(existing_graph)  = .{}`\
        /// `try new_graph.cloneIntoSelfWithAllocator(&existing_graph, alloc);`\
        /// `defer new_graph.deinit();`
        ///
        /// **@IMPORTANT:** To avoid allocation issues, `cloneIntoSelfWithAllocator()`
        /// **erases the content of *self*** (if was called on the existing graph),
        /// **changing it for the *other***. Therefore, the best practice for this
        /// method is as shown above.
        pub fn cloneIntoSelfWithAllocator(self: *Self, other: *Self, alloc: Allocator) !void {
            self.deinit();

            self.outGoing = try other.outGoing.cloneWithAllocator(alloc);
            self.inComing = if (directed) try other.inComing.cloneWithAllocator(alloc) else &self.outGoing;
            self.alloc = alloc;

            var self_values = self.outGoing.values();
            var other_values = self.outGoing.values();
            for (self_values, other_values) |*sv, cv| {
                sv.* = try cv.clone();
            }

            if (directed) {
                self_values = self.inComing.values();
                other_values = self.inComing.values();
                for (self_values, other_values) |*sv, cv| {
                    sv.* = try cv.clone();
                }
            }
        }

        /// Merges the *other* graph of the same type into *self*.
        /// Since the graph is a map, it will overwrite all those vertices
        /// of *self* that are present in the *other* as well as their respected edges.
        ///
        /// This process is iterative, bounding to the size of the *other*.
        pub fn mergeIntoSelf(self: *Self, other: *Self) !void {
            var other_entries = other.outGoing.iterator();

            for (other.vertices()) |vtx| {
                try self.insertVertexIfVertex(vtx);
            }

            while (other_entries.next()) |entry1| {
                var dest_entries = entry1.value_ptr.iterator();
                while (dest_entries.next()) |entry2| {
                    const edge = entry2.value_ptr.*;
                    try self.insertEdgeIfEdge(edge);
                }
            }
        }

        /// Makes a new Vertex instance outside the graph.
        /// For inserting it into the graph see `insertVertexIfVertex()`.
        pub fn makeVertex(self: *Self, id: vertexId) Vertex {
            _ = self;
            return Vertex.init(id);
        }

        /// Creates, inserts, and returns a Vertex object with the given *id* parameter.
        ///
        /// It *will* clobber an existing vertex in the graph if one contains
        /// the same single identical ID. Yes, vertices must be unique to be not clobbered.
        /// However, if your vertexId is a struct containing multiple fields, -
        /// they might overlap, of course, but as long as one vertex is
        /// not a perfect copy of another.
        /// The rule is simple: every vertex is a key, and all keys must be unique.
        pub fn insertVertex(self: *Self, id: vertexId) GraphError!Vertex {
            var vtx = Vertex.init(id);

            var m1 = HashMap2.init(self.alloc);
            try self.outGoing.put(vtx, m1);

            if (directed) {
                var m2 = HashMap2.init(self.alloc);
                try self.inComing.put(vtx, m2);
            }
            return vtx;
        }

        /// Inserts the given vertex into the graph.
        /// Clobbers any exciting if the graph contains such.
        /// It is for vertices created outside the graph with `makeVertex()`.
        pub fn insertVertexIfVertex(self: *Self, vtx: Vertex) GraphError!void {
            try self.outGoing.put(vtx, HashMap2.init(self.alloc));
            if (directed) try self.inComing.put(vtx, HashMap2.init(self.alloc));
        }

        /// Makes a new Edge instance outside the graph.
        /// For inserting it into the graph, see `insertEdgeIfEdge()`.
        /// If the graph was initiated as `.unweighted`,
        /// than the `weight` parameter must be set to void, `{}`.
        /// This method does not check attempts to make the edge
        /// with missing or previously deleted vertices.
        pub fn makeEdge(self: *Self, origin: Vertex, dest: Vertex, edge_id: edgeId, weight: edgeWt) Edge {
            if (edgeWt != void and weight < 0)
                return GraphError.NegativeWeight;

            _ = self;
            return Edge.init(origin, dest, edge_id, weight);
        }

        /// Inserts and returns an Edge with the given *id* and non-negative *weight*
        /// in between two given vertices *origin* and *destination*.
        /// If the graph was initiated as `.unweighted`, then the *weight* parameter must be set to void, `{}`.
        /// This method checks negative weight if one was passed, as well as missing vertices
        /// if previously deleted vertices have been passed, and returns runtime errors in both cases.
        pub fn insertEdge(self: *Self, origin: Vertex, dest: Vertex, edge_id: edgeId, weight: edgeWt) GraphError!Edge {
            if (edgeWt != void and weight < 0)
                return GraphError.NegativeWeight;

            var edge = Edge.init(origin, dest, edge_id, weight);

            if (self.outGoing.getPtr(origin)) |origin_| {
                try origin_.put(dest, edge);
            } else return GraphError.MissingVertex;

            if (self.inComing.getPtr(dest)) |dest_| {
                try dest_.put(origin, edge);
            } else return GraphError.MissingVertex;

            return edge;
        }

        /// Inserts the given edge into the graph. Clobbers if the graph already
        /// contains an edge in between origin and destination vertices.
        /// It is for edges created outside with `makeEdge()`. This method does not check
        /// attempts to insert an edge with missing vertices.
        pub fn insertEdgeIfEdge(self: *Self, edge: Edge) GraphError!void {
            try self.outGoing.getPtr(edge.origin).?.put(edge.destination, edge);
            try self.inComing.getPtr(edge.destination).?.put(edge.origin, edge);
        }

        /// Removes the given edge from the graph.
        /// Returns *true* upon success, and *false* otherwise.
        /// False might indicate that one of the vertices (of both)
        /// and the edge are missing or were not initialized.
        pub fn removeEdge(self: *Self, edge: Edge) bool {
            var origin_result: bool = false;
            var dest_result: bool = false;

            if (self.outGoing.getPtr(edge.origin)) |origin_vtx|
                origin_result = origin_vtx.swapRemove(edge.destination);

            if (self.inComing.getPtr(edge.destination)) |dest_vtx|
                dest_result = dest_vtx.swapRemove(edge.origin);

            return if (origin_result and dest_result) true else false;
        }

        /// Removes the given vertex from the graph.
        /// Returns *true* upon success, and *false* otherwise.
        /// False might indicate that the vertex is missing or was not initialized.
        pub fn removeVertex(self: *Self, vtx: Vertex) bool {
            if (self.outGoing.getPtr(vtx)) |adjacent| {
                for (adjacent.keys()) |vtx_| {
                    assert(self.inComing.getPtr(vtx_).?.swapRemove(vtx));
                }
                adjacent.deinit();
                assert(self.outGoing.swapRemove(vtx));

                if (directed) {
                    if (self.inComing.getPtr(vtx)) |adjacent_| {
                        for (adjacent_.keys()) |vtx_| {
                            assert(self.outGoing.getPtr(vtx_).?.swapRemove(vtx));
                        }
                        adjacent_.deinit();
                        assert(self.inComing.swapRemove(vtx));
                    }
                }
                return true;
            }
            return false;
        }

        /// Returns *true* if the graph contains the given vertex,
        /// otherwise returns *false*.
        pub fn gotVertex(self: *Self, vtx: Vertex) bool {
            return self.outGoing.contains(vtx);
        }

        /// Returns the total number of all vertices in the graph.
        pub fn vertexCount(self: *Self) usize {
            return self.outGoing.count();
        }

        /// Returns an array, with pointer and len, of all vertices present in the graph.
        pub fn vertices(self: *Self) []Vertex {
            return self.outGoing.keys();
        }

        /// Returns *true* if the graph has an edge in between vertices *origin* and *destination*.
        /// To check whether the graph contains a concrete edge, use `gotEdgeIfEdge()`.
        pub fn gotEdge(self: *Self, origin: Vertex, dest: Vertex) bool {
            if (self.outGoing.get(origin)) |out| {
                return out.contains(dest);
            } else return false;
        }

        /// Returns *true* if the graph has the given edge.
        pub fn gotEdgeIfEdge(self: *Self, edge: Edge) bool {
            return self.gotEdge(edge.origin, edge.destination);
        }

        /// Returns the edge in between *origin* and *destination*, or *null* if no such edge.
        pub fn getEdge(self: *Self, origin: Vertex, dest: Vertex) ?Edge {
            if (self.outGoing.get(origin)) |out| {
                if (out.get(dest)) |edge| {
                    return edge;
                } else return null;
            } else return null;
        }

        /// Return the total number of edges in the graph.
        pub fn edgeCount(self: *Self) usize {
            var out_entries = self.outGoing.iterator();
            var numberOfEdges: usize = 0;
            while (out_entries.next()) |entry| {
                numberOfEdges += entry.value_ptr.*.count();
            }
            return if (directed) numberOfEdges else @divFloor(numberOfEdges, 2);
        }

        /// Controller enum.
        pub const EdgeClass = enum { outgoing, incoming };
        /// Returns number of all *outgoing* edges incident to
        /// the given vertex *vtx* if `edge_class` parameter is set to `.outgoing`.
        /// If it is set to `.incoming`, the function returns the number
        /// of all *incoming* edges incident to the given vertex.
        /// Returns *null* if and only if the given vertex is missing, but not an error.
        pub fn degree(self: *Self, vtx: Vertex, edge_class: EdgeClass) ?usize {
            if (edge_class == .incoming) {
                if (self.inComing.getPtr(vtx)) |in| {
                    return in.count();
                } else return null;
            }
            if (self.outGoing.getPtr(vtx)) |out| {
                return out.count();
            } else return null;
        }

        /// Returns an array, with pointer and len, of all *outgoing* edges incident to
        /// the given vertex *vtx* if `edge_class` parameter is set to `.outgoing`.
        /// If it is set to `.incoming`, the function returns an array
        /// of all *incoming* edges incident to the given vertex.
        /// Returns *null* if and only if the given vertex is missing, but not an error.
        pub fn incidentEdges(self: *Self, vtx: Vertex, edge_class: EdgeClass) ?[]Edge {
            if (edge_class == .incoming) {
                if (self.inComing.getPtr(vtx)) |incident| {
                    return incident.values();
                } else return null;
            }
            if (self.outGoing.getPtr(vtx)) |incident| {
                return incident.values();
            } else return null;
        }
        /// Returns an array, with pointer and len, of all *outgoing* vertices adjacent to
        /// the given vertex *vtx* if `edge_class` parameter is set to `.outgoing`.
        /// If it is set to `.incoming`, the function returns an array
        /// of all *incoming* vertices adjacent to the given vertex.
        /// Returns *null* if and only if the given vertex is missing, but not an error.
        pub fn adjacentVertices(self: *Self, vtx: Vertex, edge_class: EdgeClass) ?[]Vertex {
            if (edge_class == .incoming) {
                if (self.inComing.getPtr(vtx)) |adjacent| {
                    return adjacent.keys();
                } else return null;
            }
            if (self.outGoing.getPtr(vtx)) |adjacent| {
                return adjacent.keys();
            } else return null;
        }

        // AUX API //

        const SetE = std.ArrayHashMap(PairV, Edge, AutoContext(PairV), false);

        /// Returns a struct containing a set variable `edges` with all edges
        /// present in the graph and accessed by calling `list()` method.
        /// Edges set can be iterated over. Edges set can be modified.
        ///
        /// Requires a `deinit()` call.
        pub fn edgesIntoSet(self: *Self) GraphError!AllEdges {
            var edges = SetE.init(self.alloc);
            var out_entries = self.outGoing.iterator();

            while (out_entries.next()) |entry1| {
                var dest_entries = entry1.value_ptr.iterator();
                while (dest_entries.next()) |entry2| {
                    const edge = entry2.value_ptr.*;
                    try edges.put(.{ edge.origin, edge.destination }, edge);
                }
            }
            return .{ .edges = edges };
        }
        pub const AllEdges = struct {
            edges: SetE,

            /// Clears the edges set and releases backing allocation.
            pub fn deinit(self: *AllEdges) void {
                self.edges.deinit();
            }
            /// Returns edges set as an array.
            pub fn list(self: AllEdges) []Edge {
                return self.edges.values();
            }
            /// Returns the number of edges in the edges set.
            pub fn count(self: AllEdges) usize {
                return self.edges.count();
            }
            /// Returns *true* if the edges set contains the given edge.
            pub fn gotEdge(self: AllEdges, edge: Edge) bool {
                const endpoints = edge.endpoints();
                return self.edges.contains(endpoints);
            }
            /// Removes the given edge from the edges set. Invalidates pointers for any
            /// array of edges obtained previously by calling `list()`.
            pub fn deletedEdge(self: *AllEdges, edge: Edge) bool {
                var pair_vtx = PairV{ edge.origin, edge.destination };
                return self.edges.swapRemove(pair_vtx);
            }
        };

        pub const SetV = std.ArrayHashMap(Vertex, void, AutoContext(Vertex), false);

        /// Returns a struct containing a set variable with all vertices present
        /// in the graph and accessed by calling `list()` method.
        /// Vertices set can be iterated over. Vertices set can be modified.
        ///
        /// Requires a `deinit()` call.
        pub fn verticesIntoSet(self: *Self) GraphError!AllVertices {
            var vertices_ = SetV.init(self.alloc);
            const all_vtx = self.outGoing.keys();
            for (all_vtx) |vtx| {
                try vertices_.put(vtx, {});
            }
            return .{ .vertices = vertices_ };
        }

        pub const AllVertices = struct {
            vertices: SetV,

            /// Clears the vertices set and releases backing allocation.
            pub fn deinit(self: *AllVertices) void {
                self.vertices.deinit();
            }

            /// Returns the vertices set as an array.
            pub fn list(self: AllVertices) []Vertex {
                return self.vertices.keys();
            }

            /// Returns the number of vertices in the vertices set.
            pub fn count(self: AllVertices) usize {
                return self.vertices.count();
            }

            /// Returns *true* if the vertices set contains the given vertex.
            pub fn gotVertex(self: AllVertices, vtx: Vertex) bool {
                return self.vertices.contains(vtx);
            }

            /// Removes the given vertex from the vertices set. Invalidates pointers for any
            /// array of vertices obtained previously by calling `list()`.
            pub fn deleteVertex(self: *AllVertices, vtx: Vertex) bool {
                return self.vertices.swapRemove(vtx);
            }
        };

        ///////////////// CLASSIC TRAVERSING ALGORITHMS /////////////////

        pub const Path = std.ArrayList(Vertex);
        const Tq = std.TailQueue(Vertex);
        pub const TreeTraverseAlg = enum { bfs, pre, post, ino };

        /// Traverses the graph, assuming it was created as a left-right tree.
        /// Returns an ArrayList with a result that needs a separate deinit() call.
        /// Gives a compiler error if `.undirected` graph. Supports four algorithms,
        /// one of which should be passed as `algorithm` param: \
        /// `.bfs` - breadth-first. `.pre` - preorder. `.post`- postorder \
        /// `.ino` - inorder; reports: child1, parent, child2.
        ///
        /// @WARNING: `traverseTree()` has no limits.
        /// All algorithms except `.ino` will cause you an (1) infinite loop if your tree has cycles
        /// or (2) unfeasible task if your tree is *NOT BINARY AND HAS MANY EDGES
        /// CONNECTED AT RANDOM DOWN THE TREE*. A small 100-node general tree with less than
        /// 20 downstream edges for each node connecting randomly to other nodes
        /// is one such unfeasible example.
        pub fn traverseTree(self: *Self, origin: Vertex, algorithm: TreeTraverseAlg) !Path {
            if (graphType == .undirected)
                @compileError("traverseTree works only with .directed, acyclic graphs");

            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            var path = Path.init(self.alloc);
            var result = TreeTraverse{ .path = &path, .cnt = self, .target = null };
            switch (algorithm) {
                .bfs => try result.breadthFirst(origin),
                .pre => try result.preorder(origin),
                .post => try result.postorder(origin),
                .ino => try result.inorder(origin),
            }
            return path;
        }
        /// The same rules as for `traverseTree()` function, except you should give a target
        /// vertex. The function will stop traversing once the target has been found.
        /// The target will be the last element in the returned ArrayList.
        pub fn traverseTreeIfTarget(self: *Self, origin: Vertex, target: Vertex, algorithm: TreeTraverseAlg) !Path {
            if (graphType == .undirected)
                @compileError("traverseTree works only with .directed, acyclic graphs");

            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;
            if (!self.gotVertex(target))
                return GraphError.MissingVertex;

            var path = Path.init(self.alloc);
            var result = TreeTraverse{ .path = &path, .cnt = self, .target = target };
            switch (algorithm) {
                .bfs => try result.breadthFirst(origin),
                .pre => try result.preorderIfTarget(origin),
                .post => try result.postorderIfTarget(origin),
                .ino => try result.inorderIfTarget(origin),
            }
            return path;
        }

        const TreeTraverse = struct {
            path: *Path,
            cnt: *Self,
            target: ?Vertex,
            stop: bool = false,

            pub fn breadthFirst(self: *TreeTraverse, origin: Vertex) !void {
                var queue = Tq{};
                var cache = Cache(Tq.Node).init(self.cnt.alloc);
                defer cache.deinit();

                var fringe = try cache.new();
                fringe.*.data = origin;
                queue.append(fringe);

                if (self.target == null) {
                    while (queue.len > 0) {
                        var elder = queue.popFirst().?;
                        try self.path.append(elder.*.data); //*
                        if (self.cnt.adjacentVertices(elder.*.data, .outgoing)) |children| {
                            for (children) |vtx| {
                                fringe = try cache.new();
                                fringe.*.data = vtx;
                                queue.append(fringe);
                            }
                        }
                        cache.delete(elder);
                    }
                } else {
                    while (queue.len > 0) {
                        var elder = queue.popFirst().?;
                        try self.path.append(elder.*.data); //*
                        if (eql(elder.*.data, self.target.?))
                            break;
                        if (self.cnt.adjacentVertices(elder.*.data, .outgoing)) |children| {
                            for (children) |vtx| {
                                fringe = try cache.new();
                                fringe.*.data = vtx;
                                queue.append(fringe);
                            }
                        }
                        cache.delete(elder);
                    }
                }
            }
            fn preorder(self: *TreeTraverse, origin: Vertex) !void {
                try self.path.append(origin); //*
                if (self.cnt.adjacentVertices(origin, .outgoing)) |children| {
                    for (children) |vtx| {
                        try self.preorder(vtx);
                    }
                }
            }
            fn preorderIfTarget(self: *TreeTraverse, origin: Vertex) !void {
                if (!eql(origin, self.target.?) and !self.stop) {
                    try self.path.append(origin); //*
                } else if (eql(origin, self.target.?)) {
                    try self.path.append(origin); //*
                    self.stop = true;
                    return;
                } else return;

                if (self.cnt.adjacentVertices(origin, .outgoing)) |children| {
                    for (children) |vtx| {
                        try self.preorderIfTarget(vtx);
                    }
                }
            }
            fn postorder(self: *TreeTraverse, origin: Vertex) !void {
                if (self.cnt.adjacentVertices(origin, .outgoing)) |children| {
                    for (children) |vtx| {
                        try self.postorder(vtx);
                    }
                    try self.path.append(origin); //*
                }
            }
            fn postorderIfTarget(self: *TreeTraverse, origin: Vertex) !void {
                if (self.cnt.adjacentVertices(origin, .outgoing)) |children| {
                    for (children) |vtx| {
                        try self.postorderIfTarget(vtx);
                    }

                    if (!eql(origin, self.target.?) and !self.stop) {
                        try self.path.append(origin); //*
                    } else if (eql(origin, self.target.?)) {
                        try self.path.append(origin); //*
                        self.stop = true;
                        return;
                    } else return;
                }
            }
            fn inorder(self: *TreeTraverse, origin: Vertex) !void {
                if (self.cnt.adjacentVertices(origin, .outgoing)) |kids| {
                    if (kids.len > 0)
                        try self.inorder(kids[0]);
                    try self.path.append(origin); //*
                    if (kids.len > 1)
                        try self.inorder(kids[1]);
                }
            }
            fn inorderIfTarget(self: *TreeTraverse, origin: Vertex) !void {
                if (self.cnt.adjacentVertices(origin, .outgoing)) |kids| {
                    if (kids.len > 0)
                        try self.inorderIfTarget(kids[0]);

                    if (!eql(origin, self.target.?) and !self.stop) {
                        try self.path.append(origin); //*
                    } else if (eql(origin, self.target.?)) {
                        try self.path.append(origin); //*
                        self.stop = true;
                        return;
                    } else return;

                    if (kids.len > 1)
                        try self.inorderIfTarget(kids[1]);
                }
            }
        };

        // AUX TYPES and TOOLS //

        // Define infinity
        pub const inf = if (@typeInfo(edgeWt) == .Float)
            @as(f64, @bitCast(std.math.inf(f64)))
        else
            @as(u64, @bitCast(std.math.inf(f64)));
        const INF = @TypeOf(inf);

        pub const WeightAndEdge = struct { weight: INF, edge: Edge };
        pub const HashMap3 = std.ArrayHashMap(Vertex, WeightAndEdge, AutoContext(Vertex), false);

        fn compareWeights(isMin: bool, a: INF, b: INF) bool {
            _ = isMin;
            return a <= b;
        }
        pub const mH = pieQ(INF, Vertex, .min, compareWeights);

        pub const MultiMap = union(enum(u2)) { bfsDfs: HashMap2, dij: HashMap3 };

        pub const SearchAlg = enum { bfs, dfsA, dfsB, dfsC, dij };
        pub const TraverseControl = enum { all, except, through };

        pub const Target = union(enum(u2)) { yes: Vertex, no: void };
        pub const noTarget = Target{ .no = {} };

        // Connection functions //

        /// A result of the graph traversing process, struct containing a connection tree,
        /// and methods to work on the result.
        pub const Connections = struct {
            /// Origin vertex from which the connection tree was computed
            origin: Vertex,
            /// ArrayList holding the path from the origin to the given destination vertex
            /// and obtained by using `getPathTo()` method. Refreshed each time
            /// `getPathTo()` is called on the new destination vertex.
            path: Path,
            /// ArrayHAshMap containing the connection tree.
            /// The map is packed into a MultiMap union depending on the algorithm used: \
            ///  `MultiMap = union(enum(u2)) { bfsDfs: HashMap2, dij: HashMap3 };`
            discovered: MultiMap,
            /// Every time the discovered tree is recalculated on the new destination,
            /// we reset last_lookup variable.
            last_lookup: Vertex = undefined,

            /// Call this method to release the Connections tree's backing allocations
            /// and free memory.
            pub fn deinit(self: *Connections) void {
                self.origin = undefined;
                self.last_lookup = undefined;
                self.path.deinit();
                const map_tag = Tag(self.discovered);
                switch (map_tag) {
                    .bfsDfs => self.discovered.bfsDfs.deinit(),
                    .dij => self.discovered.dij.deinit(),
                }
            }
            /// Checks wether the give destination vertex is present in the Connections tree.
            pub fn connectedTo(self: Connections, dest: Vertex) bool {
                const map_tag = Tag(self.discovered);
                return switch (map_tag) {
                    .bfsDfs => self.discovered.bfsDfs.contains(dest),
                    .dij => if (self.getDistanceTo(dest) > 0) true else false,
                };
            }
            /// Returns all vertices part of the Connections tree.
            pub fn getAllConnected(self: *Connections) []Vertex {
                const map_tag = Tag(self.discovered);
                return switch (map_tag) {
                    .bfsDfs => self.discovered.bfsDfs.keys(),
                    .dij => self.discovered.dij.keys(),
                };
            }
            /// Returns the numeric distance from the origin that was used
            /// to calculate the Connections tree to the given destination vertex.
            ///
            /// @IMPORTANT: This method only works if the Connections tree
            /// was built using the Dijkstra algorithm. Otherwise, it returns `0` (zero),
            /// but not an error.
            pub fn getDistanceTo(self: *Connections, dest: Vertex) INF {
                return switch (Tag(self.discovered)) {
                    .dij => if (self.discovered.dij.get(dest)) |d| d.weight else 0,
                    else => 0,
                };
            }
            /// Returns a Slice of Vertices from the origin that was used to
            /// calculate the Connections tree to the destination.
            ///
            /// @IMPORTANT: the path is built on reverse as it was
            /// reconstructed by the algorithm.
            /// The first "tile" on the path is the last element of this slice.
            /// If you need to walk it, use `walkPathTo()` or `popPathTo()` methods.
            /// If you need only the first "tile", use `getPathTakeFirst()`.
            pub fn getPathTo(
                self: *Connections,
                dest: Vertex,
            ) GraphError![]Vertex {
                try self.checkPath(dest);
                return self.path.items;
            }
            /// Returns the first vertex on the path or null if the path is empty.
            pub fn getPathTakeFirstTo(
                self: *Connections,
                dest: Vertex,
            ) GraphError!?Vertex {
                try self.checkPath(dest);
                return self.path.getLastOrNull();
            }
            /// Walks the path from the origin that was used to calculate the
            /// Connections tree to the destination.
            /// Akin to a regular Iterator. Engage with `next()`.
            pub fn walkPathTo(self: *Connections, dest: Vertex) GraphError!WalkPath {
                try self.checkPath(dest);
                return WalkPath{ .cnt = self, .idx = self.path.items.len };
            }
            const WalkPath = struct {
                cnt: *Connections,
                // path: Path,
                idx: u64,
                /// Moves forward the path till null.
                pub fn next(self: *WalkPath) ?Vertex {
                    if (self.idx == 0) return null;
                    while (self.idx > 0) {
                        self.idx -= 1;
                        return self.cnt.path.items[self.idx];
                    }
                    return null;
                }
                pub fn reset(self: *WalkPath) void {
                    self.idx = self.cnt.path.*.items.len;
                }
            };
            /// Pops the path from the origin that was used to calculate the
            /// Connections tree to the destination. The difference with `walkPathTo`
            /// is that this method removes the vertices as it walks, sort of erasing
            /// the path behind.
            pub fn popPathTo(self: *Connections, dest: Vertex) GraphError!PopPath {
                try self.checkPath(dest);
                // But we will reset the last_lookup to undefined.
                // The popped path might be claimed in the future,
                // so it must be reconstructed.
                self.last_lookup = undefined;
                return PopPath{ .cnt = self, .dest = dest };
            }
            const PopPath = struct {
                cnt: *Connections,
                dest: Vertex,
                /// Moves forward the path till null.
                pub fn next(self: *PopPath) ?Vertex {
                    if (eql(self.cnt.last_lookup, self.dest)) self.cnt.last_lookup = undefined;
                    const vtx = self.cnt.path.popOrNull();
                    while (vtx) |exist| return exist;
                    return null;
                }
                /// Rebuilds the path again with previously supplied origin and dest;
                pub fn reset(self: *PopPath) !void {
                    try self.cnt.checkPath(self.dest);
                }
            };
            fn checkPath(self: *Connections, dest: Vertex) GraphError!void {
                if (!eql(self.last_lookup, dest)) {
                    self.path.clearRetainingCapacity();
                    try self.constructPath(dest);
                    self.last_lookup = dest;
                }
            }
            fn constructPath(
                self: *Connections,
                dest: Vertex,
            ) GraphError!void {
                var contains: bool = false;
                var wae: ?*WeightAndEdge = null;
                const map_tag = Tag(self.discovered);
                switch (map_tag) {
                    .bfsDfs => contains = self.discovered.bfsDfs.contains(dest),
                    .dij => wae = self.discovered.dij.getPtr(dest),
                }
                if (contains or wae != null) {
                    var edge: Edge = undefined;
                    try self.path.append(dest);
                    var walk = dest;
                    switch (map_tag) {
                        .dij => {
                            if (wae.?.*.weight != 0 and wae.?.*.weight != inf) {
                                // TODO check against inf loop ?
                                while (!eql(walk, self.origin)) {
                                    edge = self.discovered.dij.get(walk).?.edge;
                                    var parent = edge.opposite(walk);
                                    try self.path.append(parent);
                                    walk = parent;
                                }
                            } else wae.?.*.weight = 0;
                        },
                        // One call for bfs and dfsA/B/C algorithms
                        .bfsDfs => while (!eql(walk, self.origin)) {
                            edge = self.discovered.bfsDfs.get(walk).?;
                            var parent = edge.opposite(walk);
                            try self.path.append(parent);
                            walk = parent;
                        },
                    }
                }
                _ = self.path.popOrNull();
            }
        };

        pub const Connection = struct {
            /// *true* If the target has been found during the traversing process.
            /// *false* otherwise.
            found: bool,

            /// In explore, the user can obtains the path
            /// and methods to work with the Connections tree.
            explore: Connections,

            /// Call this method to release the Connections tree' backing allocations
            /// and free memory.
            pub fn deinit(self: *Connection) void {
                self.explore.deinit();
            }
        };
        /// The same rules as for `connectionTree()` method.
        /// If *depth* is given but not *target*, the traversing will explore
        /// the graph up to the given depth, setting the return struct's field target
        /// variable `found` to *false*. If the target is given, it will stop traversing
        /// as soon as it finds the target. The `found` variable will be set to *true*.
        ///
        /// `.explore` field contains methods to work on the result.
        /// In case if the target and depth are not needed at all,
        /// give them `null` and `0` (zero) respectively. However, consider
        /// `connectionTree()` method as a shortcut. The result needs a `deinit()` call.
        ///
        /// @IMPORTANT. depth parameter has different objectives for each algorithm.
        ///
        /// `.dij` and `.dfsC` -> how many vertices it needs to explore
        /// to find the target, if any. \
        /// `.bfs`, `.dfsA`, and `dfsB` -> how many layers of neighboring vertices
        /// it needs to explore to find the target, if any.
        ///
        /// The result *Connection* needs a `deinit()` call.
        pub fn connectionTreeIfTarget(self: *Self, origin: Vertex, target: ?Vertex, depth: u64, algorithm: SearchAlg) GraphError!Connection {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            const target_ = if (target) |t| Target{ .yes = t } else noTarget;
            var connections = switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, null, target_, depth, algorithm, .all),
                .dij => try self.dijkstraAll(origin, null, target_, depth, .all),
                .dfsC => try self.dfsCAll(origin, null, target_, depth, .all),
            };
            var found: bool = false;
            if (target) |t| {
                found = if ((try connections.getPathTo(t)).len > 0) true else false;
            }
            return .{ .found = found, .explore = connections };
        }
        /// The same rules as for `connectionTreeIfTarget()` method, except you shall
        /// give a knockout set of vertices that will be omitted
        /// during the traversing process, as they are not present in the graph at all.
        ///
        /// The result *Connections* needs a `deinit()` call.
        pub fn connectionTreeIfTargetExcept(self: *Self, origin: Vertex, knockout: *SetV, target: ?Vertex, depth: u64, algorithm: SearchAlg) GraphError!Connection {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            const target_ = if (target) |t| Target{ .yes = t } else noTarget;
            var connections = switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, knockout, target_, depth, algorithm, .except),
                .dij => try self.dijkstraAll(origin, knockout, target_, depth, .except),
                .dfsC => try self.dfsCAll(origin, knockout, target_, depth, .except),
            };
            var found: bool = false;
            if (target) |t| {
                found = if ((try connections.getPathTo(t)).len > 0) true else false;
            }
            return .{ .found = found, .explore = connections };
        }

        /// The same rules as for `connectionTreeIfTarget()` method, except you shall
        /// give a knockout set of vertices that will be used exclusively
        /// during the traversing process, as only they are present in the graph,
        /// and no other vertices.
        ///
        /// The result *Connections* needs a `deinit()` call.
        pub fn connectionTreeIfTargetThrough(self: *Self, origin: Vertex, knockout: *SetV, target: ?Vertex, depth: u64, algorithm: SearchAlg) GraphError!Connection {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            const target_ = if (target) |t| Target{ .yes = t } else noTarget;
            var connections = switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, knockout, target_, depth, algorithm, .through),
                .dij => try self.dijkstraAll(origin, knockout, target_, depth, .through),
                .dfsC => try self.dfsCAll(origin, knockout, target_, depth, .through),
            };
            var found: bool = false;
            if (target) |t| {
                found = if ((try connections.getPathTo(t)).len > 0) true else false;
            }
            return .{ .found = found, .explore = connections };
        }
        /// Calculates a Connections tree in between the given *origin* vertex and
        /// all other vertices in the graph via the supplied traversing algorithm.
        /// The returned struct *Connections* has methods to work on the result.
        /// Algorithms' enums are as follows:
        ///
        /// `.bfs` -> breadth-first search. \
        /// `.dfsA` -> simple iterative depth-first search.\
        /// `.dfsB` -> iterative depth-first search, emulates true recursion (.undirected graphs only!). \
        /// `.dfsC` -> recursive depth-first search. \
        /// `.dij` -> Dijkstra algorithm (.weighted graphs only!)
        ///
        /// The result *Connections* needs a `deinit()` call.
        pub fn connectionTree(self: *Self, origin: Vertex, algorithm: SearchAlg) GraphError!Connections {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            return switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, null, noTarget, 0, algorithm, .all),
                .dij => try self.dijkstraAll(origin, null, noTarget, 0, .all),
                .dfsC => try self.dfsCAll(origin, null, noTarget, 0, .all),
            };
        }
        /// The same rules as for `connectionTree()` method, except you shall
        /// give a knockout set of vertices that will be omitted
        /// during the traversing process, as they are not present in the graph at all.
        ///
        /// The result *Connections* needs a `deinit()` call.
        pub fn connectionTreeExcept(self: *Self, origin: Vertex, knockout: *SetV, algorithm: SearchAlg) GraphError!Connections {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            return switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, knockout, noTarget, 0, algorithm, .except),
                .dij => try self.dijkstraAll(origin, knockout, noTarget, 0, .except),
                .dfsC => try self.dfsCAll(origin, knockout, noTarget, 0, .except),
            };
        }
        /// The same rules as for `connectionTree()` method, except you shall
        /// give a knockout set of vertices that will be used exclusively
        /// during the traversing process, as only they are present in the graph,
        /// and no other vertices.
        ///
        /// The result *Connections* needs a `deinit()` call.
        pub fn connectionTreeThrough(self: *Self, origin: Vertex, knockout: *SetV, algorithm: SearchAlg) GraphError!Connections {
            if (!self.gotVertex(origin))
                return GraphError.MissingVertex;

            return switch (algorithm) {
                .bfs, .dfsA, .dfsB => try self.bfsDfsAll(origin, knockout, noTarget, 0, algorithm, .through),
                .dij => try self.dijkstraAll(origin, knockout, noTarget, 0, .through),
                .dfsC => try self.dfsCAll(origin, knockout, noTarget, 0, .through),
            };
        }

        // NON-PUBLIC API //

        // DIJKSTRA SHORTEST PATH //
        fn dijkstraAll(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            target: Target,
            depth: u64,
            control: TraverseControl,
        ) GraphError!Connections {
            switch (control) {
                .all => return self.dijkstra(origin, null, false, target, depth),
                .except => return self.dijkstra(origin, knockout.?, false, target, depth),
                .through => return self.dijkstra(origin, knockout.?, true, target, depth),
            }
        }
        fn dijkstraLoop(
            self: *Self,
            fringe: *Vertex,
            knockout: ?*SetV,
            discovered: *HashMap3,
            reflect: bool,
            minHeap: *mH,
        ) GraphError!void {
            var items = self.outGoing.get(fringe.*).?.iterator();
            while (items.next()) |Item| {
                var v = Item.key_ptr.*;
                if (knockout != null and !reflect and knockout.?.contains(v)) continue;
                if (knockout != null and reflect and !knockout.?.contains(v)) continue;

                var edge = Item.value_ptr.*;
                var dist = discovered.get(fringe.*).?.weight + edge.weight;

                var wae: *WeightAndEdge = discovered.getPtr(v).?;
                if (wae.*.weight > dist) {
                    wae.*.weight = dist;

                    wae.*.edge = edge; // Collect the edge
                    try minHeap.push(.{ .key = dist, .val = v });
                }
            }
        }
        fn dijkstra(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            reflect: bool,
            target: Target,
            depth: u64,
        ) GraphError!Connections {
            if (graphMode == .unweighted)
                return GraphError.DijkstraAlg_With_Unweighted_Graph;

            // Compute discovered map
            var discovered = HashMap3.init(self.alloc);
            var path = Path.init(self.alloc);

            var minHeap = mH.init(self.alloc);
            defer minHeap.deinit();

            for (self.outGoing.keys()) |vtx| {
                if (knockout != null and !reflect and knockout.?.contains(vtx)) continue;
                if (knockout != null and reflect and !knockout.?.contains(vtx)) continue;

                try discovered.put(vtx, .{ .weight = inf, .edge = undefined });
            }
            try discovered.put(origin, .{ .weight = 0, .edge = undefined });
            try minHeap.push(.{ .key = 0, .val = origin });

            switch (Tag(target)) {
                .yes => {
                    switch (depth) {
                        0 => while (!minHeap.isEmpty()) {
                            var fringe = (try minHeap.pop()).val;

                            try self.dijkstraLoop(&fringe, knockout, &discovered, reflect, &minHeap);
                            if (eql(target.yes, fringe)) break;
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (!minHeap.isEmpty() and depth_ > 0) : (depth_ -= 1) {
                                var fringe = (try minHeap.pop()).val;
                                try self.dijkstraLoop(&fringe, knockout, &discovered, reflect, &minHeap);
                                if (eql(target.yes, fringe)) break;
                            }
                        },
                    }
                },
                .no => {
                    switch (depth) {
                        0 => while (!minHeap.isEmpty()) {
                            var fringe = (try minHeap.pop()).val;
                            try self.dijkstraLoop(&fringe, knockout, &discovered, reflect, &minHeap);
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (!minHeap.isEmpty() and depth_ > 0) : (depth_ -= 1) {
                                var fringe = (try minHeap.pop()).val;
                                try self.dijkstraLoop(&fringe, knockout, &discovered, reflect, &minHeap);
                            }
                        },
                    }
                },
            }
            _ = discovered.swapRemove(origin);
            return .{ .origin = origin, .path = path, .discovered = .{ .dij = discovered } };
        }

        // BREADTH-FIRST SEARCH adn DEPTH-FIRST A, B combined //
        fn bfsDfsAll(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            target: Target,
            depth: u64,
            algorithm: SearchAlg,
            control: TraverseControl,
        ) GraphError!Connections {
            if (algorithm == .dfsB and graphType == .directed)
                return GraphError.dfsB_With_Directed_Graph;

            switch (control) {
                .all => return switch (algorithm) {
                    .bfs, .dfsA => try self.bfsDfsA(origin, null, false, target, depth, algorithm),
                    .dfsB => try self.dfsB(origin, null, false, target, depth),
                    else => GraphError.UnknownAlgorithm,
                },
                .except => return switch (algorithm) {
                    .bfs, .dfsA => try self.bfsDfsA(origin, knockout.?, false, target, depth, algorithm),
                    .dfsB => try self.dfsB(origin, knockout, false, target, depth),
                    else => GraphError.UnknownAlgorithm,
                },
                .through => return switch (algorithm) {
                    .bfs, .dfsA => try self.bfsDfsA(origin, knockout.?, true, target, depth, algorithm),
                    .dfsB => try self.dfsB(origin, knockout, true, target, depth),
                    else => GraphError.UnknownAlgorithm,
                },
            }
        }

        // BREADTH-FIRST and DEPTH-FIRST SEARCH A combined //
        fn bfsDfsLoop(self: *Self, origin: *Vertex, discovered: *HashMap2, knockout: ?*SetV, reflect: bool, cache: anytype, found: anytype) !void {
            var edges = self.incidentEdges(origin.*, .outgoing).?;
            for (edges) |edge| {
                var v = edge.opposite(origin.*);
                if (knockout != null and !reflect and knockout.?.contains(v))
                    continue;
                if (knockout != null and reflect and !knockout.?.contains(v))
                    continue;

                if (!discovered.contains(v)) {
                    try discovered.put(v, edge);

                    var fringe = try cache.new();
                    fringe.*.data = v;
                    found.append(fringe);
                }
            }
        }
        fn bfsDfsA(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            reflect: bool,
            target: Target,
            depth: u64,
            algorithm: SearchAlg,
        ) GraphError!Connections {
            var path = Path.init(self.alloc);
            var discovered = HashMap2.init(self.alloc);

            var found = Tq{};

            var cache = Cache(Tq.Node).init(self.alloc);
            defer cache.deinit();

            var fringe = try cache.new();
            fringe.*.data = origin;
            found.append(fringe);

            switch (Tag(target)) {
                .yes => {
                    switch (depth) {
                        0 => while (found.len > 0) {
                            var origin_ = switch (algorithm) {
                                .dfsA => found.pop().?,
                                else => found.popFirst().?,
                            };
                            if (eql(target.yes, origin_.*.data)) break;

                            try self.bfsDfsLoop(&origin_.*.data, &discovered, knockout, reflect, &cache, &found);

                            cache.delete(origin_);
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (found.len > 0 and depth_ > 0) : (depth_ -= 1) {
                                var origin_ = switch (algorithm) {
                                    .dfsA => found.pop().?,
                                    else => found.popFirst().?,
                                };
                                if (eql(target.yes, origin_.*.data)) break;

                                try self.bfsDfsLoop(&origin_.*.data, &discovered, knockout, reflect, &cache, &found);

                                cache.delete(origin_);
                            }
                        },
                    }
                },
                .no => {
                    switch (depth) {
                        0 => while (found.len > 0) {
                            var origin_ = switch (algorithm) {
                                .dfsA => found.pop().?,
                                else => found.popFirst().?,
                            };
                            try self.bfsDfsLoop(&origin_.*.data, &discovered, knockout, reflect, &cache, &found);

                            cache.delete(origin_);
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (found.len > 0 and depth_ > 0) : (depth_ -= 1) {
                                var origin_ = switch (algorithm) {
                                    .dfsA => found.pop().?,
                                    else => found.popFirst().?,
                                };
                                try self.bfsDfsLoop(&origin_.*.data, &discovered, knockout, reflect, &cache, &found);

                                cache.delete(origin_);
                            }
                        },
                    }
                },
            }
            _ = discovered.swapRemove(origin);
            return .{ .origin = origin, .path = path, .discovered = .{ .bfsDfs = discovered } };
        }

        // DEPTH-FIRST SEARCH B, true recursion emulation //
        fn dfsBLoop(self: *Self, origin: *Vertex, discovered: *HashMap2, knockout: ?*SetV, reflect: bool, found: *Path, leftover: *Path) !void {
            for (self.incidentEdges(origin.*, .outgoing).?) |edge| {
                var v = edge.opposite(origin.*);
                if (knockout != null and !reflect and knockout.?.contains(v))
                    continue;
                if (knockout != null and reflect and !knockout.?.contains(v))
                    continue;

                if (!discovered.contains(v)) {
                    try discovered.put(v, edge);
                    try found.append(v);

                    for (self.incidentEdges(v, .outgoing).?) |edge_| {
                        var v_ = edge_.opposite(v);
                        if (discovered.contains(v_)) {
                            try leftover.append(v_);
                        }
                    }
                    break;
                }
            }
        }
        fn dfsB(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            reflect: bool,
            target: Target,
            depth: u64,
        ) GraphError!Connections {
            var path = Path.init(self.alloc);
            var discovered = HashMap2.init(self.alloc);

            var found = Path.init(self.alloc);
            defer found.deinit();
            var leftover = Path.init(self.alloc);
            defer leftover.deinit();

            try found.append(origin);

            const vtx_count = self.vertexCount();

            switch (Tag(target)) {
                .yes => {
                    switch (depth) {
                        0 => while (found.items.len > 0 and discovered.count() < vtx_count) {
                            var origin_ = found.pop();

                            if (eql(target.yes, origin_)) break;

                            try self.dfsBLoop(&origin_, &discovered, knockout, reflect, &found, &leftover);

                            if (found.items.len == 0) {
                                try found.resize(leftover.items.len);
                                std.mem.swap([]Vertex, &found.items, &leftover.items);
                                leftover.clearRetainingCapacity();
                            }
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (found.items.len > 0 and depth_ > 0) : (depth_ -= 1) {
                                var origin_ = found.pop();

                                if (eql(target.yes, origin_)) break;

                                try self.dfsBLoop(&origin_, &discovered, knockout, reflect, &found, &leftover);

                                if (found.items.len == 0) {
                                    try found.resize(leftover.items.len);
                                    std.mem.swap([]Vertex, &found.items, &leftover.items);
                                    leftover.clearRetainingCapacity();
                                }
                            }
                        },
                    }
                },
                .no => {
                    switch (depth) {
                        0 => while (found.items.len > 0 and discovered.count() < vtx_count) {
                            var origin_ = found.pop();

                            try self.dfsBLoop(&origin_, &discovered, knockout, reflect, &found, &leftover);

                            if (found.items.len == 0) {
                                try found.resize(leftover.items.len);
                                std.mem.swap([]Vertex, &found.items, &leftover.items);
                                leftover.clearRetainingCapacity();
                            }
                        },
                        else => {
                            var depth_: u64 = depth;
                            while (found.items.len > 0 and depth_ > 0) : (depth_ -= 1) {
                                var origin_ = found.pop();

                                try self.dfsBLoop(&origin_, &discovered, knockout, reflect, &found, &leftover);

                                if (found.items.len == 0) {
                                    try found.resize(leftover.items.len);
                                    std.mem.swap([]Vertex, &found.items, &leftover.items);
                                    leftover.clearRetainingCapacity();
                                }
                            }
                        },
                    }
                },
            }
            _ = discovered.swapRemove(origin);
            return .{ .origin = origin, .path = path, .discovered = .{ .bfsDfs = discovered } };
        }

        // DEPTH-FIRST SEARCH C, recursive//
        fn dfsCAll(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            target: Target,
            depth: u64,
            control: TraverseControl,
        ) GraphError!Connections {
            var path = Path.init(self.alloc);
            var discovered = HashMap2.init(self.alloc);
            switch (control) {
                .all => try self.dfsC(origin, null, &discovered, false, target, depth),
                .except => try self.dfsC(origin, knockout.?, &discovered, false, target, depth),
                .through => try self.dfsC(origin, knockout.?, &discovered, true, target, depth),
            }
            _ = discovered.swapRemove(origin);
            return .{ .origin = origin, .path = path, .discovered = .{ .bfsDfs = discovered } };
        }
        // dfsC is a switch functions on 4 others
        fn dfsC(
            self: *Self,
            origin: Vertex,
            knockout: ?*SetV,
            discovered: *HashMap2,
            reflect: bool,
            target: Target,
            depth: u64,
        ) GraphError!void {
            var depth_: u64 = depth;
            var run = dfcCALL_{
                .self = self,
                .knockout = knockout,
                .discovered = discovered,
                .reflect = reflect,
                .target = target,
                .depth = &depth_,
            };
            switch (Tag(target)) {
                .yes => {
                    switch (depth) {
                        0 => return try run.recur(origin),
                        else => return try run.recur1(origin),
                    }
                },
                .no => {
                    switch (depth) {
                        0 => return try run.recur2(origin),
                        else => return try run.recur3(origin),
                    }
                },
            }
        }
        const dfcCALL_ = struct {
            self: *Self,
            knockout: ?*SetV,
            discovered: *HashMap2,
            reflect: bool,
            target: Target,
            depth: *u64,

            const Self_ = @This();

            /// Target, no depth
            fn recur(self_: *Self_, origin: Vertex) GraphError!void {
                const edges = self_.self.incidentEdges(origin, .outgoing).?;
                for (edges) |edge| {
                    var v = edge.opposite(origin);
                    if (self_.knockout != null and !self_.reflect and self_.knockout.?.contains(v))
                        continue;
                    if (self_.knockout != null and self_.reflect and !self_.knockout.?.contains(v))
                        continue;

                    if (!self_.discovered.contains(v)) {
                        try self_.discovered.put(v, edge);
                        if (eql(self_.target.yes, v)) return;

                        try self_.recur(v);
                    }
                }
            }
            /// Target, depth
            fn recur1(self_: *Self_, origin: Vertex) GraphError!void {
                if (self_.depth.* == 0) return;

                const edges = self_.self.incidentEdges(origin, .outgoing).?;
                for (edges) |edge| {
                    var v = edge.opposite(origin);
                    if (self_.knockout != null and !self_.reflect and self_.knockout.?.contains(v))
                        continue;
                    if (self_.knockout != null and self_.reflect and !self_.knockout.?.contains(v))
                        continue;

                    if (!self_.discovered.contains(v)) {
                        try self_.discovered.put(v, edge);
                        if (eql(self_.target.yes, v)) return;

                        self_.depth.* -|= 1;
                        try self_.recur1(v);
                    }
                }
            }
            /// No target, no depth
            fn recur2(self_: *Self_, origin: Vertex) GraphError!void {
                const edges = self_.self.incidentEdges(origin, .outgoing).?;
                for (edges) |edge| {
                    var v = edge.opposite(origin);
                    if (self_.knockout != null and !self_.reflect and self_.knockout.?.contains(v))
                        continue;
                    if (self_.knockout != null and self_.reflect and !self_.knockout.?.contains(v))
                        continue;

                    if (!self_.discovered.contains(v)) {
                        try self_.discovered.put(v, edge);
                        try self_.recur2(v);
                    }
                }
            }
            /// No target, depth
            fn recur3(self_: *Self_, origin: Vertex) GraphError!void {
                if (self_.depth.* == 0) return;

                const edges = self_.self.incidentEdges(origin, .outgoing).?;
                for (edges) |edge| {
                    var v = edge.opposite(origin);
                    if (self_.knockout != null and !self_.reflect and self_.knockout.?.contains(v))
                        continue;
                    if (self_.knockout != null and self_.reflect and !self_.knockout.?.contains(v))
                        continue;

                    if (!self_.discovered.contains(v)) {
                        try self_.discovered.put(v, edge);

                        self_.depth.* -|= 1;
                        try self_.recur3(v);
                    }
                }
            }
        };

        // TOPOLOGICAL SORT //
        const Degree = std.ArrayHashMap(Vertex, usize, AutoContext(Vertex), false);

        /// Sorts vertices topologically into an ArrayList container
        /// and returns a struct with methods.
        /// You might get extra benefits from the fact that Musubi remembers the insertion
        /// order, so the return list will largely resemble the graph
        /// the way it was created or processed.
        /// The built-in `walk()` function walks the list.
        /// You can also query the location by using `getAny()`, giving it an order number.
        /// Also, there are `getFirst()` and `getLast()` functions.
        ///
        /// @IMPORTANT: Works only on directed graphs with no cycles (acyclic).
        /// If the graph is not acyclic, the returned struct' variable `acyclic`
        /// will be set to *false*; the returned list might miss vertices or be empty.
        ///
        /// Since it *creates* a list, it needs a `deinit()` call!
        pub fn topologicalSort(self: *Self) GraphError!Topo {
            if (graphType == .undirected)
                @compileError("Topological Sort works only on directed graphs!");

            var cache = Cache(Tq.Node).init(self.alloc);
            defer cache.deinit();

            // List containing the sorting result
            var topo = Path.init(self.alloc);
            // Linked list, collector for vertices with no incoming edges
            var freeVtx = Tq{};
            // Map, collector of degree numbers for each vertex
            var degrees = Degree.init(self.alloc);
            defer degrees.deinit();

            for (self.vertices()) |vtx| {
                // Get the number of incoming edges from every vertex in the graph
                var degree_: usize = self.degree(vtx, .incoming).?;
                try degrees.put(vtx, degree_);

                var fringe = try cache.new();
                fringe.*.data = vtx;
                // Take vertices with no incoming edges
                if (degree_ == 0) freeVtx.append(fringe);
            }

            while (freeVtx.len > 0) {
                var vtx = freeVtx.popFirst().?;
                // Add vertex to the result
                try topo.append(vtx.*.data);

                // Take all outgoing neighbors of just appended vertex
                for (self.incidentEdges(vtx.*.data, .outgoing).?) |edge| {
                    var vtx_ = edge.opposite(vtx.*.data);
                    var degree_ = degrees.getPtr(vtx_).?;
                    degree_.* -= 1; // reduce their degree by one

                    var fringe = try cache.new();
                    fringe.*.data = vtx_;
                    // If it happens to be a vertex with zero degree, put it into queue
                    if (degree_.* == 0) freeVtx.append(fringe);
                }
                cache.delete(vtx);
            }
            // TODO perhaps we need this condition too: (topo.items.len == 0 or ...
            const acyclic: bool = if (topo.items.len != self.vertexCount()) false else true;

            return Topo{ .topo = topo, .acyclic = acyclic };
        }
        pub const Topo = struct {
            /// ArrayList container which holds the `topologicalSort()` result.
            topo: Path,
            /// Variable indicating whether the graph is acyclic or not.
            acyclic: bool,
            /// Call to release the allocated memory used by topologicalSort().
            pub fn deinit(self: *Topo) void {
                self.topo.deinit();
            }
            /// Get the result in the form of a sorted array of vertices.
            pub fn getAll(self: *Topo) ?[]Vertex {
                if (self.topo.items.len == 0) return null;
                return self.topo.items;
            }
            /// Query a vertex of the sorted result by its position number.
            pub fn getPosition(self: *Topo, task_number: usize) ?Vertex {
                if (task_number >= self.topo.items.len or self.topo.items.len == 0) return null;
                return self.topo.items[task_number];
            }
            /// Get the first vertex of the result.
            pub fn getFirst(self: *Topo) ?Vertex {
                if (self.topo.items.len == 0) return null;
                return self.topo.items[0];
            }
            /// Get the last vertex of the result.
            pub fn getLast(self: *Topo) ?Vertex {
                if (self.topo.items.len == 0) return null;
                return self.topo.getLastOrNull();
            }
            /// Walk the result.
            pub fn walk(self: *Topo) WalkTopo {
                return .{ .cnt = self };
            }
            const WalkTopo = struct {
                cnt: *Topo,
                idx: u64 = 0,
                /// Moves forward the result.
                pub fn next(self: *WalkTopo) ?Vertex {
                    if (self.cnt.topo.items.len == 0) return null;

                    while (self.idx < self.cnt.topo.items.len) {
                        const answer = self.cnt.topo.items[self.idx];
                        self.idx += 1;
                        return answer;
                    }
                    return null;
                }
                /// Resets the walker to the first position of the result
                pub fn reset(self: *WalkTopo) void {
                    self.idx = 0;
                }
            };
        };

        // MINIMUM SPAN-TREE • Prim-Jarnik algorithm //

        const HashMap4 = std.ArrayHashMap(PairV, Edge, AutoContext(PairV), false);
        const WeightAndEdge2 = struct { weight: INF, edge: ?Edge };
        const HashMap5 = std.ArrayHashMap(Vertex, WeightAndEdge2, AutoContext(Vertex), false);

        /// Computes a Minimum Spanning Tree of the Graph using the Prim-Jarnik algorithm.
        /// The returned struct *MST* has methods to get access to vertices pairs
        /// and their corresponding edges.
        ///
        /// Works only on .undirected and .weighted graphs. *MST* requires a `deinit()` call.
        pub fn primJarnikMST(self: *Self) GraphError!MST {
            if (graphType != .undirected)
                @compileError("Prim-Jarnik MST algorithm works only on .undirected .weighted graphs\nCurrent graph is ." ++ @tagName(graphType));
            if (graphMode != .weighted)
                @compileError("Prim-Jarnik MST algorithm works only on .undirected .weighted graphs\nCurrent graph is ." ++ @tagName(graphMode));

            // Largely resembles Dijkstra algorithm

            // Init the tree and its cost
            var tree = HashMap4.init(self.alloc);
            // TODO to solve: the cost of the MST Tree might exceed the INF type's capacity!
            var cost: INF = 0;

            // Compute discovered map
            var discovered = HashMap5.init(self.alloc);
            defer discovered.deinit();

            var minHeap = mH.init(self.alloc);
            defer minHeap.deinit();

            var all_vertices = self.outGoing.keys();
            try discovered.put(all_vertices[0], .{ .weight = 0, .edge = null });

            for (all_vertices[1..]) |vtx| {
                try discovered.put(vtx, .{ .weight = inf, .edge = null });
            }
            try minHeap.push(.{ .key = 0, .val = all_vertices[0] });

            while (!minHeap.isEmpty()) {
                var fringe = (try minHeap.pop()).val;
                // Since we remove yet another vertex from the heap,
                // we remove it from the discovered map as well
                // to ensure that every subsequent call to discovered
                // runs only on relevant data. That gives almost x2 speed!
                var wae_ = discovered.fetchSwapRemove(fringe);
                if (wae_) |wae__| {
                    var wae: WeightAndEdge2 = wae__.value;
                    if (wae.edge) |edge| {
                        const endpoints = edge.endpoints();
                        var gop = try tree.getOrPut(endpoints);
                        if (!gop.found_existing) {
                            gop.value_ptr.* = edge;
                            cost += wae.weight;
                        }
                    }
                }

                var destinations = self.outGoing.get(fringe).?.iterator();
                while (destinations.next()) |Item| {
                    var v = Item.key_ptr.*;

                    var edge: *Edge = Item.value_ptr;
                    var wae2_ = discovered.getPtr(v);

                    if (wae2_) |wae2| {
                        if (edge.*.weight < wae2.*.weight) {
                            wae2.* = WeightAndEdge2{ .weight = edge.*.weight, .edge = edge.* };

                            try minHeap.push(.{ .key = edge.*.weight, .val = v });
                        }
                    }
                }
            }
            return MST{ .cost = cost, .tree = tree, .len = tree.count() };
        }
        const MST = struct {
            /// Total cost of the Minimum Spanning Tree, or the sum of all of its edges.
            cost: INF = 0,
            /// Minimum Spanning Tree in the form of a map.
            ///
            /// keys: vertex pairs.
            /// values: edges
            tree: HashMap4,
            /// Length of the Minimum Spanning Tree.
            len: u64,

            /// Call deinit to release the MST backing allocation and free memory.
            pub fn deinit(self: *MST) void {
                self.tree.deinit();
            }
            /// Returns the edges of the Minimum Spanning Tree.
            pub fn getEdges(self: *MST) []Edge {
                return self.tree.values();
            }
            /// Returns vertex pairs corresponded to edges of the Minimum Spanning Tree.
            pub fn getVertexPairs(self: *MST) []PairV {
                return self.tree.keys();
            }
            /// Check wether the Minimum Spanning Three contains the give vertex pair.
            pub fn gotVertexPair(self: *MST, endpoints: PairV) bool {
                return self.tree.contains(endpoints);
            }
        };

        // MINIMUM SPAN-TREE • Kruskal algorithm //

        fn compareWeights2(isMin: bool, a: edgeWt, b: edgeWt) bool {
            _ = isMin;
            return a <= b;
        }
        const mH2 = pieQ(edgeWt, Edge, .min, compareWeights2);

        /// Computes a Minimum Spanning Tree of the Graph using the Kruskal algorithm.
        /// The returned struct *MST* has method to work with result
        /// and requires a `deinit()` call.
        /// A part of the algorithm is a recursive function. Might cause overflows
        /// on graphs with more than 20M edges, but it is rather system/machine dependable.
        ///
        /// The algorithm needs x2 more aux memory than Prim-Jarnik and runs
        /// non-linearly slower as the graph's size increases.
        /// Another peculiarity over Prim-Jarnik is that
        /// it reports all edges in the non-decreasing order of their weights.
        ///
        /// Works only on .undirected and .weighted graphs.
        pub fn kruskalMST(self: *Self) GraphError!MST {
            if (graphType != .undirected)
                @compileError("Kruskal MST algorithm works only on .undirected .weighted graphs\nCurrent graph is ." ++ @tagName(graphType));
            if (graphMode != .weighted)
                @compileError("Kruskal MST algorithm works only on .undirected .weighted graphs\nCurrent graph is ." ++ @tagName(graphMode));

            // Type that manages clusters of vertices and finds a min-edge in between them.
            const Partition = struct {
                const Par = @This();

                const Position = struct {
                    size: u64,
                    parent: *Position,
                };
                fn find(self_: *Par, p: *Position) *Position {
                    if (!eql(p.parent, p)) {
                        p.parent = self_.find(p.parent);
                    }
                    return p.parent;
                }
                fn union_(self_: *Par, p: *Position, q: *Position) void {
                    _ = self_;

                    if (!eql(p, q)) {
                        if (p.size < q.size) {
                            q.parent = p;
                            p.size += q.size;
                        } else {
                            p.parent = q;
                            q.size += p.size;
                        }
                    }
                }
            };
            // Positions themselves are stored in the cache,
            var cache = Cache(Partition.Position).init(self.alloc);
            defer cache.deinit();

            var tree = HashMap4.init(self.alloc);

            var pq = mH2.init(self.alloc);
            defer pq.deinit();

            var forest = Partition{};

            // yet here, we only store pointers to positions located in the cache.
            var position = std.ArrayHashMap(Vertex, *Partition.Position, AutoContext(Vertex), false).init(self.alloc);
            defer position.deinit();

            var cost: INF = 0;

            for (self.vertices()) |vtx| {
                var pos = try cache.new();
                pos.*.size = 1;
                pos.*.parent = pos;

                try position.put(vtx, pos);

                for (self.incidentEdges(vtx, .outgoing).?) |edge| {
                    try pq.push(.{ .key = edge.weight, .val = edge });
                }
            }

            var size = self.vertexCount();
            while (tree.count() != size - 1 and !pq.isEmpty()) {
                const item = try pq.pop();
                var weight = item.key;
                var edge: Edge = item.val;

                const u = edge.origin;
                const v = edge.destination;

                var p = forest.find(position.get(u).?);
                var q = forest.find(position.get(v).?);
                if (!eql(p, q)) {
                    try tree.put(edge.endpoints(), edge);
                    cost += weight;
                    forest.union_(p, q);
                }
            }
            return MST{ .cost = cost, .tree = tree, .len = tree.count() };
        }
    };
}

const testing = std.testing;
const expect = testing.expect;
const print = std.debug.print;
var allocatorT = std.testing.allocator;

test "Musubi: basics" {
    const VertexId = []const u8;
    const EdgeId = []const u8;
    const EdgeWt = void;

    const Graph = Musubi(VertexId, EdgeId, EdgeWt, .undirected, .unweighted);
    var graph: Graph = .{};
    graph.init(allocatorT);
    defer graph.deinit();

    // Inserts and stats
    var lax = try graph.insertVertex("LAX");
    var sfo = try graph.insertVertex("SFO");
    try expect(eql(lax.id, "LAX"));
    try expect(eql(sfo.id, "SFO"));
    try expect(graph.vertexCount() == 2);
    try expect(graph.edgeCount() == 0);

    var a0000 = try graph.insertEdge(lax, sfo, "A0000", {});
    try expect(graph.degree(lax, .outgoing) == 1);
    try expect(graph.gotEdge(lax, sfo));
    try expect(eql(graph.getEdge(lax, sfo), a0000));
    try expect(graph.edgeCount() == 1);

    var tex = try graph.insertVertex("TEX");
    var por = try graph.insertVertex("POR");
    try expect(graph.vertexCount() == 4);

    var b1111 = try graph.insertEdge(tex, por, "B1111", {});
    try expect(graph.degree(tex, .outgoing) == 1);
    try expect(graph.gotEdge(tex, por));
    try expect(eql(graph.getEdge(tex, por), b1111));
    try expect(graph.edgeCount() == 2);

    var c2222 = try graph.insertEdge(lax, tex, "C2222", {});
    try expect(graph.degree(lax, .outgoing) == 2);
    try expect(graph.gotEdge(lax, tex));
    try expect(eql(graph.getEdge(lax, tex), c2222));
    try expect(eql(c2222.origin, lax));
    try expect(eql(c2222.destination, tex));
    try expect(eql(c2222.id, "C2222"));
    try expect(eql(c2222.opposite(lax), tex));

    try expect(eql(graph.getEdge(lax, por), null));
    try expect(graph.edgeCount() == 3);

    // You can own a vertex outside the graph
    var han = graph.makeVertex("HAN");
    var lwo = graph.makeVertex("LWO");

    // You can own an edge outside the graph as well
    var han_lvl = graph.makeEdge(han, lwo, "HAN-LWO", {});

    // Insert them back into the graph with immediate testing
    try graph.insertVertexIfVertex(han);
    try graph.insertVertexIfVertex(lwo);
    try graph.insertEdgeIfEdge(han_lvl);

    try expect(eql(graph.gotVertex(han), graph.gotVertex(lwo)));
    try expect(eql(graph.gotEdgeIfEdge(han_lvl), true));
    try expect(graph.vertexCount() == 6);
    try expect(graph.edgeCount() == 4);

    // Gather all edges into a Set
    var edges = try graph.edgesIntoSet();
    defer edges.deinit();

    try expect(edges.count() == graph.edgeCount());
    try expect(edges.gotEdge(a0000));
    try expect(edges.gotEdge(b1111));
    try expect(edges.gotEdge(c2222));
    try expect(edges.gotEdge(han_lvl));

    // Iterate over all edges of the graph.
    for (edges.list()) |edge| {
        try expect(graph.gotEdgeIfEdge(edge));
    }

    // Gather all vertices into a Set
    var vertices = try graph.verticesIntoSet();
    defer vertices.deinit();

    try expect(vertices.count() == graph.vertexCount());
    try expect(vertices.gotVertex(lax));
    try expect(vertices.gotVertex(sfo));
    try expect(vertices.gotVertex(tex));
    try expect(vertices.gotVertex(por));
    try expect(vertices.gotVertex(han));
    try expect(vertices.gotVertex(lwo));

    // Check adjacent vertices
    var adjacent_to_lax = graph.adjacentVertices(lax, .outgoing).?;
    try expect(adjacent_to_lax.len == graph.degree(lax, .outgoing).?);

    // Iterate over all vertices of the graph.
    for (vertices.list()) |vertex| {
        try expect(graph.gotVertex(vertex));
    }

    // Iterate over edges incident to a given vertex only.
    var edges_of_lax = graph.incidentEdges(lax, .outgoing).?;
    for (edges_of_lax) |edge| {
        try expect(edges.gotEdge(edge));
    }
    try expect(edges_of_lax.len == graph.degree(lax, .outgoing).?);

    // Remove the edge from the graph
    try expect(graph.removeEdge(c2222));
    try expect(graph.removeEdge(b1111));
    try expect(graph.removeEdge(a0000));
    try expect(eql(graph.getEdge(lax, tex), null));
    try expect(eql(graph.getEdge(tex, por), null));
    try expect(eql(graph.getEdge(lax, sfo), null));
    try expect(graph.edgeCount() == 1);

    // Remove the vertex from the graph
    try expect(graph.vertexCount() == 6);
    try expect(graph.removeVertex(lax));
    try expect(graph.vertexCount() == 5);
    try expect(graph.removeVertex(sfo));
    try expect(graph.removeVertex(tex));
    try expect(graph.removeVertex(por));
    try expect(graph.vertexCount() == 2);

    graph.clearRetainingCapacity();

    try expect(graph.vertexCount() == 0);
    try expect(graph.edgeCount() == 0);

    // Get new graph2
    var graph2: Graph = .{};
    graph2.init(allocatorT);
    defer graph2.deinit();

    // Insert vertices and edges
    var lax2 = try graph2.insertVertex("LAX2");
    var sfo2 = try graph2.insertVertex("SFO2");
    var tex2 = try graph2.insertVertex("TEX2");
    var por2 = try graph2.insertVertex("POR2");

    var a0000_2 = try graph2.insertEdge(lax2, sfo2, "A0000_2", {});
    var b1111_2 = try graph2.insertEdge(tex2, por2, "B1111_2", {});
    var c2222_2 = try graph2.insertEdge(lax2, tex2, "C2222_2", {});

    try expect(graph2.vertexCount() == 4);
    try expect(graph2.edgeCount() == 3);

    // Get new graph3
    var graph3: Graph = .{};
    graph3.init(allocatorT);
    defer graph3.deinit();

    // Insert vertices and edges
    var lax3 = try graph3.insertVertex("LAX3");
    var sfo3 = try graph3.insertVertex("SFO3");
    var tex3 = try graph3.insertVertex("TEX3");
    var por3 = try graph3.insertVertex("POR3");

    var a0000_3 = try graph3.insertEdge(lax3, sfo3, "A0000_3", {});
    var b1111_3 = try graph3.insertEdge(tex3, por3, "B1111_3", {});
    var c2222_3 = try graph3.insertEdge(lax3, tex3, "C2222_3", {});

    try expect(graph3.vertexCount() == 4);
    try expect(graph3.edgeCount() == 3);

    // Get new allocator instance
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer if (gpa.deinit() == .leak) {
        @panic("memory leak ...");
    };
    const allocatorG = gpa.allocator();

    // Clone graph2 with allocatorG
    var graph4: Graph = .{};
    try graph4.cloneIntoSelfWithAllocator(&graph2, allocatorG);
    defer graph4.deinit();

    try expect(graph4.vertexCount() == 4);
    try expect(graph4.edgeCount() == 3);

    try expect(eql(graph4.getEdge(lax2, sfo2), a0000_2));
    try expect(eql(graph4.getEdge(tex2, por2), b1111_2));
    try expect(eql(graph4.getEdge(lax2, tex2), c2222_2));

    // Merge graph3 into graph4
    try graph4.mergeIntoSelf(&graph3);

    try expect(graph4.vertexCount() == 8);
    try expect(graph4.edgeCount() == 6);

    try expect(eql(graph4.getEdge(lax3, sfo3), a0000_3));
    try expect(eql(graph4.getEdge(tex3, por3), b1111_3));
    try expect(eql(graph4.getEdge(lax3, tex3), c2222_3));

    // Add new edge to graph4
    var lax4 = try graph4.insertVertex("LAX4");
    var sfo4 = try graph4.insertVertex("SFO4");
    try expect(graph4.gotVertex(lax4));
    try expect(graph4.gotVertex(sfo4));

    var A0000_4 = try graph4.insertEdge(lax4, sfo4, "A0000_4", {});
    try expect(eql(graph4.getEdge(lax4, sfo4), A0000_4));

    try expect(graph4.vertexCount() == 10);
    try expect(graph4.edgeCount() == 7);

    // Clone graph2
    // The graph4 now becomes graph2, using graph2' allocator
    try graph4.cloneIntoSelf(&graph2);
    try expect(graph4.vertexCount() == 4);
    try expect(graph4.edgeCount() == 3);

    try expect(eql(graph4.getEdge(lax2, sfo2), a0000_2));
    try expect(eql(graph4.getEdge(tex2, por2), b1111_2));
    try expect(eql(graph4.getEdge(lax2, tex2), c2222_2));

    try graph4.insertVertexIfVertex(lax4);
    try graph4.insertVertexIfVertex(sfo4);
    try graph4.insertEdgeIfEdge(A0000_4);
    try expect(eql(graph4.getEdge(lax4, sfo4), A0000_4));

    try expect(graph4.vertexCount() == 6);
    try expect(graph4.edgeCount() == 4);

    // Clear graph4
    graph4.clearAndFree();

    try expect(graph4.vertexCount() == 0);
    try expect(graph4.edgeCount() == 0);
}
