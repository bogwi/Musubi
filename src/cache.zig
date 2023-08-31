const std = @import("std");

/// A wrapper around an allocator for keeping the memory
/// of previously deleted items aside but serving it when needed
/// instead of allocating new bytes.
/// https://zig.news/xq/cool-zig-patterns-gotta-alloc-fast-23h  (c)
pub fn Cache(comptime T: type) type {
    return struct {
        const List = std.TailQueue(T);

        arena: std.heap.ArenaAllocator,
        free: List = .{},

        pub fn init(allocator: std.mem.Allocator) @This() {
            return .{ .arena = std.heap.ArenaAllocator.init(allocator) };
        }
        pub fn deinit(self_: *@This()) void {
            self_.arena.deinit();
        }
        pub fn new(self_: *@This()) !*T {
            const obj = if (self_.free.popFirst()) |item|
                item
            else
                try self_.arena.allocator().create(List.Node);
            return &obj.data;
        }
        pub fn delete(self_: *@This(), obj: *T) void {
            const node = @fieldParentPtr(List.Node, "data", obj);
            self_.free.append(node);
        }
    };
}
