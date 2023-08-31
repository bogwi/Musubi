// Adapted from Zig's Standard Library //

const std = @import("std");

pub fn getAutoHashStratFn(comptime K: type, comptime Context: type, comptime strategy: std.hash.Strategy) (fn (Context, K) u32) {
    return struct {
        fn hash(ctx: Context, key: K) u32 {
            _ = ctx;
            var hasher = std.hash.Wyhash.init(0);
            std.hash.autoHashStrat(&hasher, key, strategy);
            return @truncate(hasher.final());
        }
    }.hash;
}
pub fn getAutoEqlFn(comptime K: type, comptime Context: type) (fn (Context, K, K, usize) bool) {
    return struct {
        fn eql(ctx: Context, a: K, b: K, b_index: usize) bool {
            _ = b_index;
            _ = ctx;
            return std.meta.eql(a, b);
        }
    }.eql;
}
pub fn AutoContext(comptime K: type) type {
    return struct {
        pub const hash = getAutoHashStratFn(K, @This(), .Deep);
        pub const eql = getAutoEqlFn(K, @This());
    };
}
