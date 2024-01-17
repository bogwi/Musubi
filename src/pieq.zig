// MIT (c) bogwi@rakumail.jp 2023 //

const std = @import("std");

/// Describes Queue's general orientation.
/// Pass one of the values as parameter when creating the Queue.
///
/// `.min`, makes the Queue min-oriented.
/// `.max`, makes the Queue max-oriented.
const Orientation = enum { min, max };

/// **Binary Heap based Priority Queue** for large data.
///
/// The Queue stores entries as structs of type ` Item{ key: Key, value: Value }`.
/// Key and Value can be of any type. Queue builds its structure by comparing keys.
/// The parameter `orientation` set to `.min`,
/// turns the Queue into *min-oriented*;
/// set `.max`, to *max-oriented*.
///
/// Since any type can be a Key, like a struct, array, or a function,
/// not always a real number, and if it is your case,
/// you need to craft a custom parameter `compareFn()`,
/// a function which performs a comparison on keys, to reflect the structure of your
/// keys accordingly, as it is shown in the testing section to this code.
/// Such mods are usually simple.
/// But it allows multiple Queues running on any keys you like with minimal extra coding.
pub fn PieQ(comptime Key: type, comptime Value: type, comptime orientation: Orientation, comptime compareFn: fn (isMin: bool, a: Key, b: Key) bool) type {
    return struct {
        const Self = @This();
        const Item = struct { key: Key, val: Value };

        const QueueError = error{
            OutOfSize,
            OutOfMemory,
            QueueIsEmpty,
            AttemptToModifyLockedRoot,
        };

        const Dashboard = struct {
            root_locked: bool = false,
            orientation: Orientation = orientation,
            allocator: std.mem.Allocator,
        };

        // non-public interface start...........
        inline fn parent(j: usize) usize {
            return if (j > 0) @divFloor((j - 1), 2) else 0;
        }
        inline fn left(j: usize) usize {
            return 2 * j + 1;
        }
        inline fn right(j: usize) usize {
            return 2 * j + 2;
        }
        inline fn has_left(self: *Self, j: usize) bool {
            return 2 * j + 1 < self.data.items.len;
        }
        inline fn has_right(self: *Self, j: usize) bool {
            return 2 * j + 2 < self.data.items.len;
        }
        inline fn swap(self: *Self, i: usize, j: usize) void {
            const temp = self.data.items[i];
            self.data.items[i] = self.data.items[j];
            self.data.items[j] = temp;
        }
        inline fn heapify(self: *Self) void {
            var j = parent(self.data.items.len - 1);
            while (j > 0) : (j -= 1) self.downheap(j);
            self.downheap(j);
        }
        fn upheap(self: *Self, j: usize) void {
            var par = parent(j);
            var j_ = j;
            while (j_ > 0 and compareFn(self.isMin(), self.data.items[j_].key, self.data.items[par].key)) {
                self.swap(j_, par);
                j_ = par;
                par = parent(j_);
            }
        }
        fn downheap(self: *Self, j: usize) void {
            var j_ = j;
            while (self.has_left(j_)) {
                var smallest = if (self.has_right(j_) and compareFn(self.isMin(), self.data.items[right(j_)].key, self.data.items[left(j_)].key))
                    right(j_)
                else
                    left(j_);

                if (compareFn(self.isMin(), self.data.items[smallest].key, self.data.items[j_].key)) {
                    self.swap(j_, smallest);
                    j_ = smallest;
                    continue;
                }
                break;
            }
        }
        // ..........non-public interface end.

        const Data = std.ArrayList(Item);
        data: Data = undefined,
        dashboard: Dashboard,

        /// Creates the Queue with the given Allocator.
        pub fn init(allocator: std.mem.Allocator) Self {
            return .{ .data = Data.init(allocator), .dashboard = .{ .allocator = allocator } };
        }
        /// Releases all Queue's allocated memory.
        pub fn deinit(self: *Self) void {
            self.data.deinit();
        }
        /// Returns `true` if the Queue has at least 1 `Item`.
        pub fn isEmpty(self: *Self) bool {
            return self.data.items.len == 0;
        }
        /// Returns `true` if the Queue is .min oriented.
        /// Otherwise returns `false`.
        pub fn isMin(self: *Self) bool {
            return self.dashboard.orientation == .min;
        }
        /// Returns the number of items the Queue contains.
        pub fn count(self: *Self) usize {
            return self.data.items.len;
        }
        /// Returns true if the root is locked.
        /// See `lockRoot()` method.
        pub fn isRootLocked(self: *Self) bool {
            return self.dashboard.root_locked == true;
        }
        /// Returns the enum literal *.min* if the Queue was initiated as min-oriented,
        /// or *.max*, if the Queue was initiated as max-oriented Queue.
        pub fn minOrMax(self: *Self) Orientation {
            return self.dashboard.orientation;
        }
        /// Pushes the given `Item` into the Queue, keeping the Queue invariant.
        pub fn push(self: *Self, item: Item) QueueError!void {
            try self.data.append(item);
            self.upheap(self.data.items.len - 1);
        }

        /// Removes and returns the root of the Queue.
        /// If Queue is .min oriented, the root Item has the *smallest* key, often called *best key*.
        ///
        /// *Important:* If the root is locked, it will remove
        /// and return the Item with the 2nd best key, then with the 3rd, and so on.
        /// Root will stay intact until you unlock it.
        /// See `lockRoot()` and `isRootLocked()` methods.
        pub fn pop(self: *Self) QueueError!Item {
            while (!self.isEmpty() and !self.dashboard.root_locked) {
                const item = self.data.swapRemove(0);
                self.downheap(0);
                return item;
            }
            while (!self.isEmpty() and self.dashboard.root_locked) {
                switch (self.count()) {
                    else => {
                        if (compareFn(self.isMin(), self.data.items[1].key, self.data.items[2].key)) {
                            var item = self.data.swapRemove(1);
                            self.downheap(1);
                            return item;
                        } else {
                            var item = self.data.swapRemove(2);
                            self.downheap(2);
                            return item;
                        }
                    },
                    2 => return self.data.pop(),

                    1 => return self.data.items[0],
                }
            }
            return QueueError.QueueIsEmpty;
        }
        /// Returns the current root of the Queue but not deletes it.
        pub fn getRoot(self: *Self) QueueError!Item {
            while (!self.isEmpty()) {
                return self.data.items[0];
            }
            return QueueError.QueueIsEmpty;
        }
        /// Returns the Item with the second-best key, one after the root.
        /// If only one item has left, the root itself, it will return the root.
        pub fn getAfterRoot(self: *Self) QueueError!Item {
            while (!self.isEmpty()) {
                switch (self.count()) {
                    else => {
                        return if (compareFn(self.isMin(), self.data.items[1].key, self.data.items[2].key))
                            self.data.items[1]
                        else
                            self.data.items[2];
                    },
                    2 => return self.data.items[1],

                    1 => return self.data.items[0],
                }
            }
            return QueueError.QueueIsEmpty;
        }

        /// Guarantees the root will never be pooped if called with a parameter *lock*
        /// as `true`. Likewise, called `false` will unlock the root back; this is the default.
        pub fn lockRoot(self: *Self, lock: bool) void {
            self.dashboard.root_locked = lock;
        }
        /// Changes the root of the Queue. The call might be considered O(1)
        /// if the new root **is** best for the Queue. If the new root's key **is not** the best for the Queue,
        /// it will modify the Queue's structure bringing the newly inserted Item down
        /// to its proper place to preserve Queue's invariance.
        ///
        /// *Important*: If the Queue has only one element, the root, `changeRoot()`
        /// is not akin to `push()` method; it is only for changing the root,
        /// and it will not increase Queue's size, only displace the root Item
        /// with the new one provided.
        ///
        /// However, if the root is locked, it will give you an error.
        /// A locked root can't be modified. See `lockRoot()`, `isRootLocked()` methods.
        pub fn changeRoot(self: *Self, item: Item) QueueError!void {
            if (self.isEmpty()) return QueueError.QueueIsEmpty;

            if (!self.dashboard.root_locked) {
                self.data.items[0] = item;
                self.downheap(0);
                return;
            } else return QueueError.AttemptToModifyLockedRoot;
        }

        /// Returns the exact copy of the current Queue using the same allocator.
        /// Separate call deinit() on the clone is required as well.
        pub fn clone(self: *Self) !Self {
            return .{ .data = try self.data.clone(), .dashboard = .{ .root_locked = self.dashboard.root_locked, .orientation = self.dashboard.orientation, .allocator = self.dashboard.allocator } };
        }
        /// Returns the same Queue only with the opposite orientation, using the same allocator.
        /// If the Queue was initiated as `.min`, it will return max Queue.
        /// Contrary, if the Queue was built as `.max`, it will return min Queue.
        /// The original stays intact. A separate call deinit() on the clone is required as well.
        pub fn cloneAsOpposite(self: *Self) QueueError!Self {
            var opposite = init(self.dashboard.allocator);

            opposite.dashboard.orientation = if (self.dashboard.orientation == Orientation.min) Orientation.max else Orientation.min;

            try opposite.meldIntoSelf(self.data.items);
            return opposite;
        }
        /// Wipes the Queue out clean, reducing its size to 0(zero).
        /// Think of it as an instant eraser.
        pub fn clear(self: *Self) void {
            self.data.clearAndFree();
        }

        /// Exports the Queue as a slice of key-value pairs that can be iterated.
        /// Common usage is for melding one Queue into another.
        /// See `meldIntoSelf()` and `meldIntoNew()` methods.
        pub fn exportAsIterable(self: *Self) []Item {
            return self.data.items;
        }
        /// Melds into itself from any indexable iterable type containing
        /// key-value pairs of the same type as initiated Queue, that is,
        /// struct `.{ key: Key, val: Value }`.
        /// Triggers a gentle compile error if:
        /// >(1) your iterable is *not* indexable,
        /// (2) key-value pairs are *not* of the format mentioned in the header and correct type.
        ///
        /// Common usage is melding one Queue into another.
        /// Write:
        /// > `firstQueue.meldIntoSelf(secondQueue.exportAsIterable())`.
        ///
        /// Regardless of the secondQueue orientation, melding into Self
        /// preserves firstQueue's orientation. See `Orientation` enum, `exportAsIterable()` method.
        /// If you need to meld an iterable array of key-value pairs or tuples,
        /// see `meldIntoSelfFromIndexable()` method;
        pub fn meldIntoSelf(self: *Self, iterable: anytype) QueueError!void {
            // compatibility test; yields a compile error if fails
            _ = iterable[0].key;
            _ = iterable[0].val;

            for (iterable) |item| try self.data.append(.{ .key = item.key, .val = item.val });
            self.heapify();
        }
        /// Melds into a new Queue from any indexable iterable type containing
        /// key-value pairs of the same type as initiated Queue, that is,
        /// struct `.{ key: Key, val: Value }`.
        /// Triggers a gentle compile error if:
        /// >(1) your iterable is *not* indexable,
        /// (2) key-value pairs are *not* of the format mentioned in the header and correct type.
        ///
        /// Common usage is melding two Queues together into a third.
        /// Write:
        /// > `var newQueue = firstQueue.meldIntoNew(secondQueue.exportAsIterable())`.
        ///
        /// Regardless of the secondQueue orientation, melding into a new
        /// gives the new firstQueue's orientation. See `Orientation` enum `exportAsIterable()` method.
        /// If you need to meld a Queue with an iterable array of key-value pairs
        /// into a new Queue, see `meldIntoNewFromIndexable()` method;
        pub fn meldIntoNew(self: *Self, iterable: anytype) QueueError!Self {
            var new = try self.clone();
            try new.meldIntoSelf(iterable);
            return new;
        }
        /// Melds into a new Queue an indexable iterable which contains
        /// key-value pairs without specified field names but accessible
        /// through indexing, like tuples `.{ Key, Value }`;
        /// `tuple[0]` means key, `tuple[1]` means value.
        /// Triggers a gentle compile error if:
        /// >(1) your iterable is *not* indexable,
        /// (2) key-value pairs are *not* of the format mentioned in the header and correct type.
        ///
        /// Common usage is melding a Queue with an iterable array of key-value pairs into a new Queue.
        /// Write:
        /// > `var newQueue = yourQueue.meldIntoNewFromIndexable(&yourArrayOfKeyValuePairs)`
        ///
        /// Regardless of your array's orientation, melding into a new Queue
        /// gives Self's orientation. See `Orientation` enum.
        /// If you need to meld a Queue with another Queue into a new Queue,
        /// see `meldIntoNewFromIndexable()` method;
        pub fn meldIntoNewFromIndexable(self: *Self, iterable: anytype) QueueError!Self {
            var new = try self.clone();
            try new.meldIntoSelfFromIndexable(iterable);
            return new;
        }
        /// Melds into itself an indexable iterable which contains
        /// key-value pairs without specified field names but accessible
        /// through indexing, like tuples `.{ Key, Value }`;
        /// `tuple[0]` means key, `tuple[1]` means value.
        /// Triggers a gentle compile error if:
        /// >(1) your iterable is *not* indexable,
        /// (2) key-value pairs are *not* of the format mentioned in the header and correct type.
        ///
        /// Common usage is melding an iterable array of key-value pairs into Self.
        /// Write:
        /// > `yourQueue.meldIntoSelfFromIndexable(&yourArrayOfKeyValuePairs)`
        ///
        /// Regardless of your array's orientation, melding into Self
        /// preserves Self's orientation. See `Orientation` enum.
        /// If you need to meld one Queue into another, see `meldIntoSelf()` method.
        pub fn meldIntoSelfFromIndexable(self: *Self, iterable: anytype) QueueError!void {
            // compatibility test; yields a compile error if fails
            _ = iterable[0][0];
            _ = iterable[0][1];

            for (iterable) |item| try self.data.append(.{ .key = item[0], .val = item[1] });
            if (self.data.items.len > 1) self.heapify();
        }
    };
}

// TESTING BLOCK IS DISABLED //

// const expect = std.testing.expect;

// /// This is an example of a generic comparison function on unsigned 8-bit integers,
// /// which you must pass as compareFn parameter when initializing PieQ ADT.
// /// You will want to change the a and b types to reflect the keys your Queue is using.
// /// Leave the parameter isMin as it is, to boolean.
// /// Other variants of this function are presented in the testing section below.
// fn compareU8(isMin: bool, a: u8, b: u8) bool {
//     // a and b must match Queue's Key type! These are the keys we compare.
//     // But still, mismatch triggers compile stage panic: zig compiler bug: GenericPoison.
//     // This is for the bug identification if you will get it.
//     while (isMin)
//         return a <= b;
//     return a >= b;
// }
// test "Queue basics: understand minQueue and maxQueue" {
//     // create minQueue
//     var minQueue = PieQ(u8, u8, .min, compareU8).init(std.testing.allocator);
//     try expect(minQueue.minOrMax() == .min);
//     defer minQueue.deinit();

//     // push without declaring an Item type
//     const numbers = [_]u8{ 0, 3, 1, 4, 2, 5, 0 };
//     for (numbers) |num| try minQueue.push(.{ .key = num, .val = num });

//     // declaring an Item type might more explicit
//     // yet, this way Item will fit only this particular minQueue
//     // from which it was extracted
//     const Item = PieQ(u8, u8, .min, compareU8).Item;

//     // build an array of pairs using the declared type Item
//     var array_of_pairs = blk: {
//         var buffer: [10]Item = undefined;
//         for (&buffer, 0..) |*ptr, i| {
//             ptr.* = .{ .key = @intCast(i), .val = @intCast(i) };
//         }
//         break :blk buffer; // { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }
//     };

//     // meld array_of_pairs into minQueue and test whether the Queue is .min oriented
//     try minQueue.meldIntoSelf(&array_of_pairs);
//     const order = [_]u8{ 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8, 9 };
//     for (order) |digit| {
//         var item = try minQueue.pop();
//         // check minQueue against <order> array whether it pops items sorted from left to right
//         try expect(item.key == digit);
//     }
//     try expect(minQueue.count() == 0);
//     // popping an empty Queue, gives you an error
//     _ = minQueue.pop() catch |err| try std.testing.expectEqual(err, PieQ(u8, u8, .min, compareU8).QueueError.QueueIsEmpty);

//     // create maxQueue
//     var maxQueue = PieQ(u8, u8, .max, compareU8).init(std.testing.allocator);
//     try expect(maxQueue.minOrMax() == .max);
//     defer maxQueue.deinit();

//     // create an array from anonymous structs
//     var array_of_pairs2 = blk: {
//         var buffer: [17]struct { u8, u8 } = undefined;
//         for (&buffer, order) |*ptr, digit| {
//             ptr.* = .{ @intCast(digit), @intCast(digit) };
//         }
//         break :blk buffer;
//     };

//     // meld array_of_pairs2 into maxQueue and test whether the Queue is .max oriented
//     try maxQueue.meldIntoSelfFromIndexable(&array_of_pairs2);
//     var idx: usize = 1;
//     while (idx <= 17) : (idx += 1) {
//         var item = try maxQueue.pop();
//         // check maxQueue against <order> array but on reverse, from right to left
//         try expect(item.key == order[order.len - idx]);
//     }
//     try expect(maxQueue.count() == 0);

//     // re-build minQueue and maxQueue again
//     try minQueue.meldIntoSelfFromIndexable(&array_of_pairs2);
//     try maxQueue.meldIntoSelf(&array_of_pairs);

//     // create new Queues
//     // here we explore a built-in method exportAsIterable()
//     var thirdQueue = try minQueue.meldIntoNew(maxQueue.exportAsIterable());
//     defer thirdQueue.deinit();
//     try expect(thirdQueue.isMin());

//     var fourthQueue = try maxQueue.meldIntoNewFromIndexable(&array_of_pairs2);
//     defer fourthQueue.deinit();
//     try expect(fourthQueue.minOrMax() == .max);

//     var fifthQueue = try fourthQueue.cloneAsOpposite();
//     defer fifthQueue.deinit();
//     while (!fifthQueue.isEmpty()) {
//         try expect((try fifthQueue.pop()).key == (try thirdQueue.pop()).key);
//     }
// }

// fn compareI8(isMin: bool, a: i8, b: i8) bool {
//     while (isMin)
//         return a <= b;
//     return a >= b;
// }
// test "Queue: clone, cloneAsOpposite " {
//     var maxQueue = PieQ(i8, i32, .max, compareI8).init(std.testing.allocator);
//     defer maxQueue.deinit();

//     const numbers = [_]i8{ -8, 0, -7, 3, 1, 4, -9, 2, -11, 5, 6, -10 };
//     for (numbers) |num| try maxQueue.push(.{ .key = num, .val = num * num });
//     // got { 6, 5, 4, 3, 2, 1, 0, -7, -8, -9, -10, -11 }

//     const order = [_]i8{ 6, 5, 4, 3, 2, 1, 0, -7, -8, -9, -10, -11 };
//     for (order[0..4]) |digit| {
//         var item = try maxQueue.pop();
//         try expect(item.key == digit);
//     } // left { 2, 1, 0, -7, -8, -9, -10, -11 }

//     var maxQueue2 = try maxQueue.clone();
//     defer maxQueue2.deinit();

//     try expect(maxQueue2.count() == maxQueue.count());

//     for (order[4..8]) |digit| {
//         var item = try maxQueue2.pop();
//         try expect(item.key == digit);
//     } // left { -8, -9, -10, -11 }

//     var minQueue = try maxQueue2.cloneAsOpposite(); // got mirror { -11, -10, -9, -8 }
//     defer minQueue.deinit();

//     for (order[8..]) |_| {
//         var item = try minQueue.pop();
//         try expect(item.key == order[@intCast(try std.math.absInt(item.key))]);
//     }

//     try expect(minQueue.count() == 0);
// }

// test "Queue: getRoot, getAfterRoot, lockRoot" {
//     const BestErrorEver = error{
//         konnichiwa,
//         sayonara,
//         kochi,
//     };

//     var maxQueue = PieQ(u8, BestErrorEver, .max, compareU8).init(std.testing.allocator);
//     defer maxQueue.deinit();

//     const alphabet = [_]u8{
//         'H', 'F', 'A', 'G', 'E', 'B', 'C', 'D',
//     };
//     for (alphabet) |char| {
//         switch (@mod(72, char)) {
//             0 => try maxQueue.push(.{ .key = char, .val = BestErrorEver.konnichiwa }),
//             7 => try maxQueue.push(.{ .key = char, .val = BestErrorEver.sayonara }),
//             else => try maxQueue.push(.{ .key = char, .val = BestErrorEver.kochi }),
//         }
//         var root = try maxQueue.getRoot();
//         try expect(root.val == BestErrorEver.konnichiwa);
//     }

//     const ordered_alphabet = [_]u8{
//         'H', 'G', 'F', 'E', 'D', 'C', 'B', 'A',
//     };

//     maxQueue.lockRoot(true); // root protection is on
//     try expect(maxQueue.isRootLocked());

//     for (ordered_alphabet[1..8]) |char| {
//         var after_root = try maxQueue.getAfterRoot();
//         var item = try maxQueue.pop();

//         try expect(item.val == blk: {
//             if (char == 65) {
//                 var expected = BestErrorEver.sayonara;
//                 break :blk expected;
//             } else {
//                 var expected = BestErrorEver.kochi;
//                 break :blk expected;
//             }
//         });
//         try expect(std.meta.eql(item, after_root));

//         var root = try maxQueue.getRoot();
//         try expect(root.val == BestErrorEver.konnichiwa);
//     }
//     try expect(maxQueue.count() == 1);

//     maxQueue.lockRoot(false); // root protection is off
//     try expect(!maxQueue.isRootLocked());

//     var root = try maxQueue.pop();
//     try expect(root.val == BestErrorEver.konnichiwa);
//     try expect(maxQueue.count() == 0);
// }

// fn compareU32(isMin: bool, a: u32, b: u32) bool {
//     while (isMin)
//         return a <= b;
//     return a >= b;
// }
// test "Queue: understand changeRoot method" {
//     const ValuePlaceHolder = enum { JohnDoe };

//     var maxQueue = PieQ(u32, ValuePlaceHolder, .max, compareU32).init(std.testing.allocator);
//     defer maxQueue.deinit();

//     const powers = [_]u32{ 128, 2, 8, 4096, 4, 16, 256, 512, 1, 32, 64, 2048, 1024, 8192 };
//     for (powers) |power| {
//         try maxQueue.push(.{ .key = power, .val = .JohnDoe });
//     }
//     try expect((try maxQueue.getRoot()).key == 8192);

//     // change the root to some larger key than is present in powers
//     const new_root = .{ .key = 16384, .val = .JohnDoe };
//     try maxQueue.changeRoot(new_root);
//     try expect((try maxQueue.getRoot()).key == 16384);

//     // change the root to the Item with a zero key; since it is .max Queue,
//     // and giving the initial conditions, it will cause the Item with a zero key
//     // to fall down so some other Item with the largest key already present in the Queue might take its place.
//     const new_root2 = .{ .key = 0, .val = .JohnDoe };
//     try maxQueue.changeRoot(new_root2);
//     try expect((try maxQueue.getRoot()).key == 4096);

//     // reinsert key 16384 again!
//     // Queue will cause it to climb up and take the root's place.
//     const last = .{ .key = 16384, .val = .JohnDoe };
//     try maxQueue.push(last);
//     try expect((try maxQueue.getRoot()).key == 16384);

//     // attempt to change the root while it is locked will cause an error!
//     maxQueue.lockRoot(true);
//     try expect(maxQueue.isRootLocked());

//     const new_root3 = .{ .key = 32768, .val = .JohnDoe };
//     maxQueue.changeRoot(new_root3) catch |err| try std.testing.expectEqual(err, PieQ(u32, u32, .max, compareU32).QueueError.AttemptToModifyLockedRoot);
//     try expect((try maxQueue.getRoot()).key == 16384); // stays the same

//     // while the root is locked, pop will remove all but the root,
//     // starting from the 2nd best, then 3rd best, and so on.
//     var item = try maxQueue.pop();
//     try expect(item.key == 4096);
//     while (maxQueue.count() > 1) {
//         item = try maxQueue.pop();
//         try expect((try maxQueue.getRoot()).key == 16384);
//     }
//     try expect(item.key == 0); // the last one popped was the smallest present

//     // if nothing has left save the root, and root is still locked,
//     // all subsequent calls to pop() will return again, nothing but the root
//     try expect((try maxQueue.pop()).key == 16384);
//     try expect((try maxQueue.pop()).key == 16384);
//     try expect((try maxQueue.pop()).key == 16384);
//     try expect((try maxQueue.getRoot()).key == 16384);

//     // unlocking the root allows you to free the Queue entirely
//     maxQueue.lockRoot(false);
//     try expect(!maxQueue.isRootLocked());
//     try expect((try maxQueue.pop()).key == 16384);
//     try expect(maxQueue.count() == 0);

//     // attempt to change the root while the Queue is empty, will cause and error, too!
//     maxQueue.changeRoot(new_root3) catch |err| try std.testing.expectEqual(err, PieQ(u32, u32, .max, compareU32).QueueError.QueueIsEmpty);
// }

// const Vec = struct {
//     const T = u32;
//     vec: @Vector(4, T) = undefined,

//     const Self = @This();

//     fn prod(self: Self) T {
//         var idx: usize = 0;
//         var result: T = 1;
//         while (idx < 4) : (idx += 1) result *= self.vec[idx];
//         return result;
//     }
//     fn init(self: *Self, slice: [4]T) Self {
//         _ = self;
//         return .{ .vec = slice };
//     }
// };

// fn compareVecProd(isMin: bool, a: Vec, b: Vec) bool {
//     while (isMin)
//         return a.prod() <= b.prod();
//     return a.prod() >= b.prod();
// }
// test "Queue: use struct as key, sort vectors by their self-product" {
//     var rng = std.rand.DefaultPrng.init(0);
//     const random = rng.random();

//     var keys_ = [_][4]u32{
//         [4]u32{ 1, 2, 3, 4 },
//         [4]u32{ 5, 6, 7, 8 },
//         [4]u32{ 9, 10, 11, 12 },
//         [4]u32{ 13, 14, 15, 16 },
//         [4]u32{ 17, 18, 19, 20 },
//         [4]u32{ 21, 22, 23, 24 },
//     };

//     random.shuffle([4]u32, &keys_);

//     var maxQueue = PieQ(Vec, u32, .max, compareVecProd).init(std.testing.allocator);
//     defer maxQueue.deinit();

//     for (keys_) |item| {
//         var vec: Vec = .{};
//         vec = vec.init(item);
//         try maxQueue.push(.{ .key = vec, .val = vec.prod() });
//     }

//     try expect((try maxQueue.pop()).val == 255024);
//     try expect((try maxQueue.pop()).val == 116280);
//     try expect((try maxQueue.pop()).val == 43680);
//     try expect((try maxQueue.pop()).val == 11880);
//     try expect((try maxQueue.pop()).val == 1680);
//     try expect((try maxQueue.pop()).val == 24);
// }

// fn compareOddEven(isMin: bool, a: u32, b: u32) bool {
//     while (isMin)
//         return a % 2 == 0 and b % 2 == 1;
//     return a % 2 == 1 and b % 2 == 0;
// }

// test "Queue: sort binary keys" {
//     var rng = std.rand.DefaultPrng.init(0);
//     const random = rng.random();

//     const Value = enum { even, odd };

//     var keys = std.ArrayList(u32).init(std.testing.allocator);
//     defer keys.deinit();

//     var int: u32 = 0;
//     while (int < 16) : (int += 1) {
//         // fill the array with binary keys
//         keys.append(int % 2) catch unreachable;
//     }
//     random.shuffle(u32, keys.items); // scatter keys

//     var minQueue = PieQ(u32, Value, .min, compareOddEven).init(std.testing.allocator);
//     defer minQueue.deinit();

//     for (keys.items) |key| try minQueue.push(.{ .key = key, .val = if (key == 0) .even else .odd });

//     var maxQueue = try minQueue.cloneAsOpposite(); // reverse the Queue
//     defer maxQueue.deinit();

//     while (int > 0) : (int -= 1) {
//         var item = try minQueue.pop();
//         var item2 = try maxQueue.pop();
//         try expect(if (int > 8) item.val == .even else item.val == .odd);
//         try expect(if (int > 8) item2.val == .odd else item2.val == .even);
//     }
// }

// const Colors = enum { violet, blue, yellow, red };
// fn compareColors(isMin: bool, a: Colors, b: Colors) bool {
//     while (isMin)
//         return a == .violet or b == .yellow and a == .blue or b == .red;
//     return a == .red or b == .blue and a == .yellow or b == .violet;
// }
// test "Queue: sort enum literals" {
//     var rng = std.rand.DefaultPrng.init(0);
//     const random = rng.random();

//     var colors = std.ArrayList(Colors).init(std.testing.allocator);
//     defer colors.deinit();

//     var int: u32 = 0;
//     while (int < 4) : (int += 1) {
//         colors.append(Colors.violet) catch unreachable;
//         colors.append(Colors.blue) catch unreachable;
//         colors.append(Colors.yellow) catch unreachable;
//         colors.append(Colors.red) catch unreachable;
//     }
//     random.shuffle(Colors, colors.items); // scatter colors

//     var minQueue = PieQ(Colors, u2, .min, compareColors).init(std.testing.allocator);
//     defer minQueue.deinit();

//     for (colors.items) |color| try minQueue.push(.{ .key = color, .val = blk: {
//         var result: u2 = undefined;
//         switch (color) {
//             .violet => result = 0,
//             .blue => result = 1,
//             .yellow => result = 2,
//             .red => result = 3,
//         }
//         break :blk result;
//     } });

//     var maxQueue = try minQueue.cloneAsOpposite(); // get reverse of the Queue
//     defer maxQueue.deinit();

//     int = 4 * 4;
//     while (int > 0) : (int -= 1) {
//         var item = try minQueue.pop();
//         var item2 = try maxQueue.pop();
//         try expect(item.val + item2.val == 3); // a proof that items are sorted when popped out
//     }
// }
