const std = @import("std");

const b = @import("base.zig");

const Self = @This();

points: std.ArrayList(b.Vec2),
bound: b.Mat22,
alloc: std.mem.Allocator,
