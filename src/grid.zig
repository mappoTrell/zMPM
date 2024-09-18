const std = @import("std");
const b = @import("base.zig");

const Vec2 = b.Vec2;

const Grid_Point_2d = struct {
    mass: f64 = 0,
    pos: Vec2 = .{ 0, 0 },
    mom: Vec2 = .{ 0, 0 },
    mom0: Vec2 = .{ 0, 0 },
    force: Vec2 = .{ 0, 0 },
    fixed: [2]bool = .{ false, false },
    changes: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),
};

pub const Grid_2d = struct {
    alloc: std.mem.Allocator,
    lenght_grid: Vec2,
    nodesD: [2]u64,
    nodes: u64,
    lenght_cell: Vec2,
    lenght_cell_I: Vec2,
    grid_points: std.MultiArrayList(Grid_Point_2d),
    gp: std.MultiArrayList(Grid_Point_2d).Slice,

    pub fn indexFrom2d(Self: Grid_2d, x: u64, y: u64) u64 {
        const i: u64 = Self.nodesD[1] * x + y;

        std.debug.assert(i < Self.nodes);
        return i;
    }

    pub fn deinit(Self: *Grid_2d) void {
        Self.grid_points.deinit(Self.alloc);
    }

    pub fn init(g_x: f64, g_y: f64, n_x: u64, n_y: u64, alloc: std.mem.Allocator) !Grid_2d {
        var t = Grid_2d{
            .alloc = alloc,
            .lenght_grid = .{ g_x, g_y },
            .nodesD = .{ n_x, n_y },
            .nodes = ((n_x) * (n_y)),
            .lenght_cell = .{ 0, 0 },
            .lenght_cell_I = .{ 0, 0 },
            .grid_points = .{},
            .gp = undefined,
        };

        t.lenght_cell = t.lenght_grid / Vec2{ @as(f64, @floatFromInt(n_x - 1)), @as(f64, @floatFromInt(n_y - 1)) };
        t.lenght_cell_I = b.invert2(t.lenght_cell);

        try t.grid_points.setCapacity(alloc, t.nodes);

        for (0..n_x) |x| {
            for (0..n_y) |y| {
                var p = Grid_Point_2d{};
                p.pos = t.lenght_cell * Vec2{ @as(f64, @floatFromInt(x)), @as(f64, @floatFromInt(y)) };
                // if (x == 0 or x == n_x) {
                //     p.fixed[0] = true;
                // }
                // if (y == 0 or y == n_y) {
                //     p.fixed[1] = true;
                // }

                t.grid_points.insertAssumeCapacity(indexFrom2d(t, x, y), p);
            }
        }

        t.gp = t.grid_points.slice();

        return t;
    }

    pub fn adjacentGridPoints(Self: Grid_2d, pos: Vec2) [4]u64 {
        const p = pos / Self.lenght_cell;
        const x: u64 = @intFromFloat(p[0]);
        const y: u64 = @intFromFloat(p[1]);

        return .{
            Self.indexFrom2d(x, y),
            Self.indexFrom2d(x + 1, y),
            Self.indexFrom2d(x, y + 1),
            Self.indexFrom2d(x + 1, y + 1),
        };
    }
};

// pub const Grid_2d_2 = struct {
//     alloc: std.mem.Allocator,
//     lenght_grid: Vec2,
//     nodesD: [2]u64,
//     nodes: u64,
//     lenght_cell: Vec2,
//     lenght_cell_I: Vec2,
//     grid_points: std.ArrayList(Grid_Point_2d),

//     pub fn indexFrom2d(Self: Grid_2d_2, x: u64, y: u64) u64 {
//         const i: u64 = Self.nodesD[1] * x + y;

//         std.debug.assert(i < Self.nodes);
//         return i;
//     }

//     pub fn deinit(Self: *Grid_2d_2) void {
//         Self.grid_points.deinit();
//     }

//     pub fn init(g_x: f64, g_y: f64, n_x: u64, n_y: u64, alloc: std.mem.Allocator) !Grid_2d_2 {
//         var t = Grid_2d_2{
//             .alloc = alloc,
//             .lenght_grid = .{ g_x, g_y },
//             .nodesD = .{ n_x + 1, n_y + 1 },
//             .nodes = ((n_x + 1) * (n_y + 1)),
//             .lenght_cell = .{ 0, 0 },
//             .lenght_cell_I = .{ 0, 0 },
//             .grid_points = std.ArrayList(Grid_Point_2d).init(alloc),
//         };

//         t.lenght_cell = t.lenght_grid / Vec2{ @as(f64, @floatFromInt(n_x)), @as(f64, @floatFromInt(n_y)) };
//         t.lenght_cell_I = b.invert2(t.lenght_cell);

//         try t.grid_points.ensureTotalCapacity(t.nodes);

//         for (0..n_x + 1) |x| {
//             for (0..n_y + 1) |y| {
//                 var p = Grid_Point_2d{};
//                 p.pos = t.lenght_cell * Vec2{ @as(f64, @floatFromInt(x)), @as(f64, @floatFromInt(y)) };
//                 // if (x == 0 or x == n_x) {
//                 //     p.fixed[0] = true;
//                 // }
//                 // if (y == 0 or y == n_y) {
//                 //     p.fixed[1] = true;
//                 // }
//                 t.grid_points.insertAssumeCapacity(indexFrom2d(t, x, y), p);
//             }
//         }

//         return t;
//     }

//     pub fn adjacentGridPoints(Self: Grid_2d_2, pos: Vec2) [4]u64 {
//         const p = pos / Self.lenght_cell;
//         const x: u64 = @intFromFloat(p[0]);
//         const y: u64 = @intFromFloat(p[1]);

//         return .{
//             Self.indexFrom2d(x, y),
//             Self.indexFrom2d(x + 1, y),
//             Self.indexFrom2d(x, y + 1),
//             Self.indexFrom2d(x + 1, y + 1),
//         };
//     }
// };

test "Grid_2d" {
    const all = std.testing.allocator;
    var grid = try Grid_2d.init(1, 1, 3, 3, all);
    defer grid.deinit();

    try std.testing.expectEqual(.{ 0.5, 0.5 }, grid.lenght_cell);
    try std.testing.expectEqual(.{ 1.0, 1.0 }, grid.grid_points.items(.pos)[grid.grid_points.len - 1]);
    try std.testing.expectEqual(.{ 0.5, 0.5 }, grid.grid_points.items(.pos)[grid.indexFrom2d(1, 1)]);

    const sl = grid.grid_points.slice();

    // var changes = std.atomic.Value(bool).init(false);

    // if(changes.cmpxchgStrong(expected_value: T, new_value: T, comptime success_order: AtomicOrder, comptime fail_order: AtomicOrder))

    const tes = grid.adjacentGridPoints(Vec2{ 0.2, 0.2 });
    try std.testing.expectEqual(Vec2{ 0.0, 0.0 }, sl.items(.pos)[tes[0]]);
    try std.testing.expectEqual(Vec2{ 0.5, 0.5 }, sl.items(.pos)[tes[3]]);
}
