const std = @import("std");
const g = @import("grid.zig");
const bs = @import("bsplines.zig");

pub const Vec2 = @Vector(2, f64);

const z2: Vec2 = .{ 0, 0 };

pub fn dot2(a: Vec2, b: Vec2) f64 {
    return @reduce(.Add, a * b);
}

pub fn scalar2(a: Vec2, s: f64) Vec2 {
    return @as(Vec2, @splat(s)) * a;
}

pub fn invert2(a: Vec2) Vec2 {
    return @as(Vec2, @splat(1)) / a;
}

pub const Mat22 = [2]Vec2;

pub fn initMat22(e1: f64, e2: f64, e3: f64, e4: f64) Mat22 {
    return .{
        .{ e1, e2 },
        .{ e3, e4 },
    };
}

pub fn add22(a: Mat22, b: Mat22) Mat22 {
    return .{
        a[0] + b[0],
        a[1] + b[1],
    };
}

pub fn scalar22(a: Mat22, s: f64) Mat22 {
    return .{
        scalar2(a[0], s),
        scalar2(a[1], s),
    };
}

pub fn transpose22(a: Mat22) Mat22 {
    return .{
        @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 0), ~@as(i32, 0) }),
        @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 1), ~@as(i32, 1) }),
    };
}

pub fn det22(a: Mat22) f64 {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

pub fn getShapeValueGradient_R(cL: Vec2, mp_p: Vec2, gp_p: Vec2) [3]f64 {
    const dist: Vec2 = mp_p - gp_p;

    //std.debug.print("{d}\n", .{dist});

    //const sValue: Vec2 = @max(Vec2{ 1, 1 } - Vec2{ @abs(dist[0]) / cL[0], @abs(dist[1]) / cL[1] }, Vec2{ 0, 0 });
    const sValue: Vec2 = Vec2{ 1, 1 } - @abs(dist) / cL;
    //std.debug.print("{d} {d}\n", .{ @abs(dist), cL });

    return .{ sValue[0] * sValue[1], -sValue[1] * std.math.sign(dist[0]) / cL[0], -sValue[0] * std.math.sign(dist[1]) / cL[1] };
}

pub fn matMult22(a: Mat22, b: Mat22) Mat22 {
    const c = transpose22(b);
    return initMat22(dot2(a[0], c[0]), dot2(a[0], c[1]), dot2(a[1], c[0]), dot2(a[1], c[1]));
}

test transpose22 {
    const c = initMat22(1, 2, 3, 4);
    try std.testing.expectEqual(initMat22(1, 3, 2, 4), transpose22(c));
}

test matMult22 {
    const c = initMat22(1, 2, 3, 4);
    try std.testing.expectEqual(initMat22(7, 10, 15, 22), matMult22(c, c));
}

pub const basis = enum { linear_basis, quad_bspline_basis, cubic_bspline_basis };

pub const Linear_Basis = struct {
    nP: [4]u64 = .{0} ** 4,
    func: [4]f64 = .{0} ** 4,
    ders: [4]Vec2 = .{z2} ** 4,

    pub fn getShapeValueGradient(Self: *Linear_Basis, mp_p: Vec2, grid: g.Grid_2d) void {
        Self.nP = grid.adjacentGridPoints(mp_p);
        inline for (Self.nP, 0..) |nP_I, i| {
            const gp_p = grid.grid_points.items(.pos)[nP_I];
            const cL = grid.lenght_cell;
            const dist: Vec2 = mp_p - gp_p;
            const sValue: Vec2 = Vec2{ 1, 1 } - @abs(dist) / cL;
            Self.func[i] = sValue[0] * sValue[1];
            Self.ders[i] = .{ -sValue[1] * std.math.sign(dist[0]) / cL[0], -sValue[0] * std.math.sign(dist[1]) / cL[1] };
        }
    }
};

pub const Quad_Bspline_Basis = struct {
    nP: [9]u64 = .{0} ** 9,
    func: [9]f64 = .{0} ** 9,
    ders: [9]Vec2 = .{z2} ** 9,

    pub fn getShapeValueGradient(Self: *Quad_Bspline_Basis, mp_p: Vec2, grid: g.Grid_2d) void {
        const p = mp_p * grid.lenght_cell_I;

        if (p[0] > 2 or p[0] < 0) std.debug.print("{}\n", .{p[0]});

        const bottom_left: [2]u64 = .{ @intFromFloat(p[0]), @intFromFloat(p[1]) };

        const remv = mp_p - Vec2{ @as(f64, @floatFromInt(bottom_left[0])), @as(f64, @floatFromInt(bottom_left[1])) };

        var np: [2][3]u64 = undefined;

        for (bottom_left, 0.., grid.nodesD) |bL, i, nodes| {
            const rem = remv[i];

            if (bL == 0) {
                np[i] = .{ bL, bL + 1, bL + 2 };
            } else if (bL == nodes - 2) {
                np[i] = .{ bL - 1, bL, bL + 1 };
            } else {
                if (rem < 0.5) {
                    np[i] = .{ bL - 1, bL, bL + 1 };
                } else {
                    np[i] = .{ bL, bL + 1, bL + 2 };
                }
            }
        }

        var index: u32 = 0;
        var nW_x: Vec2 = z2;
        var nW_y: Vec2 = z2;

        for (np[0]) |np_x| {
            const gp_x: f64 = @as(f64, @floatFromInt(np_x)) * grid.lenght_cell[0];
            const r_x: f64 = (mp_p[0] - gp_x) * grid.lenght_cell_I[0];

            if (np_x == 0 or np_x == grid.nodesD[0] - 1) {
                nW_x = bs.quad_bespline_type1(r_x);
            } else if (np_x == 1) {
                nW_x = bs.quad_bespline_type2(r_x);
            } else if (np_x == grid.nodesD[0] - 2) {
                nW_x = bs.quad_bespline_type4(r_x);
            } else {
                nW_x = bs.quad_bespline_type3(r_x);
            }
            for (np[1]) |np_y| {
                const gp_y: f64 = @as(f64, @floatFromInt(np_y)) * grid.lenght_cell[1];
                const r_y: f64 = (mp_p[1] - gp_y) * grid.lenght_cell_I[1];

                if (np_y == 0 or np_y == grid.nodesD[1] - 1) {
                    nW_y = bs.quad_bespline_type1(r_y);
                } else if (np_y == 1) {
                    nW_y = bs.quad_bespline_type2(r_y);
                } else if (np_x == grid.nodesD[1] - 2) {
                    nW_y = bs.quad_bespline_type4(r_y);
                } else {
                    nW_y = bs.quad_bespline_type3(r_y);
                }

                Self.nP[index] = grid.indexFrom2d(np_x, np_y);

                Self.func[index] = nW_x[0] * nW_y[0];

                Self.ders[index] = .{ nW_x[1] * grid.lenght_cell_I[0] * nW_y[0], nW_y[1] * grid.lenght_cell_I[1] * nW_x[0] };

                index += 1;
            }
        }
    }
};

test "Quad Bespline" {
    const all = std.testing.allocator;
    var grid = try g.Grid_2d.init(2, 2, 21, 21, all);
    defer grid.deinit();

    const sl = grid.grid_points.slice();

    try std.testing.expectEqual(.{ 1, 1 }, grid.lenght_cell);

    var base = Quad_Bspline_Basis{};
    std.debug.print("{d},\n{d},\n", .{ base.nP, base.func });

    base.getShapeValueGradient(Vec2{ 2.6, 2.4 }, grid);

    for (base.nP) |g_p| {
        std.debug.print("{d}", .{sl.items(.pos)[g_p]});
    }

    std.debug.print("\n{d},\n{d},\n", .{ base.nP, base.func });
    for (base.ders) |der| {
        std.debug.print("{d}", .{der});
    }
}

pub const Cubic_Bspline_Basis = struct {
    nP: [16]u64 = .{0} ** 16,
    func: [16]f64 = .{0} ** 16,
    ders: [16]Vec2 = .{z2} ** 16,

    pub fn getShapeValueGradient(Self: *Cubic_Bspline_Basis, mp_p: Vec2, grid: g.Grid_2d) void {
        const p = mp_p * grid.lenght_cell_I;

        const bottom_left: [2]u64 = .{ @intFromFloat(p[0]), @intFromFloat(p[1]) };

        var np: [2][4]u64 = undefined;

        for (bottom_left, 0.., grid.nodesD) |bL, i, nodes| {
            if (bL == 0) {
                np[i] = .{ 0, 1, 2, 3 };
            } else if (bL == nodes - 2) {
                np[i] = .{ bL - 2, bL - 1, bL, bL + 1 };
            } else {
                np[i] = .{ bL - 1, bL, bL + 1, bL + 2 };
            }
        }

        var index: u32 = 0;
        var nW_x: Vec2 = z2;
        var nW_y: Vec2 = z2;

        for (np[0]) |np_x| {
            const gp_x: f64 = @as(f64, @floatFromInt(np_x)) * grid.lenght_cell[0];
            const r_x: f64 = (mp_p[0] - gp_x) * grid.lenght_cell_I[0];

            if (np_x == 0 or np_x == grid.nodesD[0] - 1) {
                nW_x = bs.cubic_bespline_type1(r_x);
            } else if (np_x == 1) {
                nW_x = bs.cubic_bespline_type2(r_x);
            } else if (np_x == grid.nodesD[0] - 2) {
                nW_x = bs.cubic_bespline_type4(r_x);
            } else {
                nW_x = bs.cubic_bespline_type3(r_x);
            }
            for (np[1]) |np_y| {
                const gp_y: f64 = @as(f64, @floatFromInt(np_y)) * grid.lenght_cell[1];
                const r_y: f64 = (mp_p[1] - gp_y) * grid.lenght_cell_I[1];

                if (np_y == 0 or np_y == grid.nodesD[1] - 1) {
                    nW_y = bs.cubic_bespline_type1(r_y);
                } else if (np_y == 1) {
                    nW_y = bs.cubic_bespline_type2(r_y);
                } else if (np_x == grid.nodesD[1] - 2) {
                    nW_y = bs.cubic_bespline_type4(r_y);
                } else {
                    nW_y = bs.cubic_bespline_type3(r_y);
                }

                Self.nP[index] = grid.indexFrom2d(np_x, np_y);

                Self.func[index] = nW_x[0] * nW_y[0];

                Self.ders[index] = .{ nW_x[1] * grid.lenght_cell_I[0] * nW_y[0], nW_y[1] * grid.lenght_cell_I[1] * nW_x[0] };

                index += 1;
            }
        }
    }
};

// pub const Basis = union(enum) {
//     linear_basis: Linear_Basis,
//     linear_basis_2: Linear_Basis_2,

//     pub fn getShapeValueGradient(Self: *Basis, mp_p: Vec2, grid: g.Grid_2d) void {
//         switch (Self) {
//             inline else => |case| return case.getShapeValueGradient(mp_p, grid),
//         }
//     }
// };
