const std = @import("std");
const g = @import("grid.zig");
const bs = @import("bsplines.zig");

pub const Vec2 = @Vector(2, f64);

pub const Vec3 = @Vector(3, f64);

const z2: Vec2 = .{ 0, 0 };

pub fn Vec(comptime i: u8) type {
    return @Vector(i, f64);
}

pub fn dot2_old(a: Vec2, b: Vec2) f64 {
    return dot(a, b);
}

pub fn dot(a: anytype, b: anytype) f64 {
    const T = @TypeOf(a, b);
    switch (@typeInfo(T)) {
        .Vector => {
            return @reduce(.Add, a * b);
        },
        else => @compileError("no array of vectors"),
    }
}

pub fn scalar2_old(a: Vec2, s: f64) Vec2 {
    return scalar(a, s);
}

pub fn scalar(a: anytype, s: f64) @TypeOf(a) {
    const T = @TypeOf(a);
    switch (@typeInfo(T)) {
        .Vector => {
            return @as(@TypeOf(a), @splat(s)) * a;
        },
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => {
                    var res: T = undefined;
                    inline for (0..arr.len) |i| {
                        res[i] = scalar(a[i], s);
                    }
                    return res;
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => {
            @compileError("type no vector or array of Vector");
        },
    }
}

test "scalar" {
    var x = Vec2{ 1, 1 };
    var z = initMat22(1, 1, 1, 1);
    x = scalar(x, 3);
    z = scalar(z, 3);

    //const y = [4]Vec2;
    //std.debug.print("{any}\n\n", .{z});
    //std.debug.print("{any}\n\n", .{@TypeOf(x, z)});
}

pub fn invert2_old(a: Vec2) Vec2 {
    return invert(a);
}

pub fn invert(a: anytype) @TypeOf(a) {
    const T = @TypeOf(a);
    switch (@typeInfo(T)) {
        .Vector => {
            return @as(Vec2, @splat(1)) / a;
        },
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => {
                    var res: T = undefined;
                    inline for (0..arr.len) |i| {
                        res[i] = invert(a[i]);
                    }
                    return res;
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => {
            @compileError("type no vector or array of Vector");
        },
    }
}

pub const Mat22 = [2]Vec2;

pub fn initMat22(e1: f64, e2: f64, e3: f64, e4: f64) Mat22 {
    return .{
        .{ e1, e2 },
        .{ e3, e4 },
    };
}

pub fn addMat_old(a: Mat22, b: Mat22) Mat22 {
    return addMat(a, b);
}

pub fn addMat(a: anytype, b: anytype) @TypeOf(a, b) {
    const T = @TypeOf(a, b);
    switch (@typeInfo(T)) {
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => {
                    var res: T = undefined;
                    inline for (0..arr.len) |i| {
                        res[i] = a[i] + b[i];
                    }
                    return res;
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => @compileError("no array of vectors"),
    }
}

test "addMat" {
    const z = initMat22(1, 1, 1, 1);
    try std.testing.expectEqual(initMat22(2, 2, 2, 2), addMat(z, z));
}

pub fn scalar22_old(a: Mat22, s: f64) Mat22 {
    return scalar(a, s);
}

pub fn transpose(a: anytype) @TypeOf(a) {
    const T = @TypeOf(a);
    switch (@typeInfo(T)) {
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => |vec_t| {
                    if (vec_t.len == arr.len) {
                        var res: T = undefined;

                        switch (arr.len) {
                            0 => @compileError("matrix cannot have n = 0"),
                            1 => res = a,
                            2 => {
                                res = .{
                                    @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 0), ~@as(i32, 0) }),
                                    @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 1), ~@as(i32, 1) }),
                                };
                            },
                            3 => {
                                const temp0 = @shuffle(f64, a[0], a[1], [4]i32{ 0, 1, -1, -2 });
                                const temp1 = @shuffle(f64, a[0], a[1], [4]i32{ 2, -3, 2, -3 });

                                res[0] = @shuffle(f64, temp0, a[2], [3]i32{ 0, 2, -1 });
                                res[1] = @shuffle(f64, temp0, a[2], [3]i32{ 1, 3, -2 });
                                res[2] = @shuffle(f64, temp1, a[2], [3]i32{ 0, 1, -3 });
                            },
                            else => @compileError("not yet implemented"),
                        }

                        return res;
                    } else @compileError("matrix not square");
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => @compileError("no array of vectors"),
    }
}

test "transpose" {
    const x = initMat22(1, 2, 3, 4);
    try std.testing.expectEqual(initMat22(1, 3, 2, 4), transpose(x));

    const y = [3]@Vector(3, f64){
        .{ 1, 2, 3 },
        .{ 4, 5, 6 },
        .{ 7, 8, 9 },
    };
    const y_t = [3]@Vector(3, f64){
        .{ 1, 4, 7 },
        .{ 2, 5, 8 },
        .{ 3, 6, 9 },
    };
    try std.testing.expectEqual(y_t, transpose(y));
}

pub fn transpose22_old(a: Mat22) Mat22 {
    return .{
        @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 0), ~@as(i32, 0) }),
        @shuffle(f64, a[0], a[1], @Vector(2, i32){ @as(i32, 1), ~@as(i32, 1) }),
    };
}

pub fn det22_old(a: Mat22) f64 {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

pub fn det(a: anytype) f64 {
    const T = @TypeOf(a);
    switch (@typeInfo(T)) {
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => |vec_t| {
                    if (vec_t.len == arr.len) {
                        var res: f64 = 0;

                        switch (arr.len) {
                            0 => @compileError("matrix cannot have n = 0"),
                            1 => res = a[0],
                            2 => {
                                res = a[0][0] * a[1][1] - a[0][1] * a[1][0];
                            },
                            3 => {
                                var sub: [vec_t.len - 1]Vec(vec_t.len - 1) = undefined;

                                inline for (0..vec_t.len) |i| {
                                    var sub_i: u8 = 0;
                                    inline for (0..vec_t.len) |j| {
                                        if (i != j) {
                                            sub[sub_i] = @as(Vec(arr.len - 1), @as([vec_t.len]f64, a[j])[1..].*);
                                            sub_i += 1;
                                        }
                                    }
                                    res += std.math.pow(f64, -1, 2 + i) * a[i][0] * det(sub);
                                }
                            },
                            else => @compileError("not yet implemented"),
                        }

                        return res;
                    } else @compileError("matrix not square");
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => @compileError("no array of vectors"),
    }
}

test "det" {
    const x = initMat22(2, 4, 3, 1);
    try std.testing.expectEqual(-10, det(x));

    const y = [3]@Vector(3, f64){
        .{ 2, 1, 4 },
        .{ -3, 0, 2 },
        .{ 1, 0, 1 },
    };

    try std.testing.expectEqual(5, det(y));
}

pub fn getShapeValueGradient_R(cL: Vec2, mp_p: Vec2, gp_p: Vec2) [3]f64 {
    const dist: Vec2 = mp_p - gp_p;

    //std.debug.print("{d}\n", .{dist});

    //const sValue: Vec2 = @max(Vec2{ 1, 1 } - Vec2{ @abs(dist[0]) / cL[0], @abs(dist[1]) / cL[1] }, Vec2{ 0, 0 });
    const sValue: Vec2 = Vec2{ 1, 1 } - @abs(dist) / cL;
    //std.debug.print("{d} {d}\n", .{ @abs(dist), cL });

    return .{ sValue[0] * sValue[1], -sValue[1] * std.math.sign(dist[0]) / cL[0], -sValue[0] * std.math.sign(dist[1]) / cL[1] };
}

pub fn matMult22_old(a: Mat22, b: Mat22) Mat22 {
    const c = transpose(b);
    return initMat22(dot(a[0], c[0]), dot(a[0], c[1]), dot(a[1], c[0]), dot(a[1], c[1]));
}

pub fn matMult(a: anytype, b: anytype) @TypeOf(a, b) {
    const T = @TypeOf(a, b);
    switch (@typeInfo(T)) {
        .Array => |arr| {
            switch (@typeInfo(arr.child)) {
                .Vector => {
                    var res: T = undefined;
                    const bT = transpose(b);
                    inline for (&res, a) |*r_r, a_r| {
                        inline for (bT, 0..) |b_r, i| {
                            r_r[i] = dot(a_r, b_r);
                        }
                    }
                    return res;
                },
                else => {
                    //@compileLog(.{(@TypeOf(a))});
                    @compileError("no array of vectors");
                },
            }
        },
        else => @compileError("no array of vectors"),
    }
}

test transpose {
    const c = initMat22(1, 2, 3, 4);
    try std.testing.expectEqual(initMat22(1, 3, 2, 4), transpose(c));
}

test "matMult" {
    const c = initMat22(1, 2, 3, 4);
    try std.testing.expectEqual(initMat22(7, 10, 15, 22), matMult(c, c));
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

///doesn't work
pub const Quad_Bspline_Basis = struct {
    nP: [9]u64 = .{0} ** 9,
    func: [9]f64 = .{0} ** 9,
    ders: [9]Vec2 = .{z2} ** 9,

    pub fn getShapeValueGradient(Self: *Quad_Bspline_Basis, mp_p: Vec2, grid: g.Grid_2d) void {
        const p = mp_p * grid.lenght_cell_I;

        //if (p[0] > 2 or p[0] < 0) std.debug.print("{}\n", .{p[0]});

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
    var grid = try g.Grid_2d.init(4, 4, 21, 21, all);
    defer grid.deinit();

    const my_log = std.log.scoped(.test_quad);

    const sl = grid.grid_points.slice();

    //try std.testing.expectEqual(.{ 1, 1 }, grid.lenght_cell);

    var base = Quad_Bspline_Basis{};
    my_log.debug("{d},\n{d},\n", .{ base.nP, base.func });

    base.getShapeValueGradient(Vec2{ 2.6, 2.4 }, grid);

    for (base.nP) |g_p| {
        my_log.debug("{d}", .{sl.items(.pos)[g_p]});
    }

    std.log.debug("\n{d},\n{d},\n", .{ base.nP, base.func });
    for (base.ders) |der| {
        my_log.debug("{d}", .{der});
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
