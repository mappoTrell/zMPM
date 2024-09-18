const std = @import("std");

const b = @import("base.zig");
const m = @import("materials.zig");

const Vec2 = b.Vec2;
const Mat22 = b.Mat22;

const Ident: Mat22 = b.initMat22(1, 0, 0, 1);
const zero22: Mat22 = b.initMat22(0, 0, 0, 0);
const z2: Vec2 = .{ 0, 0 };

// pub const nG = struct {
// nI: [4]u64 = .{ 0, 0, 0, 0 },
// func: [4]f64 = .{ 0, 0, 0, 0 },
// derv: [4]Vec2 = .{ z2, z2, z2, z2 },
// };

pub const MP_2d_Classic = struct {
    mass: f64 = 0,
    vol0: f64 = 0,
    vol: f64 = 0,
    pos: Vec2 = .{ 0, 0 },
    vel: Vec2 = .{ 0, 0 },
    restraint: Vec2 = .{ 0, 0 },
    def_grad: Mat22 = Ident,
    //vel_grad: Mat22 = zero22,
    strain: Mat22 = zero22,
    stress: Mat22 = zero22,
};

pub const Shape_2d_Classic = struct {
    alloc: std.mem.Allocator,
    material_points: std.MultiArrayList(MP_2d_Classic),
    mp_slice: std.MultiArrayList(MP_2d_Classic).Slice,
    len: u64,
    mat: m.Elastic_Material,

    pub fn deinit(Self: *Shape_2d_Classic) void {
        Self.material_points.deinit(Self.alloc);
    }

    pub fn createCircle(center: Vec2, r: f64, offset: f64, mat: m.Elastic_Material, vel0: Vec2, alloc: std.mem.Allocator) !Shape_2d_Classic {
        var x: f64 = -r;

        var temp = std.MultiArrayList(MP_2d_Classic){};
        errdefer temp.deinit(alloc);

        var i: u64 = 0;

        while (x <= r + r / 10) {
            var y: f64 = -r;
            while (y <= r) {
                if (x * x + y * y <= r * r) {
                    var p = MP_2d_Classic{};
                    //std.debug.print("{}\n", .{i});
                    p.pos = Vec2{ x, y } + center;
                    p.vol = offset * offset;
                    p.vel = vel0;
                    p.vol0 = p.vol;
                    p.mass = p.vol * mat.density;
                    try temp.append(alloc, p);

                    i += 1;
                }
                y += offset;
            }
            x += offset;
        }

        return .{
            .alloc = alloc,
            .material_points = temp,
            .mp_slice = temp.slice(),
            .mat = mat,
            .len = i,
        };
    }
};
