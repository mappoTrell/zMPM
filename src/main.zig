const std = @import("std");
const Thread = std.Thread;
const b = @import("base.zig");
const grid = @import("grid.zig");
const mp = @import("materialPoint.zig");
const mat = @import("materials.zig");
const s = @import("solver.zig");
const s2 = @import("solver2.zig");
const s3 = @import("solver3.zig");
const ThreadPool = @import("ThreadPool.zig");

///usr/lib/linux-tools/5.15.0-119-generic/perf
const time = std.time;
const Instant = time.Instant;
const Timer = time.Timer;

const Vec2 = b.Vec2;
const Mat22 = b.Mat22;

const Ident: Mat22 = b.initMat22(1, 0, 0, 1);
const zero22: Mat22 = b.initMat22(0, 0, 0, 0);

pub const env = struct {
    alloc: std.mem.Allocator,
    grid: grid.Grid_2d,
    shapes: std.ArrayList(mp.Shape_2d_Classic),
    timeStep: f64,
    timeEnd: f64,
    ext_acc: Vec2 = .{ 0, 0 },
    out_Folder: std.fs.Dir,
};

const env_2 = struct {
    alloc: std.mem.Allocator,
    grid: grid.Grid_2d_2,
    shapes: std.ArrayList(mp.Shape_2d_Classic_2),
    timeStep: f64,
    timeEnd: f64,
    ext_acc: Vec2 = .{ 0, 0 },
    out_Folder: std.fs.Dir,
};

pub fn main() !void {
    std.debug.print("{d}\n", .{std.simd.suggestVectorLength(f64) orelse 0});

    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    const cwd = std.fs.cwd();
    cwd.makeDir("out") catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };
    const out_dir = try cwd.openDir("out", .{});

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    defer {
        if (gpa.detectLeaks()) std.debug.print("memoryleak", .{});
        _ = gpa.deinit();
    }

    var timer = try Timer.start();
    try test_2(out_dir, alloc);
    var elapsed2: f64 = @floatFromInt(timer.read());
    std.debug.print("ma: {d:.3}ms\n", .{
        elapsed2 / time.ns_per_ms,
    });

    elapsed2 = 0;

    // timer.reset();
    // try test_2(out_dir, alloc);
    // elapsed2 = @floatFromInt(timer.read());
    // std.debug.print("a: {d:.3}ms\n", .{
    //     elapsed2 / time.ns_per_ms,
    // });
}

fn test_1(out_dir: std.fs.Dir, alloc: std.mem.Allocator) !void {
    var g1 = try grid.Grid_2d.init(2, 2, 300, 300, alloc);
    defer g1.deinit();

    const mat1 = mat.Elastic_Material.init(1000, 0.1, 1000);
    const mat2 = mat.Elastic_Material.init(1000, 0.1, 1000);

    const offset = g1.lenght_cell[0] / 3;

    var shp1 = try mp.Shape_2d_Classic.createCircle(.{ 0.8, 0.8 }, 0.2, offset, mat1, .{ 0.1, 0.1 }, alloc);
    defer shp1.deinit();

    var shp2 = try mp.Shape_2d_Classic.createCircle(.{ 1.2, 1.2 }, 0.2, offset, mat2, .{ -0.1, -0.1 }, alloc);
    defer shp2.deinit();

    var l = std.ArrayList(mp.Shape_2d_Classic).init(alloc);
    defer l.deinit();

    try l.append(shp1);
    try l.append(shp2);

    var env1 = env{
        .alloc = alloc,
        .grid = g1,
        .shapes = l,
        .timeEnd = 1,
        .timeStep = 0.0001,
        .out_Folder = out_dir,
        .ext_acc = .{ 0, 0 },
    };

    var thread_safe_alloc: std.heap.ThreadSafeAllocator = .{
        .child_allocator = alloc,
    };
    const ts_alloc = thread_safe_alloc.allocator();

    _ = ts_alloc;

    // var solv = try s.Solver_2d_2.init(env1, alloc, .cubic_bspline_basis);
    // defer {
    //     solv.deinit();
    // }

    // //try solve_2d_Classic_3(env1, .cubic_bspline_basis);

    // const y = try solv.run();

    var t_p: std.Thread.Pool = undefined;
    try t_p.init(.{ .allocator = alloc });
    defer t_p.deinit();

    const y = try s2.run(&env1, &t_p, .cubic_bspline_basis);

    std.log.info("Hä {}\n ", .{y});
}

fn test_2(out_dir: std.fs.Dir, alloc: std.mem.Allocator) !void {
    var g1 = try grid.Grid_2d.init(2, 2, 300, 300, alloc);
    defer g1.deinit();

    const mat1 = mat.Elastic_Material.init(1000, 0.1, 1000);
    const mat2 = mat.Elastic_Material.init(1000, 0.1, 1000);

    const offset = g1.lenght_cell[0] / 3;

    var shp1 = try mp.Shape_2d_Classic.createCircle(.{ 0.8, 0.8 }, 0.2, offset, mat1, .{ 0.1, 0.1 }, alloc);
    defer shp1.deinit();

    var shp2 = try mp.Shape_2d_Classic.createCircle(.{ 1.2, 1.2 }, 0.2, offset, mat2, .{ -0.1, -0.1 }, alloc);
    defer shp2.deinit();

    var l = std.ArrayList(mp.Shape_2d_Classic).init(alloc);
    defer l.deinit();

    try l.append(shp1);
    try l.append(shp2);

    var env1 = env{
        .alloc = alloc,
        .grid = g1,
        .shapes = l,
        .timeEnd = 1,
        .timeStep = 0.0001,
        .out_Folder = out_dir,
        .ext_acc = .{ 0, 0 },
    };

    var thread_safe_alloc: std.heap.ThreadSafeAllocator = .{
        .child_allocator = alloc,
    };
    const ts_alloc = thread_safe_alloc.allocator();

    // var solv = try s.Solver_2d_2.init(env1, alloc, .cubic_bspline_basis);
    // defer {
    //     solv.deinit();
    // }

    // //try solve_2d_Classic_3(env1, .cubic_bspline_basis);

    // const y = try solv.run();

    var tp = ThreadPool.init(.{ .max_threads = 24 });
    defer tp.deinit();
    defer tp.shutdown();

    const y = try s3.run(&env1, &tp, .cubic_bspline_basis, ts_alloc);

    std.log.info("Hä {}\n ", .{y});
}

//usl
fn solve_2d_Classic(e: s.env) !void {
    var t: u64 = 0;
    var it: u64 = 1;

    const gp = e.grid.grid_points.slice();
    std.debug.print("test\n", .{});

    var ts = @Vector(4, u64){ 0, 0, 0, 0 };

    var timer = try Timer.start();

    while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {

        //reset grid
        timer.reset();
        // for (gp.items(@enumFromInt(0)), gp.items(.force), gp.items(.mom0)) |*mass, *force, *momentum| {
        //     //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });

        //     mass.* = 0;
        //     force.* = .{ 0, 0 };
        //     momentum.* = .{ 0, 0 };
        // }

        // @prefetch(gp.items(.mass), .{ .cache = .data, .rw = .write, .locality = 3 });
        // @prefetch(gp.items(.force), .{ .cache = .data, .rw = .write, .locality = 3 });
        // @prefetch(gp.items(.mom0), .{ .cache = .data, .rw = .write, .locality = 3 });

        @memset(gp.items(.mass), 0);

        // var k: u64 = 0;
        // while (k < gp.len) {
        //     if (k + 4 < gp.len) {
        //         gp.items(.mass)[k..][0..4].* = @Vector(4, f64){ 0, 0, 0, 0 };
        //         k += 4;
        //     } else {
        //         gp.items(.mass)[k] = 0;
        //         k += 1;
        //     }
        // }

        @memset(gp.items(.force), .{ 0, 0 });

        @memset(gp.items(.mom0), .{ 0, 0 });

        ts[0] += timer.lap();
        //point to grid
        for (e.shapes.items) |shape| {
            //std.debug.print("{}\n", .{shape.len});
            // @prefetch(shape.mp_slice.items(.pos), .{ .cache = .data, .rw = .read, .locality = 3 });
            // @prefetch(shape.mp_slice.items(.mass), .{ .cache = .data, .rw = .read, .locality = 3 });
            // @prefetch(shape.mp_slice.items(.stress), .{ .cache = .data, .rw = .read, .locality = 3 });

            for (0..shape.len) |ip| {
                const p = shape.mp_slice.get(ip);
                const near_points_I = e.grid.adjacentGridPoints(p.pos);
                //if (ip == 1) std.debug.print("mom:{}\n", .{p.mom});

                inline for (near_points_I) |gp_I| {
                    const shapeVal = b.getShapeValueGradient(e.grid.lenght_cell, p.pos, gp.items(.pos)[gp_I]);
                    gp.items(.mass)[gp_I] += p.mass * shapeVal[0];
                    gp.items(.mom0)[gp_I] += b.scalar2(p.vel, shapeVal[0] * p.mass);
                    gp.items(.force)[gp_I] += b.scalar2(.{
                        p.stress[0][0] * shapeVal[1] + p.stress[0][1] * shapeVal[2],
                        p.stress[1][0] * shapeVal[1] + p.stress[1][1] * shapeVal[2],
                    }, -p.vol) + b.scalar2(e.ext_acc, p.mass * shapeVal[0]);
                }
            }
        }
        ts[1] += timer.lap();
        //update grid
        for (gp.items(.force), gp.items(.mom0), gp.items(.mom), gp.items(.fixed)) |force, *momentum0, *momentum, fixed| {
            momentum.* = momentum0.* + b.scalar2(force, e.timeStep);

            if (fixed[0]) {
                momentum.*[0] = 0;
                momentum0.*[0] = 0;
            }

            if (fixed[1]) {
                momentum.*[1] = 0;
                momentum0.*[1] = 0;
            }
        }
        ts[2] += timer.lap();
        //grid to particle
        for (e.shapes.items) |*shape| {
            for (0..shape.len) |ip| {
                var p = shape.mp_slice.get(ip);
                const near_points_I = e.grid.adjacentGridPoints(p.pos);
                var vel_grad = zero22;
                inline for (near_points_I) |gp_I| {
                    const shapeVal = b.getShapeValueGradient(e.grid.lenght_cell, p.pos, gp.items(.pos)[gp_I]);

                    //if (shapeVal[0] != 0) std.debug.print("{d}", .{p.pos});
                    if (gp.items(.mass)[gp_I] > 1.0e-8) {
                        const vI = b.scalar2(gp.items(.mom)[gp_I], 1 / gp.items(.mass)[gp_I]);
                        p.vel += b.scalar2((gp.items(.mom)[gp_I] - gp.items(.mom0)[gp_I]), (shapeVal[0] / gp.items(.mass)[gp_I]));
                        p.pos += b.scalar2((gp.items(.mom)[gp_I]), ((shapeVal[0] * e.timeStep) / gp.items(.mass)[gp_I]));
                        //not switched 2nd and 3rd
                        vel_grad = b.add22(vel_grad, b.initMat22(shapeVal[1] * vI[0], shapeVal[1] * vI[1], shapeVal[2] * vI[0], shapeVal[2] * vI[1]));
                    }
                }

                p.strain = b.add22(p.strain, b.scalar22(b.add22(vel_grad, b.transpose22(vel_grad)), e.timeStep * 0.5));
                p.def_grad = b.matMult22(p.def_grad, b.add22(Ident, b.scalar22(vel_grad, e.timeStep)));
                const J = b.det22(p.def_grad);
                p.vol = J * p.vol0;

                p.stress = shape.mat.update_stress(p.strain);

                shape.mp_slice.set(ip, p);
            }
        }
        ts[3] += timer.lap();
        if (t % 200 == 0) {
            std.debug.print("i:{d}\n", .{it});
            // const fP = try std.fmt.allocPrint(e.alloc, "{d}.csv", .{it});
            // defer e.alloc.free(fP);
            // const file = try e.out_Folder.createFile(fP, .{ .read = true });
            // defer file.close();

            // var buf = std.io.bufferedWriter(file.writer());
            // var w = buf.writer();

            // for (e.shapes.items, 0..) |shape, sN| {
            //     const mps = shape.material_points;
            //     for (mps.items(.pos), mps.items(.stress)) |pos, stress| {
            //         const vm: f64 = std.math.sqrt(stress[0][0] * stress[0][0] + stress[1][1] * stress[1][1] - stress[0][0] * stress[1][1] + 3 * (stress[0][1] * stress[1][0]));
            //         //const vm1: f64 = std.math.sqrt2(3);
            //         try w.print("{d},{d},{d},{d},{d},0\n", .{ pos[0], pos[1], 0, vm, sN + 1 });
            //     }
            // }

            it += 1;
        }
    }

    std.debug.print("{d}\n", .{ts / @as(@TypeOf(ts), @splat(t))});
}

fn solve_2d_Classic_3(e: s.env, comptime basis: b.basis) !void {
    var t: u64 = 0;
    var it: u64 = 1;

    var bas = switch (basis) {
        .linear_basis => b.Linear_Basis{},
        .quad_bspline_basis => b.Quad_Bspline_Basis{},
        .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
    };

    @setFloatMode(std.builtin.FloatMode.optimized);

    // _ = basis;
    // var bas = b.Quad_Bspline_Basis{};

    // if (basis == .linear_basis_2) {
    //     std.debug.print("Test:{}", .{bas.testFeld});
    // }

    const gp = e.grid.grid_points.slice();
    std.debug.print("test\n", .{});

    var ts = @Vector(4, u64){ 0, 0, 0, 0 };

    var timer = try Timer.start();
    var vel_grad: Mat22 = undefined;

    while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {
        if (t % 25 == 0) {
            std.debug.print("i:{d}\n", .{it});
            const fP = try std.fmt.allocPrint(e.alloc, "{d}.csv", .{it});
            defer e.alloc.free(fP);
            const file = try e.out_Folder.createFile(fP, .{ .read = true });
            defer file.close();

            var buf = std.io.bufferedWriter(file.writer());
            var w = buf.writer();

            for (e.shapes.items, 0..) |shape, sN| {
                const mps = shape.material_points;
                for (mps.items(.pos), mps.items(.stress)) |pos, stress| {
                    const vm: f64 = @sqrt(stress[0][0] * stress[0][0] + stress[1][1] * stress[1][1] - stress[0][0] * stress[1][1] + 3 * (stress[0][1] * stress[1][0]));
                    //const vm1: f64 = std.math.sqrt2(3);
                    try w.print("{d},{d},{d},{d},{d},0\n", .{ pos[0], pos[1], 0, vm, sN + 1 });
                }
            }

            it += 1;
        }

        //reset grid
        timer.reset();
        // for (gp.items(@enumFromInt(0)), gp.items(.force), gp.items(.mom0)) |*mass, *force, *momentum| {
        //     //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });

        //     mass.* = 0;
        //     force.* = .{ 0, 0 };
        //     momentum.* = .{ 0, 0 };
        // }

        // @prefetch(gp.items(.mass), .{ .cache = .data, .rw = .write, .locality = 3 });
        // @prefetch(gp.items(.force), .{ .cache = .data, .rw = .write, .locality = 3 });
        // @prefetch(gp.items(.mom0), .{ .cache = .data, .rw = .write, .locality = 3 });

        @memset(gp.items(.mass), 0);

        // var k: u64 = 0;
        // while (k < gp.len) {
        //     if (k + 4 < gp.len) {
        //         gp.items(.mass)[k..][0..4].* = @Vector(4, f64){ 0, 0, 0, 0 };
        //         k += 4;
        //     } else {
        //         gp.items(.mass)[k] = 0;
        //         k += 1;
        //     }
        // }

        @memset(gp.items(.force), .{ 0, 0 });

        @memset(gp.items(.mom0), .{ 0, 0 });

        ts[0] += timer.lap();
        //point to grid
        for (e.shapes.items) |shape| {
            const mps = shape.mp_slice;

            for (mps.items(.pos), mps.items(.mass), mps.items(.vel), mps.items(.stress), mps.items(.vol)) |pos, mass, vel, stress, vol| {
                bas.getShapeValueGradient(pos, e.grid);

                for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                    gp.items(.mass)[gp_I] += mass * fV;
                    gp.items(.mom0)[gp_I] += b.scalar2(vel, fV * mass);
                    gp.items(.force)[gp_I] += b.scalar2(.{
                        stress[0][0] * der[0] + stress[0][1] * der[1],
                        stress[1][0] * der[0] + stress[1][1] * der[1],
                    }, -vol) + b.scalar2(e.ext_acc, mass * fV);
                }
            }
        }
        ts[1] += timer.lap();
        //update grid
        for (gp.items(.force), gp.items(.mom0), gp.items(.mom), gp.items(.fixed)) |force, *momentum0, *momentum, fixed| {
            momentum.* = momentum0.* + b.scalar2(force, e.timeStep);

            if (fixed[0]) {
                momentum.*[0] = 0;
                momentum0.*[0] = 0;
            }

            if (fixed[1]) {
                momentum.*[1] = 0;
                momentum0.*[1] = 0;
            }
        }
        ts[2] += timer.lap();
        //grid to particle
        for (e.shapes.items) |*shape| {
            const mps = shape.mp_slice;

            for (mps.items(.pos), mps.items(.vel), mps.items(.stress), mps.items(.vol), mps.items(.vol0), mps.items(.strain), mps.items(.def_grad)) |*pos, *vel, *stress, *vol, vol0, *strain, *def_grad| {
                bas.getShapeValueGradient(pos.*, e.grid);

                vel_grad = zero22;

                for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                    //if (fV != 0) std.debug.print("{d}", .{i});
                    if (gp.items(.mass)[gp_I] > 1.0e-9) {
                        const vI = b.scalar2(gp.items(.mom)[gp_I], 1 / gp.items(.mass)[gp_I]);
                        vel.* += b.scalar2((gp.items(.mom)[gp_I] - gp.items(.mom0)[gp_I]), (fV / gp.items(.mass)[gp_I]));
                        pos.* += b.scalar2((gp.items(.mom)[gp_I]), ((fV * e.timeStep) / gp.items(.mass)[gp_I]));
                        //not switched 2nd and 3rd
                        vel_grad = b.add22(vel_grad, b.initMat22(der[0] * vI[0], der[0] * vI[1], der[1] * vI[0], der[1] * vI[1]));
                    }
                }
                strain.* = b.add22(strain.*, b.scalar22(b.add22(vel_grad, b.transpose22(vel_grad)), e.timeStep * 0.5));
                def_grad.* = b.matMult22(def_grad.*, b.add22(Ident, b.scalar22(vel_grad, e.timeStep)));
                const J = b.det22(def_grad.*);
                vol.* = J * vol0;

                stress.* = shape.mat.update_stress(strain.*);
            }
        }

        ts[3] += timer.lap();
    }

    std.debug.print("{d}\n", .{ts / @as(@TypeOf(ts), @splat(t))});
}

fn solve_2d_Classic_2(e: env_2) !void {
    var t: u64 = 0;
    var it: u64 = 1;

    const gp = e.grid.grid_points;

    std.debug.print("test\n", .{});

    while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {

        //reset grid

        for (gp.items) |*grid_p| {
            //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });
            grid_p.mass = 0;
            grid_p.force = .{ 0, 0 };
            grid_p.mom0 = .{ 0, 0 };
        }

        //point to grid
        for (e.shapes.items) |shape| {
            //std.debug.print("{}\n", .{shape.len});
            for (0..shape.len) |ip| {
                const p = shape.material_points.items[ip];
                const near_points_I = e.grid.adjacentGridPoints(p.pos);
                //if (ip == 1) std.debug.print("mom:{}\n", .{p.mom});
                inline for (near_points_I) |gp_I| {
                    const shapeVal = b.getShapeValueGradient(e.grid.lenght_cell, p.pos, gp.items[gp_I].pos);
                    gp.items[gp_I].mass += p.mass * shapeVal[0];
                    gp.items[gp_I].mom0 += b.scalar2(p.vel, shapeVal[0] * p.mass);
                    gp.items[gp_I].force += b.scalar2(.{
                        p.stress[0][0] * shapeVal[1] + p.stress[0][1] * shapeVal[2],
                        p.stress[1][0] * shapeVal[1] + p.stress[1][1] * shapeVal[2],
                    }, -p.vol) + b.scalar2(e.ext_acc, p.mass * shapeVal[0]);
                }
            }
        }

        //update grid
        for (gp.items) |*g| {
            g.mom = g.mom0 + b.scalar2(g.force, e.timeStep);

            if (g.fixed[0]) {
                g.mom[0] = 0;
                g.mom0[0] = 0;
            }

            if (g.fixed[1]) {
                g.mom[1] = 0;
                g.mom0[1] = 0;
            }
        }

        //grid to particle
        for (e.shapes.items) |*shape| {
            for (shape.material_points.items) |*p| {
                const near_points_I = e.grid.adjacentGridPoints(p.pos);
                var vel_grad = zero22;
                inline for (near_points_I) |gp_I| {
                    const shapeVal = b.getShapeValueGradient(e.grid.lenght_cell, p.pos, gp.items[gp_I].pos);

                    if (gp.items[gp_I].mass > 1.0e-8) {
                        const vI = b.scalar2(gp.items[gp_I].mom, 1 / gp.items[gp_I].mass);
                        p.vel += b.scalar2((gp.items[gp_I].mom - gp.items[gp_I].mom0), (shapeVal[0] / gp.items[gp_I].mass));
                        p.pos += b.scalar2((gp.items[gp_I].mom), ((shapeVal[0] * e.timeStep) / gp.items[gp_I].mass));
                        //not switched 2nd and 3rd
                        vel_grad = b.add22(vel_grad, b.initMat22(shapeVal[1] * vI[0], shapeVal[1] * vI[1], shapeVal[2] * vI[0], shapeVal[2] * vI[1]));
                    }
                }

                p.strain = b.add22(p.strain, b.scalar22(b.add22(vel_grad, b.transpose22(vel_grad)), e.timeStep * 0.5));
                p.def_grad = b.matMult22(p.def_grad, b.add22(Ident, b.scalar22(vel_grad, e.timeStep)));
                const J = b.det22(p.def_grad);
                p.vol = J * p.vol0;

                p.stress = shape.mat.update_stress(p.strain);
            }
        }

        if (t % 200 == 0) {
            std.debug.print("i:{d}\n", .{it});
            // const fP = try std.fmt.allocPrint(e.alloc, "{d}.csv", .{it});
            // defer e.alloc.free(fP);
            // const file = try e.out_Folder.createFile(fP, .{ .read = true });
            // defer file.close();

            // var buf = std.io.bufferedWriter(file.writer());
            // var w = buf.writer();

            // for (e.shapes.items, 0..) |shape, sN| {
            //     const mps = shape.material_points;
            //     for (mps.items) |mp_t| {
            //         const stress = mp_t.stress;
            //         const pos = mp_t.pos;
            //         const vm: f64 = std.math.sqrt(stress[0][0] * stress[0][0] + stress[1][1] * stress[1][1] - stress[0][0] * stress[1][1] + 3 * (stress[0][1] * stress[1][0]));
            //         //const vm1: f64 = std.math.sqrt2(3);
            //         try w.print("{d},{d},{d},{d},{d},0\n", .{ pos[0], pos[1], 0, vm, sN + 1 });
            //     }
            // }

            it += 1;
        }
    }
}

test "simple test" {
    _ = b;
    _ = grid;
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
