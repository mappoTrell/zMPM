const std = @import("std");
const Thread = std.Thread;
const b = @import("base.zig");
const grid = @import("grid.zig");
const mp = @import("materialPoint.zig");
const mat = @import("materials.zig");
const m = @import("main.zig");

const Vec2 = b.Vec2;
const Mat22 = b.Mat22;

const Ident: Mat22 = b.initMat22(1, 0, 0, 1);
const zero22: Mat22 = b.initMat22(0, 0, 0, 0);

pub const Solver_2d = struct {
    alloc: std.mem.Allocator,
    e: m.env = undefined,
    thread_pool: Thread.Pool = undefined,
    comptime base: b.basis = .cubic_bspline_basis,

    pub fn init(e: m.env, alloc: std.mem.Allocator, comptime base: b.basis) !@This() {
        const t = @This(){ .alloc = alloc, .e = e };
        //try t.thread_pool.init(.{ .allocator = alloc, .n_jobs = 6 });
        _ = base;
        return t;
    }

    pub fn deinit(Self: *@This()) void {
        _ = Self;
        //Self.thread_pool.deinit();
    }

    pub fn run(Self: *@This()) !u64 {
        const e = Self.e;
        var t: u64 = 0;
        var it: u64 = 1;

        const c_count = std.Thread.getCpuCount() catch 5;
        std.debug.print("{}\n", .{c_count});

        var wG = Thread.WaitGroup{};
        wG.reset();

        while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {
            if (t % 50 == 0) {
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

            // reset Grid
            var i: u64 = 0;
            var incr: u64 = 0;
            incr = e.grid.gp.len / c_count;
            while (i < e.grid.gp.len) : (i += incr) {
                const end = if (i + incr <= e.grid.gp.len) i + incr else e.grid.gp.len;
                //Self.thread_pool.spawnWg(&wG, resetGrid, .{ &(Self.e.grid), i, end });
                //const thr = try Thread.spawn(.{}, resetGrid, .{ Self, &wG, i, end });
                wG.spawnManager(resetGrid, .{ Self, i, end });
                //thr.detach();
            }

            while (!wG.isDone()) {
                std.atomic.spinLoopHint();
            }
            wG.reset();

            //resetGrid(Self, &wG, 0, e.grid.gp.len);

            //points to grid
            for (e.shapes.items, 0..) |shape, shpN| {
                i = 0;
                incr = shape.len / c_count;
                while (i < shape.len) : (i += incr) {
                    const end = if (i + incr <= shape.len) (i + incr) else shape.len;
                    wG.spawnManager(pointsToGrid, .{ Self, shpN, Self.base, i, end });
                }
            }
            while (!wG.isDone()) {
                std.atomic.spinLoopHint();
            }

            // for (e.shapes.items, 0..) |shape, shpN| {
            //     pointsToGrid(Self, &wG, shpN, Self.base, 0, shape.len);
            //
            // }

            wG.reset();

            // //update Grid
            i = 0;
            incr = e.grid.gp.len / c_count;
            while (i < e.grid.gp.len) : (i += incr) {
                const end = if (i + incr <= e.grid.gp.len) i + incr else e.grid.gp.len;
                wG.spawnManager(updateGrid, .{ Self, i, end });
            }
            while (!wG.isDone()) {
                std.atomic.spinLoopHint();
            }
            wG.reset();

            //updateGrid(Self, &wG, 0, e.grid.gp.len);

            // //grid to points
            for (e.shapes.items, 0..) |shape, shpN| {
                i = 0;
                incr = shape.len / c_count;
                while (i < shape.len) : (i += incr) {
                    const end = if (i + incr <= shape.len) i + incr else shape.len;
                    wG.spawnManager(gridToPoints, .{ Self, shpN, Self.base, i, end });
                }
            }
            while (!wG.isDone()) {
                std.atomic.spinLoopHint();
            }
            wG.reset();

            // for (e.shapes.items, 0..) |shape, shpN| {
            //     gridToPoints(Self, &wG, shpN, Self.base, 0, shape.len);
            //
            //     // }
            // }
        }
        std.debug.print("Done\n", .{});
        return 1;
    }

    fn resetGrid(Self: *@This(), start: u64, end: u64) void {
        // wg.start();
        // defer wg.finish();
        const gp = Self.e.grid.gp;
        for (gp.items(.mass)[start..end], gp.items(.force)[start..end], gp.items(.mom0)[start..end]) |*mass, *force, *momentum| {
            //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });
            mass.* = 0;
            force.* = .{ 0, 0 };
            momentum.* = .{ 0, 0 };
        }
    }

    fn updateGrid(Self: *@This(), start: u64, end: u64) void {
        const e = Self.e;

        const gp = Self.e.grid.gp;
        for (gp.items(.force)[start..end], gp.items(.mom0)[start..end], gp.items(.mom)[start..end], gp.items(.fixed)[start..end]) |force, *momentum0, *momentum, fixed| {
            momentum.* = momentum0.* + b.scalar(force, e.timeStep);

            if (fixed[0]) {
                momentum.*[0] = 0;
                momentum0.*[0] = 0;
            }

            if (fixed[1]) {
                momentum.*[1] = 0;
                momentum0.*[1] = 0;
            }
        }
    }

    fn gridToPoints(Self: *@This(), shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
        const e = Self.e;
        const gp = Self.e.grid.gp;

        var bas = switch (basis) {
            .linear_basis => b.Linear_Basis{},
            .quad_bspline_basis => b.Quad_Bspline_Basis{},
            .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
        };

        const shape = e.shapes.items[shpN];

        const mps = shape.mp_slice;

        for (mps.items(.pos)[start..end], mps.items(.vel)[start..end], mps.items(.stress)[start..end], mps.items(.vol)[start..end], mps.items(.vol0)[start..end], mps.items(.strain)[start..end], mps.items(.def_grad)[start..end]) |*pos, *vel, *stress, *vol, vol0, *strain, *def_grad| {
            bas.getShapeValueGradient(pos.*, e.grid);

            var vel_grad = zero22;

            for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                //if (fV != 0) std.debug.print("{d}", .{i});
                if (gp.items(.mass)[gp_I] > 1.0e-9) {
                    const vI = b.scalar(gp.items(.mom)[gp_I], 1 / gp.items(.mass)[gp_I]);
                    vel.* += b.scalar((gp.items(.mom)[gp_I] - gp.items(.mom0)[gp_I]), (fV / gp.items(.mass)[gp_I]));
                    pos.* += b.scalar((gp.items(.mom)[gp_I]), ((fV * e.timeStep) / gp.items(.mass)[gp_I]));
                    //not switched 2nd and 3rd
                    vel_grad = b.addMat(vel_grad, b.initMat22(der[0] * vI[0], der[0] * vI[1], der[1] * vI[0], der[1] * vI[1]));
                }
            }
            strain.* = b.addMat(strain.*, b.scalar(b.addMat(vel_grad, b.transpose(vel_grad)), e.timeStep * 0.5));
            def_grad.* = b.matMult(def_grad.*, b.addMat(Ident, b.scalar(vel_grad, e.timeStep)));
            const J = b.det(def_grad.*);
            vol.* = J * vol0;

            stress.* = shape.mat.update_stress(strain.*);
        }
    }

    fn pointsToGrid(Self: *@This(), shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
        const e = Self.e;
        const gp = Self.e.grid.gp;

        var bas = switch (basis) {
            .linear_basis => b.Linear_Basis{},
            .quad_bspline_basis => b.Quad_Bspline_Basis{},
            .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
        };

        const shape = e.shapes.items[shpN];

        const mps = shape.mp_slice;

        for (mps.items(.pos)[start..end], mps.items(.mass)[start..end], mps.items(.vel)[start..end], mps.items(.stress)[start..end], mps.items(.vol)[start..end]) |pos, mass, vel, stress, vol| {
            //std.debug.print("{d}\n", .{pos});
            bas.getShapeValueGradient(pos, e.grid);

            for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                while (gp.items(.changes)[gp_I].cmpxchgWeak(false, true, .acquire, .monotonic) orelse false) {
                    std.atomic.spinLoopHint();
                }

                gp.items(.mass)[gp_I] += mass * fV;
                gp.items(.mom0)[gp_I] += b.scalar(vel, fV * mass);
                gp.items(.force)[gp_I] += b.scalar(.{
                    stress[0][0] * der[0] + stress[0][1] * der[1],
                    stress[1][0] * der[0] + stress[1][1] * der[1],
                }, -vol) + b.scalar(e.ext_acc, mass * fV);

                gp.items(.changes)[gp_I].store(false, .release);
            }
        }
    }
};

pub const Solver_2d_2 = struct {
    alloc: std.mem.Allocator,
    e: m.env = undefined,
    thread_pool: Thread.Pool = undefined,
    comptime base: b.basis = .cubic_bspline_basis,

    pub fn init(e: m.env, alloc: std.mem.Allocator, comptime base: b.basis) !@This() {
        var t = @This(){ .alloc = alloc, .e = e };
        try t.thread_pool.init(.{ .allocator = alloc });
        _ = base;
        return t;
    }

    pub fn deinit(Self: *@This()) void {
        Self.thread_pool.deinit();
    }

    pub fn run(Self: *@This()) !u64 {
        const e = Self.e;
        var t: u64 = 0;
        var it: u64 = 1;

        const c_count = std.Thread.getCpuCount() catch 5;
        std.debug.print("{}\n", .{c_count});

        var wG = Thread.WaitGroup{};
        wG.reset();

        while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {
            if (t % 50 == 0) {
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

            // reset Grid
            var i: u64 = 0;
            var incr: u64 = 0;
            incr = e.grid.gp.len / c_count;
            while (i < e.grid.gp.len) : (i += incr) {
                const end = if (i + incr <= e.grid.gp.len) i + incr else e.grid.gp.len;
                Self.thread_pool.spawnWg(&wG, resetGrid, .{ Self, i, end });
                //const thr = try Thread.spawn(.{}, resetGrid, .{ Self, &wG, i, end });
                //wG.spawnManager(resetGrid, .{ Self, i, end });
                //thr.detach();
            }

            Self.thread_pool.waitAndWork(&wG);

            // while (!wG.isDone()) {
            //     std.atomic.spinLoopHint();
            // }
            //wG.reset();

            //resetGrid(Self, &wG, 0, e.grid.gp.len);

            //points to grid
            for (e.shapes.items, 0..) |shape, shpN| {
                i = 0;
                incr = shape.len / c_count;
                while (i < shape.len) : (i += incr) {
                    const end = if (i + incr <= shape.len) (i + incr) else shape.len;
                    Self.thread_pool.spawnWg(&wG, pointsToGrid, .{ Self, shpN, Self.base, i, end });
                }
            }
            Self.thread_pool.waitAndWork(&wG);

            // for (e.shapes.items, 0..) |shape, shpN| {
            //     pointsToGrid(Self, &wG, shpN, Self.base, 0, shape.len);
            //
            // }

            //wG.reset();

            // //update Grid
            i = 0;
            incr = e.grid.gp.len / c_count;
            while (i < e.grid.gp.len) : (i += incr) {
                const end = if (i + incr <= e.grid.gp.len) i + incr else e.grid.gp.len;
                Self.thread_pool.spawnWg(&wG, updateGrid, .{ Self, i, end });
            }
            Self.thread_pool.waitAndWork(&wG);
            //wG.reset();

            //updateGrid(Self, &wG, 0, e.grid.gp.len);

            // //grid to points
            for (e.shapes.items, 0..) |shape, shpN| {
                i = 0;
                incr = shape.len / c_count;
                while (i < shape.len) : (i += incr) {
                    const end = if (i + incr <= shape.len) i + incr else shape.len;
                    Self.thread_pool.spawnWg(&wG, gridToPoints, .{ Self, shpN, Self.base, i, end });
                }
            }
            Self.thread_pool.waitAndWork(&wG);
            //wG.reset();

            // for (e.shapes.items, 0..) |shape, shpN| {
            //     gridToPoints(Self, &wG, shpN, Self.base, 0, shape.len);
            //
            //     // }
            // }
        }
        std.debug.print("Done\n", .{});
        return 1;
    }

    pub fn resetGrid(Self: *@This(), start: u64, end: u64) void {
        // wg.start();
        // defer wg.finish();
        const gp = Self.e.grid.gp;
        for (gp.items(.mass)[start..end], gp.items(.force)[start..end], gp.items(.mom0)[start..end]) |*mass, *force, *momentum| {
            //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });
            mass.* = 0;
            force.* = .{ 0, 0 };
            momentum.* = .{ 0, 0 };
        }
    }

    pub fn updateGrid(Self: *@This(), start: u64, end: u64) void {
        const e = Self.e;

        const gp = Self.e.grid.gp;
        for (gp.items(.force)[start..end], gp.items(.mom0)[start..end], gp.items(.mom)[start..end], gp.items(.fixed)[start..end]) |force, *momentum0, *momentum, fixed| {
            momentum.* = momentum0.* + b.scalar(force, e.timeStep);

            if (fixed[0]) {
                momentum.*[0] = 0;
                momentum0.*[0] = 0;
            }

            if (fixed[1]) {
                momentum.*[1] = 0;
                momentum0.*[1] = 0;
            }
        }
    }

    pub fn gridToPoints(Self: *@This(), shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
        const e = Self.e;
        const gp = Self.e.grid.gp;

        var bas = switch (basis) {
            .linear_basis => b.Linear_Basis{},
            .quad_bspline_basis => b.Quad_Bspline_Basis{},
            .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
        };

        const shape = e.shapes.items[shpN];

        const mps = shape.mp_slice;

        for (mps.items(.pos)[start..end], mps.items(.vel)[start..end], mps.items(.stress)[start..end], mps.items(.vol)[start..end], mps.items(.vol0)[start..end], mps.items(.strain)[start..end], mps.items(.def_grad)[start..end]) |*pos, *vel, *stress, *vol, vol0, *strain, *def_grad| {
            bas.getShapeValueGradient(pos.*, e.grid);

            var vel_grad = zero22;

            for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                //if (fV != 0) std.debug.print("{d}", .{i});
                if (gp.items(.mass)[gp_I] > 1.0e-9) {
                    const vI = b.scalar(gp.items(.mom)[gp_I], 1 / gp.items(.mass)[gp_I]);
                    vel.* += b.scalar((gp.items(.mom)[gp_I] - gp.items(.mom0)[gp_I]), (fV / gp.items(.mass)[gp_I]));
                    pos.* += b.scalar((gp.items(.mom)[gp_I]), ((fV * e.timeStep) / gp.items(.mass)[gp_I]));
                    //not switched 2nd and 3rd
                    vel_grad = b.addMat(vel_grad, b.initMat22(der[0] * vI[0], der[0] * vI[1], der[1] * vI[0], der[1] * vI[1]));
                }
            }
            strain.* = b.addMat(strain.*, b.scalar(b.addMat(vel_grad, b.transpose(vel_grad)), e.timeStep * 0.5));
            def_grad.* = b.matMult(def_grad.*, b.addMat(Ident, b.scalar(vel_grad, e.timeStep)));
            const J = b.det(def_grad.*);
            vol.* = J * vol0;

            stress.* = shape.mat.update_stress(strain.*);
        }
    }

    pub fn pointsToGrid(Self: *@This(), shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
        const e = Self.e;
        const gp = Self.e.grid.gp;

        var bas = switch (basis) {
            .linear_basis => b.Linear_Basis{},
            .quad_bspline_basis => b.Quad_Bspline_Basis{},
            .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
        };

        const shape = e.shapes.items[shpN];

        const mps = shape.mp_slice;

        for (mps.items(.pos)[start..end], mps.items(.mass)[start..end], mps.items(.vel)[start..end], mps.items(.stress)[start..end], mps.items(.vol)[start..end]) |pos, mass, vel, stress, vol| {
            //std.debug.print("{d}\n", .{pos});
            bas.getShapeValueGradient(pos, e.grid);

            for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
                while (gp.items(.changes)[gp_I].cmpxchgWeak(false, true, .acquire, .monotonic) orelse false) {
                    std.atomic.spinLoopHint();
                }

                gp.items(.mass)[gp_I] += mass * fV;
                gp.items(.mom0)[gp_I] += b.scalar(vel, fV * mass);
                gp.items(.force)[gp_I] += b.scalar(.{
                    stress[0][0] * der[0] + stress[0][1] * der[1],
                    stress[1][0] * der[0] + stress[1][1] * der[1],
                }, -vol) + b.scalar(e.ext_acc, mass * fV);

                gp.items(.changes)[gp_I].store(false, .release);
            }
        }
    }
};

// fn solve_2d_Classic_3(e: env, comptime basis: b.basis) !void {
//     var t: u64 = 0;
//     var it: u64 = 1;

//     var bas = switch (basis) {
//         .linear_basis => b.Linear_Basis{},
//         .quad_bspline_basis => b.Quad_Bspline_Basis{},
//         .cubic_bspline_basis => b.Cubic_Bspline_Basis{},
//     };

//     @setFloatMode(std.builtin.FloatMode.optimized);

//     // _ = basis;
//     // var bas = b.Quad_Bspline_Basis{};

//     // if (basis == .linear_basis_2) {
//     //     std.debug.print("Test:{}", .{bas.testFeld});
//     // }

//     const gp = e.grid.grid_points.slice();
//     std.debug.print("test\n", .{});

//     var ts = @Vector(4, u64){ 0, 0, 0, 0 };

//     var timer = try Timer.start();
//     var vel_grad: Mat22 = undefined;

//     while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {
//         if (t % 25 == 0) {
//             std.debug.print("i:{d}\n", .{it});
//             const fP = try std.fmt.allocPrint(e.alloc, "{d}.csv", .{it});
//             defer e.alloc.free(fP);
//             const file = try e.out_Folder.createFile(fP, .{ .read = true });
//             defer file.close();

//             var buf = std.io.bufferedWriter(file.writer());
//             var w = buf.writer();

//             for (e.shapes.items, 0..) |shape, sN| {
//                 const mps = shape.material_points;
//                 for (mps.items(.pos), mps.items(.stress)) |pos, stress| {
//                     const vm: f64 = @sqrt(stress[0][0] * stress[0][0] + stress[1][1] * stress[1][1] - stress[0][0] * stress[1][1] + 3 * (stress[0][1] * stress[1][0]));
//                     //const vm1: f64 = std.math.sqrt2(3);
//                     try w.print("{d},{d},{d},{d},{d},0\n", .{ pos[0], pos[1], 0, vm, sN + 1 });
//                 }
//             }

//             it += 1;
//         }

//         //reset grid
//         timer.reset();
//         for (gp.items(@enumFromInt(0)), gp.items(.force), gp.items(.mom0)) |*mass, *force, *momentum| {
//             //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });

//             mass.* = 0;
//             force.* = .{ 0, 0 };
//             momentum.* = .{ 0, 0 };
//         }

//         // @prefetch(gp.items(.mass), .{ .cache = .data, .rw = .write, .locality = 3 });
//         // @prefetch(gp.items(.force), .{ .cache = .data, .rw = .write, .locality = 3 });
//         // @prefetch(gp.items(.mom0), .{ .cache = .data, .rw = .write, .locality = 3 });

//         //@memset(gp.items(.mass), 0);

//         // var k: u64 = 0;
//         // while (k < gp.len) {
//         //     if (k + 4 < gp.len) {
//         //         gp.items(.mass)[k..][0..4].* = @Vector(4, f64){ 0, 0, 0, 0 };
//         //         k += 4;
//         //     } else {
//         //         gp.items(.mass)[k] = 0;
//         //         k += 1;
//         //     }
//         // }

//         //@memset(gp.items(.force), .{ 0, 0 });

//         //@memset(gp.items(.mom0), .{ 0, 0 });

//         ts[0] += timer.lap();
//         //point to grid
//         for (e.shapes.items) |shape| {
//             const mps = shape.mp_slice;

//             for (mps.items(.pos), mps.items(.mass), mps.items(.vel), mps.items(.stress), mps.items(.vol)) |pos, mass, vel, stress, vol| {
//                 bas.getShapeValueGradient(pos, e.grid);

//                 for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
//                     gp.items(.mass)[gp_I] += mass * fV;
//                     gp.items(.mom0)[gp_I] += b.scalar2(vel, fV * mass);
//                     gp.items(.force)[gp_I] += b.scalar2(.{
//                         stress[0][0] * der[0] + stress[0][1] * der[1],
//                         stress[1][0] * der[0] + stress[1][1] * der[1],
//                     }, -vol) + b.scalar2(e.ext_acc, mass * fV);
//                 }
//             }
//         }
//         ts[1] += timer.lap();
//         //update grid
//         for (gp.items(.force), gp.items(.mom0), gp.items(.mom), gp.items(.fixed)) |force, *momentum0, *momentum, fixed| {
//             momentum.* = momentum0.* + b.scalar2(force, e.timeStep);

//             if (fixed[0]) {
//                 momentum.*[0] = 0;
//                 momentum0.*[0] = 0;
//             }

//             if (fixed[1]) {
//                 momentum.*[1] = 0;
//                 momentum0.*[1] = 0;
//             }
//         }
//         ts[2] += timer.lap();
//         //grid to particle
//         for (e.shapes.items) |*shape| {
//             const mps = shape.mp_slice;

//             for (mps.items(.pos), mps.items(.vel), mps.items(.stress), mps.items(.vol), mps.items(.vol0), mps.items(.strain), mps.items(.def_grad)) |*pos, *vel, *stress, *vol, vol0, *strain, *def_grad| {
//                 bas.getShapeValueGradient(pos.*, e.grid);

//                 vel_grad = zero22;

//                 for (bas.nP, bas.func, bas.ders) |gp_I, fV, der| {
//                     //if (fV != 0) std.debug.print("{d}", .{i});
//                     if (gp.items(.mass)[gp_I] > 1.0e-9) {
//                         const vI = b.scalar2(gp.items(.mom)[gp_I], 1 / gp.items(.mass)[gp_I]);
//                         vel.* += b.scalar2((gp.items(.mom)[gp_I] - gp.items(.mom0)[gp_I]), (fV / gp.items(.mass)[gp_I]));
//                         pos.* += b.scalar2((gp.items(.mom)[gp_I]), ((fV * e.timeStep) / gp.items(.mass)[gp_I]));
//                         //not switched 2nd and 3rd
//                         vel_grad = b.add22(vel_grad, b.initMat22(der[0] * vI[0], der[0] * vI[1], der[1] * vI[0], der[1] * vI[1]));
//                     }
//                 }
//                 strain.* = b.add22(strain.*, b.scalar22(b.add22(vel_grad, b.transpose22(vel_grad)), e.timeStep * 0.5));
//                 def_grad.* = b.matMult22(def_grad.*, b.add22(Ident, b.scalar22(vel_grad, e.timeStep)));
//                 const J = b.det22(def_grad.*);
//                 vol.* = J * vol0;

//                 stress.* = shape.mat.update_stress(strain.*);
//             }
//         }

//         ts[3] += timer.lap();
//     }

//     std.debug.print("{d}\n", .{ts / @as(@TypeOf(ts), @splat(t))});
// }
