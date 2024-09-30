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

pub fn run(e: *m.env, t_p: *Thread.Pool, comptime base: b.basis) !u64 {
    var t: u64 = 0;
    var it: u64 = 1;

    const c_count = std.Thread.getCpuCount() catch 5;
    std.debug.print("{}\n", .{c_count});

    var wG = Thread.WaitGroup{};
    wG.reset();

    while (@as(f64, @floatFromInt(t)) <= e.timeEnd / e.timeStep) : (t += 1) {
        if (t % 250 == 0) {
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
            t_p.spawnWg(&wG, resetGrid, .{ e, i, end });
            //const thr = try Thread.spawn(.{}, resetGrid, .{ Self, &wG, i, end });
            //wG.spawnManager(resetGrid, .{ Self, i, end });
            //thr.detach();
        }

        t_p.waitAndWork(&wG);

        // while (!wG.isDone()) {
        //     std.atomic.spinLoopHint();
        // }
        wG.reset();

        //resetGrid(Self, &wG, 0, e.grid.gp.len);

        //points to grid
        for (e.shapes.items, 0..) |shape, shpN| {
            i = 0;
            incr = shape.len / c_count;
            while (i < shape.len) : (i += incr) {
                const end = if (i + incr <= shape.len) (i + incr) else shape.len;
                t_p.spawnWg(&wG, pointsToGrid, .{ e, shpN, base, i, end });
            }
        }
        t_p.waitAndWork(&wG);

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
            t_p.spawnWg(&wG, updateGrid, .{ e, i, end });
        }
        t_p.waitAndWork(&wG);
        wG.reset();

        //updateGrid(Self, &wG, 0, e.grid.gp.len);

        // //grid to points
        for (e.shapes.items, 0..) |shape, shpN| {
            i = 0;
            incr = shape.len / c_count;
            while (i < shape.len) : (i += incr) {
                const end = if (i + incr <= shape.len) i + incr else shape.len;
                t_p.spawnWg(&wG, gridToPoints, .{ e, shpN, base, i, end });
            }
        }
        t_p.waitAndWork(&wG);
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

pub fn resetGrid(e: *m.env, start: u64, end: u64) void {
    // wg.start();
    // defer wg.finish();
    const gp = e.grid.gp;
    for (gp.items(.mass)[start..end], gp.items(.force)[start..end], gp.items(.mom0)[start..end]) |*mass, *force, *momentum| {
        //if (mass.* > 0) std.debug.print("{}:{d}\n", .{ i, mass.* });
        mass.* = 0;
        force.* = .{ 0, 0 };
        momentum.* = .{ 0, 0 };
    }
}

pub fn updateGrid(e: *m.env, start: u64, end: u64) void {
    const gp = e.grid.gp;
    for (gp.items(.force)[start..end], gp.items(.mom0)[start..end], gp.items(.mom)[start..end], gp.items(.fixed)[start..end]) |force, *momentum0, *momentum, fixed| {
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
}

pub fn gridToPoints(e: *m.env, shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
    const gp = e.grid.gp;

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

pub fn pointsToGrid(e: *m.env, shpN: usize, comptime basis: b.basis, start: u64, end: u64) void {
    const gp = e.grid.gp;

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
            gp.items(.mom0)[gp_I] += b.scalar2(vel, fV * mass);
            gp.items(.force)[gp_I] += b.scalar2(.{
                stress[0][0] * der[0] + stress[0][1] * der[1],
                stress[1][0] * der[0] + stress[1][1] * der[1],
            }, -vol) + b.scalar2(e.ext_acc, mass * fV);

            gp.items(.changes)[gp_I].store(false, .release);
        }
    }
}