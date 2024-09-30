const std = @import("std");

const b = @import("base.zig");

const Vec2 = b.Vec2;
const Mat22 = b.Mat22;

const Ident: Mat22 = b.initMat22(1, 0, 0, 1);

const Rigid_Material = struct {
    density: f64,
};

pub const Elastic_Material = struct {
    density: f64,
    E: f64,
    nu: f64,
    lambda: f64,
    mu: f64,

    pub fn init(E: f64, nu: f64, density: f64) Elastic_Material {
        return .{
            .density = density,
            .E = E,
            .nu = nu,
            .lambda = E * nu / (1 + nu) / (1 - 2 * nu),
            .mu = E / 2 / (1 + nu),
        };
    }

    pub fn update_stress(Self: Elastic_Material, strain: Mat22) Mat22 {
        return b.add22(b.scalar22(Ident, Self.lambda * (strain[0][0] + strain[1][1])), b.scalar22(strain, 2 * Self.mu));
    }
};