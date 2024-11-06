const std = @import("std");

const b = @import("base.zig");

const Vec2 = b.Vec2;
const Mat22 = b.Mat22;

const Ident: Mat22 = b.initMat22(1, 0, 0, 1);

pub const Material = union(enum) {
    elastic_material: Elastic_Material,
    rigid_material: Rigid_Material,

    pub fn update_stress(Self: *Material, strain: Mat22) Mat22 {
        switch (Self.*) {
            inline else => |*case| return case.update_stress(strain),
        }
    }
    pub fn get_density(Self: Material) f64 {
        switch (Self) {
            inline else => |case| return case.density,
        }
    }
};

const Rigid_Material = struct {
    density: f64,

    pub fn update_stress(Self: Rigid_Material, _: Mat22) Mat22 {
        _ = Self;
        return Ident;
    }
};

pub const Elastic_Material = struct {
    density: f64,
    E: f64,
    nu: f64,
    lambda: f64,
    mu: f64,
    t1: i8 = -1,

    pub fn init(E: f64, nu: f64, density: f64) Elastic_Material {
        return .{
            .density = density,
            .E = E,
            .nu = nu,
            .lambda = E * nu / (1 + nu) / (1 - 2 * nu),
            .mu = E / 2 / (1 + nu),
        };
    }

    pub fn update_stress(Self: *Elastic_Material, strain: Mat22) Mat22 {
        Self.t1 *= -1;
        return b.addMat(b.scalar(Ident, Self.lambda * (strain[0][0] + strain[1][1])), b.scalar(strain, 2 * Self.mu));
    }
};

pub const Elasto_Plastic_Material = struct {
    density: f64,
    E: f64,
    nu: f64,
    lambda: f64,
    mu: f64,
    kappa: f64,

    pub fn init(E: f64, nu: f64, density: f64) Elasto_Plastic_Material {
        return .{
            .density = density,
            .E = E,
            .nu = nu,
            .lambda = E * nu / (1 + nu) / (1 - 2 * nu),
            .mu = E / 2 / (1 + nu),
        };
    }

    pub fn update_stress(Self: Elasto_Plastic_Material, strain: Mat22) Mat22 {
        return b.addMat(b.scalar(Ident, Self.lambda * (strain[0][0] + strain[1][1])), b.scalar(strain, 2 * Self.mu));
    }
};
