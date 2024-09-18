const std = @import("std");
const b = @import("base.zig");
const Vec2 = b.Vec2;

pub fn quad_bespline_type1(r: f64) Vec2 {
    if (r < -1.5 or r > 1.5) return .{ 0, 0 };
    var N: f64 = undefined;
    var dN: f64 = undefined;

    if (r >= -1.5 and r <= -0.5) {
        N = 0.5 * r * r + 1.5 * r + 1.125;
        dN = r + 1.5;
    } else if (r <= 0) {
        N = 1 + r;
        dN = 1;
    } else if (r <= 0.5) {
        N = 1 - r;
        dN = -1;
    } else {
        N = 0.5 * r * r - 1.5 * r + 1.125;
        dN = r - 1.5;
    }

    return .{ N, dN };
}

pub fn quad_bespline_type2(r: f64) Vec2 {
    if (r < -1 or r > 1.5) return .{ 0, 0 };
    var N: f64 = undefined;
    var dN: f64 = undefined;

    if (r >= -1 and r <= -0.5) {
        N = 1 + r;
        dN = 1;
    } else if (r <= 0.5) {
        N = -r * r + 0.75;
        dN = -2 * r;
    } else {
        N = 0.5 * r * r - 1.5 * r + 1.125;
        dN = r - 1.5;
    }

    return .{ N, dN };
}

pub fn quad_bespline_type3(r: f64) Vec2 {
    if (r < -1.5 or r > 1.5) return .{ 0, 0 };
    var N: f64 = undefined;
    var dN: f64 = undefined;

    if (r >= -1.5 and r <= -0.5) {
        N = 0.5 * r * r + 1.5 * r + 1.125;
        dN = r + 1.5;
    } else if (r <= 0.5) {
        N = -r * r + 0.75;
        dN = -2 * r;
    } else {
        N = 0.5 * r * r - 1.5 * r + 1.125;
        dN = r - 1.5;
    }

    return .{ N, dN };
}

pub fn quad_bespline_type4(r: f64) Vec2 {
    if (r < -1.5 or r > 1) return .{ 0, 0 };
    var N: f64 = undefined;
    var dN: f64 = undefined;

    if (r >= -1.5 and r <= -0.5) {
        N = 0.5 * r * r + 1.5 * r + 1.125;
        dN = r + 1.5;
    } else if (r <= 0.5) {
        N = -r * r + 0.75;
        dN = -2 * r;
    } else {
        N = 1 - r;
        dN = -1;
    }

    return .{ N, dN };
}

pub fn cubic_bespline_type1(r: f64) Vec2 {
    if (r < -2 or r > 2) return .{ 0, 0 };
    var N: f64 = 0;
    var dN: f64 = 0;
    var r2: f64 = 0;

    if (r >= -2 and r <= -1) {
        r2 = r * r;
        N = 0.16666666666666666 * r * r2 + r2 + 2 * r + 1.3333333333333333;
        dN = 0.5 * r2 + 2 * r + 2.0;
    } else if (r >= -1 and r <= 0) {
        r2 = r * r;
        N = -0.16666666666666666 * r * r2 + r + 1;
        dN = -0.5 * r2 + 1.0;
    } else if (r >= 0 and r <= 1) {
        r2 = r * r;
        N = 0.16666666666666666 * r * r2 - r + 1;
        dN = 0.5 * r2 - 1.0;
    } else {
        r2 = r * r;
        N = -0.16666666666666666 * r * r2 + r2 - 2 * r + 1.3333333333333333;
        dN = -0.5 * r2 + 2 * r - 2.0;
    }

    return .{ N, dN };
}

pub fn cubic_bespline_type2(r: f64) Vec2 {
    if (r < -1 or r > 2) return .{ 0, 0 };
    var N: f64 = 0;
    var dN: f64 = 0;
    var r2: f64 = 0;

    if (r >= -1 and r <= 0) {
        r2 = r * r;
        N = -0.3333333333333333 * r * r2 - r2 + 0.6666666666666666;
        dN = -r2 - 2 * r;
    } else if (r >= 0 and r <= 1) {
        r2 = r * r;
        N = 0.5 * r * r2 - r2 + 0.6666666666666666;
        dN = 1.5 * r2 - 2 * r;
    } else {
        r2 = r * r;
        N = -0.16666666666666666 * r * r2 + r2 - 2 * r + 1.3333333333333333;
        dN = -0.5 * r2 + 2 * r - 2.0;
    }

    return .{ N, dN };
}

pub fn cubic_bespline_type3(r: f64) Vec2 {
    if (r < -2 or r > 2) return .{ 0, 0 };
    var N: f64 = 0;
    var dN: f64 = 0;
    var r2: f64 = 0;

    if (r >= -2 and r <= -1) {
        r2 = r * r;
        N = 0.16666666666666666 * r * r2 + r2 + 2 * r + 1.3333333333333333;
        dN = 0.5 * r2 + 2 * r + 2.0;
    } else if (r >= -1 and r <= 0) {
        r2 = r * r;
        N = -0.5 * r * r2 - r2 + 0.6666666666666666;
        dN = -1.5 * r2 - 2 * r;
    } else if (r >= 0 and r <= 1) {
        r2 = r * r;
        N = 0.5 * r * r2 - r2 + 0.6666666666666666;
        dN = 1.5 * r2 - 2 * r;
    } else {
        r2 = r * r;
        N = -0.16666666666666666 * r * r2 + r2 - 2 * r + 1.3333333333333333;
        dN = -0.5 * r2 + 2 * r - 2.0;
    }

    return .{ N, dN };
}

pub fn cubic_bespline_type4(r: f64) Vec2 {
    if (r < -2 or r > 1) return .{ 0, 0 };
    var N: f64 = 0;
    var dN: f64 = 0;
    var r2: f64 = 0;

    if (r >= -2 and r <= -1) {
        r2 = r * r;
        N = 0.16666666666666666 * r * r2 + r2 + 2 * r + 1.3333333333333333;
        dN = 0.5 * r2 + 2 * r + 2;
    } else if (r >= -1 and r <= 0) {
        r2 = r * r;
        N = -0.5 * r * r2 - r2 + 0.6666666666666666;
        dN = -1.5 * r2 - 2 * r;
    } else if (r >= 0 and r <= 1) {
        r2 = r * r;
        N = 0.3333333333333333 * r * r2 - r2 + 0.6666666667;
        dN = r2 - 2 * r;
    }

    return .{ N, dN };
}
