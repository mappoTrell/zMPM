const std = @import("std");
const testing = std.testing;

export fn add(a: i32, b: i32) i32 {
    return a + b;
}
const y = add();
fn test1(z: u32) void {
    while (z < 10) {
        z -= 1;
    }
    return z;
}
const x = std.Thread.Pool{};
