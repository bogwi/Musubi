const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});

    // Module
    // _ = b.addModule("musubi", .{ .source_file = .{ .path = "src/musubi.zig" } });

    // Library
    const lib = b.addStaticLibrary(.{
        .name = "Musubi",
        .root_source_file = .{ .path = "src/musubi.zig" },
        .target = target,
        .optimize = .ReleaseSafe,
        .version = .{ .major = 1, .minor = 0, .patch = 0 },
    });

    b.installArtifact(lib);

    // Tests
    const test_step = b.step("test", "Run library tests");

    const main_tests = b.addTest(.{
        .root_source_file = .{ .path = "src/musubi.zig" },
        .target = target,
        .optimize = .ReleaseSafe,
    });

    const run_main_tests = b.addRunArtifact(main_tests);
    test_step.dependOn(&run_main_tests.step);
}
