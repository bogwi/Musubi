const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Module
    const mod = b.createModule(.{
        .root_source_file = b.path("src/musubi.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Library
    const lib = b.addLibrary(.{
        .name = "Musubi",
        .root_module = mod,
    });

    b.installArtifact(lib);

    // Tests
    const test_step = b.step("test", "Run library tests");

    const main_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_main_tests = b.addRunArtifact(main_tests);
    test_step.dependOn(&run_main_tests.step);
}
