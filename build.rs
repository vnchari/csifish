use std::env;
use std::path::Path;

fn main() {
    // println!("cargo:rerun-if-changed=src/inv");
    let target = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let rustflags = "-C target-cpu=native";
    let asm_path = match target.as_str() {
        "x86" | "x86_64" => Path::new("src/csifish/field_arithmetic/inv/modinv_x86.S"),
        "arm" | "aarch64" => Path::new("src/csifish/field_arithmetic/inv/modinv_arm.S"),
        _ => panic!("unsupported architecture"),
    };
    let header = Path::new("src/csifish/field_arithmetic/inv/include/");
    cc::Build::new()
        .static_flag(true)
        .include(header)
        .file(asm_path)
        .opt_level(3)
        .compile("modinv");

    println!("cargo:rustc-link-search=src/inv");
    println!("cargo:rustc-link-lib=static=modinv");

    println!("cargo:rustc-env=RUSTFLAGS={}", rustflags);
}
