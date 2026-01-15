use std::env;

fn main() {
    let mut build = cc::Build::new();
    
    build
        .cpp(true)
        .file("ripser.cpp")
        .define("RIPSER_FFI", None)
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-O3")
        .flag_if_supported("-march=native")
        // Suppress warnings from ripser.cpp (it's third-party code)
        .flag_if_supported("-w");
    
    // Platform-specific flags
    let target = env::var("TARGET").unwrap();
    if target.contains("msvc") {
        build.flag("/EHsc");
    }
    
    build.compile("ripser");
    
    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=ripser.cpp");
}
