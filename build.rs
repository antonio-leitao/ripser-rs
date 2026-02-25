use std::env;

fn main() {
    let target = env::var("TARGET").unwrap();

    let mut build = cc::Build::new();

    build
        .cpp(true)
        .file("ripser.cpp")
        .define("RIPSER_FFI", None)
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-O3")
        // Suppress warnings from ripser.cpp (it's third-party code)
        .flag_if_supported("-w");

    // Platform-specific flags
    if target.contains("msvc") {
        build.flag("/EHsc");
    }

    build.compile("ripser");

    // Build ripser128 (128-bit index version) — not available on Windows
    // because MSVC does not support __int128.
    if !target.contains("windows") {
        let mut build128 = cc::Build::new();

        build128
            .cpp(true)
            .file("ripser128.cpp")
            .define("RIPSER_FFI", None)
            // Rename symbols to avoid duplicate symbol errors when linking with ripser.cpp.
            .define("is_prime", Some("is_prime128"))
            .define(
                "multiplicative_inverse_vector",
                Some("multiplicative_inverse_vector128"),
            )
            .define("check_overflow", Some("check_overflow128"))
            .define("make_entry", Some("make_entry128"))
            .define("get_index", Some("get_index128"))
            .define("get_entry", Some("get_entry128"))
            .define("get_diameter", Some("get_diameter128"))
            .define("get_coefficient", Some("get_coefficient128"))
            .define("set_coefficient", Some("set_coefficient128"))
            .define(
                "greater_diameter_or_smaller_index",
                Some("greater_diameter_or_smaller_index128"),
            )
            .define(
                "greater_diameter_or_smaller_index_comp",
                Some("greater_diameter_or_smaller_index_comp128"),
            )
            .define("diameter_entry_t", Some("diameter_entry_128_t"))
            .define("diameter_index_t", Some("diameter_index_128_t"))
            .define("index_diameter_t", Some("index_diameter_128_t"))
            .define("entry_t", Some("entry_128_t"))
            .define(
                "compressed_distance_matrix",
                Some("compressed_distance_matrix128"),
            )
            .define(
                "compressed_lower_distance_matrix",
                Some("compressed_lower_distance_matrix128"),
            )
            .define(
                "compressed_upper_distance_matrix",
                Some("compressed_upper_distance_matrix128"),
            )
            .define(
                "compressed_sparse_matrix",
                Some("compressed_sparse_matrix128"),
            )
            .define("sparse_distance_matrix", Some("sparse_distance_matrix128"))
            .define(
                "borrowed_lower_distance_matrix",
                Some("borrowed_lower_distance_matrix128"),
            )
            .define("barcode_collector", Some("barcode_collector128"))
            .define("union_find", Some("union_find128"))
            .define("binomial_coeff_table", Some("binomial_coeff_table128"))
            .define("ripser", Some("ripser128_impl"))
            .define("get_max", Some("get_max128"))
            .flag_if_supported("-std=c++11")
            .flag_if_supported("-O3")
            .flag_if_supported("-w");

        build128.compile("ripser128");
    }

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=ripser.cpp");
    println!("cargo:rerun-if-changed=ripser128.cpp");
}
