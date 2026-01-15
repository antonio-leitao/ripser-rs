//! Basic example showing how to use ripser-rs
//!
//! Run with: cargo run --example basic

use ripser::{euclidean_distance_matrix, filter_by_dim, infinite_bars, ripser};

fn main() {
    println!("=== Example 1: Triangle ===\n");
    example_triangle();

    println!("\n=== Example 2: Point Cloud ===\n");
    example_point_cloud();

    println!("\n=== Example 3: With Threshold ===\n");
    example_threshold();
}

fn example_triangle() {
    // Three points forming an equilateral triangle with side length 1
    // This is a simple shape that demonstrates H0 (connected components) and H1 (loops)
    let distances = vec![
        1.0, // d(1,0)
        1.0, 1.0, // d(2,0), d(2,1)
    ];

    let barcodes = ripser(&distances, 3, 1, f32::INFINITY);

    println!("Triangle with side length 1:");
    println!("Number of barcodes: {}", barcodes.len());

    for bar in &barcodes {
        let death_str = if bar.is_infinite() {
            "∞".to_string()
        } else {
            format!("{:.3}", bar.death)
        };
        println!("  H{}: [{:.3}, {})", bar.dim, bar.birth, death_str);
    }

    let h0 = filter_by_dim(&barcodes, 0);
    let h1 = filter_by_dim(&barcodes, 1);

    println!("\nInterpretation:");
    println!("  H0 bars: {} (connected components)", h0.len());
    println!("    - {} merge events", h0.len() - infinite_bars(&h0).len());
    println!("    - {} final component(s)", infinite_bars(&h0).len());
    println!("  H1 bars: {} (loops)", h1.len());
}

fn example_point_cloud() {
    // A circle of 8 points
    let n = 8;
    let mut points = Vec::with_capacity(n * 2);

    for i in 0..n {
        let angle = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
        points.push(angle.cos());
        points.push(angle.sin());
    }

    // Compute distance matrix
    let distances = euclidean_distance_matrix(&points, n, 2);

    // Compute persistent homology up to dimension 1
    let barcodes = ripser(&distances, n, 1, f32::INFINITY);

    println!("Circle with {} points:", n);

    let h0 = filter_by_dim(&barcodes, 0);
    let h1 = filter_by_dim(&barcodes, 1);

    println!(
        "  H0: {} bars ({} infinite)",
        h0.len(),
        infinite_bars(&h0).len()
    );
    println!("  H1: {} bars", h1.len());

    // Find the most persistent H1 bar (the main loop)
    if let Some(main_loop) = h1
        .iter()
        .max_by(|a, b| a.persistence().partial_cmp(&b.persistence()).unwrap())
    {
        println!(
            "\n  Most persistent H1 bar: [{:.3}, {:.3})",
            main_loop.birth, main_loop.death
        );
        println!("  Persistence: {:.3}", main_loop.persistence());
    }
}

fn example_threshold() {
    // Four points in a square pattern
    // Using threshold to filter out long edges
    let sqrt2 = 2.0_f32.sqrt();

    // Square with side 1
    let distances = vec![
        1.0, // d(1,0)
        1.0, sqrt2, // d(2,0), d(2,1)
        sqrt2, 1.0, 1.0, // d(3,0), d(3,1), d(3,2)
    ];

    println!("Square (side=1, diagonal=√2):");

    // Without threshold
    let barcodes_full = ripser(&distances, 4, 1, f32::INFINITY);
    let h1_full = filter_by_dim(&barcodes_full, 1);
    println!("  Without threshold: {} H1 bars", h1_full.len());
    for bar in &h1_full {
        println!(
            "    [{:.3}, {:.3}), persistence = {:.3}",
            bar.birth,
            bar.death,
            bar.persistence()
        );
    }

    // With threshold = 1.2 (excludes diagonals)
    let barcodes_thresh = ripser(&distances, 4, 1, 1.2);
    let h1_thresh = filter_by_dim(&barcodes_thresh, 1);
    println!("  With threshold 1.2: {} H1 bars", h1_thresh.len());
    for bar in &h1_thresh {
        println!(
            "    [{:.3}, {:.3}), persistence = {:.3}",
            bar.birth,
            bar.death,
            bar.persistence()
        );
    }
}
