//! Example usage of ripser-rs
//!
//! Run with: cargo run --example basic

use ripser::{ripser, Barcode};

fn main() {
    println!("=== Example 1: Triangle ===\n");
    example_triangle();

    println!("\n=== Example 2: Point Cloud (Circle) ===\n");
    example_point_cloud();

    println!("\n=== Example 3: With Threshold ===\n");
    example_threshold();
}

// ============================================================================
// Helper functions (for building distance matrices)
// ============================================================================

/// Compute Euclidean distance matrix from a point cloud.
///
/// # Arguments
/// * `points` - Flat array: [x0, y0, z0, x1, y1, z1, ...]
/// * `n` - Number of points
/// * `dim` - Dimension of each point
///
/// # Returns
/// Lower triangular distance matrix
fn euclidean_distance_matrix(points: &[f32], n: usize, dim: usize) -> Vec<f32> {
    assert_eq!(points.len(), n * dim);

    let mut distances = Vec::with_capacity(n * (n - 1) / 2);
    for i in 1..n {
        for j in 0..i {
            let mut d2 = 0.0f32;
            for k in 0..dim {
                let diff = points[i * dim + k] - points[j * dim + k];
                d2 += diff * diff;
            }
            distances.push(d2.sqrt());
        }
    }
    distances
}

/// Filter barcodes by homological dimension.
fn filter_by_dim(barcodes: &[Barcode], dim: i32) -> Vec<&Barcode> {
    barcodes.iter().filter(|b| b.dim == dim).collect()
}

/// Get only infinite (essential) bars.
fn infinite_bars(barcodes: &[Barcode]) -> Vec<&Barcode> {
    barcodes.iter().filter(|b| b.is_infinite()).collect()
}

// ============================================================================
// Examples
// ============================================================================

fn example_triangle() {
    // Three points forming an equilateral triangle with side length 1
    let distances = vec![
        1.0, // d(1,0)
        1.0, 1.0, // d(2,0), d(2,1)
    ];

    let barcodes = ripser(&distances, 3, 1, None);

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
    println!("  H0: {} bars (connected components)", h0.len());
    println!(
        "    - {} merge events",
        h0.len()
            - infinite_bars(&barcodes)
                .iter()
                .filter(|b| b.dim == 0)
                .count()
    );
    println!(
        "    - {} final component(s)",
        infinite_bars(&barcodes)
            .iter()
            .filter(|b| b.dim == 0)
            .count()
    );
    println!("  H1: {} bars (loops)", h1.len());
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

    // Compute distance matrix from point cloud
    let distances = euclidean_distance_matrix(&points, n, 2);

    // Compute persistent homology up to dimension 1
    let barcodes = ripser(&distances, n, 1, None);

    println!("Circle with {} points:", n);

    let h0 = filter_by_dim(&barcodes, 0);
    let h1 = filter_by_dim(&barcodes, 1);

    println!(
        "  H0: {} bars ({} infinite)",
        h0.len(),
        infinite_bars(&barcodes)
            .iter()
            .filter(|b| b.dim == 0)
            .count()
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
    let sqrt2 = 2.0_f32.sqrt();

    // Square with side 1, diagonal sqrt(2)
    // Points: 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
    let distances = vec![
        1.0, // d(1,0) = 1
        1.0, sqrt2, // d(2,0) = 1, d(2,1) = sqrt(2)
        sqrt2, 1.0, 1.0, // d(3,0) = sqrt(2), d(3,1) = 1, d(3,2) = 1
    ];

    println!("Square (side=1, diagonal=√2≈{:.3}):", sqrt2);

    // Without threshold - all edges included
    let barcodes_full = ripser(&distances, 4, 1, None);
    let h1_full = filter_by_dim(&barcodes_full, 1);
    println!("\n  Without threshold: {} H1 bar(s)", h1_full.len());
    for bar in &h1_full {
        println!(
            "    [{:.3}, {:.3}), persistence = {:.3}",
            bar.birth,
            bar.death,
            bar.persistence()
        );
    }

    // With threshold = 1.2 (excludes diagonal edges)
    let barcodes_thresh = ripser(&distances, 4, 1, Some(1.2));
    let h1_thresh = filter_by_dim(&barcodes_thresh, 1);
    println!("\n  With threshold 1.2: {} H1 bar(s)", h1_thresh.len());
    for bar in &h1_thresh {
        let death_str = if bar.is_infinite() {
            "∞".to_string()
        } else {
            format!("{:.3}", bar.death)
        };
        println!(
            "    [{:.3}, {}), persistence = {:.3}",
            bar.birth,
            death_str,
            bar.persistence()
        );
    }
    println!("\n  (With threshold, diagonal edges never added, so the loop persists forever)");
}
