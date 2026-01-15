//! # ripser-rs
//!
//! Rust bindings for [Ripser](https://github.com/Ripser/ripser), a lean C++ code
//! for computation of Vietoris-Rips persistence barcodes.
//!
//! ## Example
//!
//! ```rust
//! use ripser::{ripser, Barcode};
//!
//! // Lower triangular distance matrix for 4 points
//! // Layout: d(1,0), d(2,0), d(2,1), d(3,0), d(3,1), d(3,2)
//! let distances = vec![
//!     1.0,        // d(1,0)
//!     2.0, 1.5,   // d(2,0), d(2,1)
//!     3.0, 2.5, 1.0, // d(3,0), d(3,1), d(3,2)
//! ];
//!
//! let barcodes = ripser(&distances, 4, 1, f32::INFINITY);
//!
//! for bar in &barcodes {
//!     println!("dim {}: [{}, {})", bar.dim, bar.birth, bar.death);
//! }
//! ```

use std::f32;

/// A persistence barcode (birth-death pair) with dimension.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Barcode {
    /// Homological dimension (0 = connected components, 1 = loops, 2 = voids, etc.)
    pub dim: i32,
    /// Birth time (filtration value where the feature appears)
    pub birth: f32,
    /// Death time (filtration value where the feature disappears, f32::INFINITY for essential features)
    pub death: f32,
}

impl Barcode {
    /// Returns true if this is an infinite bar (essential feature).
    #[inline]
    pub fn is_infinite(&self) -> bool {
        self.death.is_infinite()
    }

    /// Returns the persistence (death - birth). Returns f32::INFINITY for infinite bars.
    #[inline]
    pub fn persistence(&self) -> f32 {
        self.death - self.birth
    }
}

/// FFI result struct - must match C++ definition exactly.
#[repr(C)]
struct BarcodeResult {
    data: *mut Barcode,
    len: usize,
}

extern "C" {
    fn ripser_from_lower_distance_matrix(
        distances: *const f32,
        n: usize,
        max_dim: i32,
        threshold: f32,
    ) -> BarcodeResult;

    fn free_barcodes(result: BarcodeResult);
}

/// Compute persistent homology of a point cloud given as a distance matrix.
///
/// # Arguments
///
/// * `distances` - Lower triangular distance matrix in row-major order.
///   Length must be exactly `n * (n - 1) / 2`.
///   Layout: `[d(1,0), d(2,0), d(2,1), d(3,0), d(3,1), d(3,2), ...]`
/// * `n` - Number of points.
/// * `max_dim` - Maximum homological dimension to compute (0, 1, 2, ...).
///   Dimension 0 gives connected components, dimension 1 gives loops, etc.
/// * `threshold` - Maximum filtration value. Use `f32::INFINITY` for no threshold.
///
/// # Returns
///
/// A vector of `Barcode` structs representing the persistence diagram.
///
/// # Panics
///
/// Panics if `distances.len() != n * (n - 1) / 2`.
///
/// # Example
///
/// ```rust
/// use ripser::ripser;
///
/// // Triangle with vertices at distance 1 from each other
/// let distances = vec![1.0, 1.0, 1.0]; // d(1,0), d(2,0), d(2,1)
/// let barcodes = ripser(&distances, 3, 1, f32::INFINITY);
///
/// // Should have:
/// // - 3 H0 bars (one infinite, two finite that merge)
/// // - 1 H1 bar (the triangle forms a loop)
/// for bar in &barcodes {
///     println!("H{}: [{}, {})", bar.dim, bar.birth, bar.death);
/// }
/// ```
pub fn ripser(distances: &[f32], n: usize, max_dim: i32, threshold: f32) -> Vec<Barcode> {
    let expected_len = n * (n.saturating_sub(1)) / 2;
    assert_eq!(
        distances.len(),
        expected_len,
        "Distance matrix has wrong size: expected {} for n={}, got {}",
        expected_len,
        n,
        distances.len()
    );

    // Handle edge cases
    if n == 0 {
        return Vec::new();
    }

    // Call the C++ implementation
    let result = unsafe {
        ripser_from_lower_distance_matrix(
            distances.as_ptr(),
            n,
            max_dim,
            if threshold.is_infinite() {
                f32::MAX
            } else {
                threshold
            },
        )
    };

    // Convert to Vec and free the C++ memory
    let barcodes = if result.data.is_null() || result.len == 0 {
        Vec::new()
    } else {
        // Safety: result.data points to a valid array of result.len Barcodes
        // allocated with malloc by the C++ code
        let slice = unsafe { std::slice::from_raw_parts(result.data, result.len) };
        slice.to_vec()
    };

    // Free the C++ allocated memory
    unsafe { free_barcodes(result) };

    barcodes
}

/// Convenience function to filter barcodes by dimension.
pub fn filter_by_dim(barcodes: &[Barcode], dim: i32) -> Vec<Barcode> {
    barcodes.iter().filter(|b| b.dim == dim).copied().collect()
}

/// Convenience function to get only infinite bars.
pub fn infinite_bars(barcodes: &[Barcode]) -> Vec<Barcode> {
    barcodes
        .iter()
        .filter(|b| b.is_infinite())
        .copied()
        .collect()
}

/// Convenience function to get only finite bars.
pub fn finite_bars(barcodes: &[Barcode]) -> Vec<Barcode> {
    barcodes
        .iter()
        .filter(|b| !b.is_infinite())
        .copied()
        .collect()
}

/// Helper to create a lower triangular distance matrix from a full matrix.
///
/// # Arguments
///
/// * `full_matrix` - n×n distance matrix as a flat array in row-major order.
/// * `n` - Matrix dimension.
///
/// # Returns
///
/// Lower triangular part: `[d(1,0), d(2,0), d(2,1), ...]`
pub fn full_to_lower_triangular(full_matrix: &[f32], n: usize) -> Vec<f32> {
    assert_eq!(full_matrix.len(), n * n, "Expected {}x{} matrix", n, n);

    let mut lower = Vec::with_capacity(n * (n - 1) / 2);
    for i in 1..n {
        for j in 0..i {
            lower.push(full_matrix[i * n + j]);
        }
    }
    lower
}

/// Compute a distance matrix from a point cloud using Euclidean distance.
///
/// # Arguments
///
/// * `points` - Flat array of coordinates: `[x0, y0, z0, x1, y1, z1, ...]`
/// * `n` - Number of points.
/// * `dim` - Dimension of each point.
///
/// # Returns
///
/// Lower triangular distance matrix.
pub fn euclidean_distance_matrix(points: &[f32], n: usize, dim: usize) -> Vec<f32> {
    assert_eq!(
        points.len(),
        n * dim,
        "Expected {} points of dimension {}",
        n,
        dim
    );

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_point() {
        let distances: Vec<f32> = vec![];
        let barcodes = ripser(&distances, 1, 1, f32::INFINITY);

        // Single point should have exactly one infinite H0 bar
        assert_eq!(barcodes.len(), 1);
        assert_eq!(barcodes[0].dim, 0);
        assert_eq!(barcodes[0].birth, 0.0);
        assert!(barcodes[0].is_infinite());
    }

    #[test]
    fn test_two_points() {
        let distances = vec![1.0]; // d(1,0) = 1
        let barcodes = ripser(&distances, 2, 1, f32::INFINITY);

        // Two points at distance 1:
        // - Two H0 components that merge at time 1
        // - One survives to infinity
        let h0: Vec<_> = barcodes.iter().filter(|b| b.dim == 0).collect();

        // Should have exactly 2 H0 bars
        assert_eq!(h0.len(), 2);

        // One should be [0, 1) and one should be [0, inf)
        let finite: Vec<_> = h0.iter().filter(|b| !b.is_infinite()).collect();
        let infinite: Vec<_> = h0.iter().filter(|b| b.is_infinite()).collect();

        assert_eq!(finite.len(), 1);
        assert_eq!(infinite.len(), 1);
        assert_eq!(finite[0].death, 1.0);
    }

    #[test]
    fn test_equilateral_triangle() {
        // Three points forming an equilateral triangle with side length 1
        let distances = vec![1.0, 1.0, 1.0];
        let barcodes = ripser(&distances, 3, 1, f32::INFINITY);

        let h0: Vec<_> = filter_by_dim(&barcodes, 0);
        let h1: Vec<_> = filter_by_dim(&barcodes, 1);

        // H0: 3 components, 2 merge, 1 survives
        assert_eq!(h0.len(), 3);
        assert_eq!(infinite_bars(&h0).len(), 1);

        // H1: For an equilateral triangle, a loop forms at time 1 but is immediately killed
        // by the 2-simplex (since all edges have length 1). So persistence = 0 and it's filtered out.
        // This is correct behavior - no persistent H1 feature.
        assert_eq!(h1.len(), 0);
    }

    #[test]
    fn test_square() {
        // Four points forming a square
        // 0 -- 1
        // |    |
        // 2 -- 3
        // Side length 1, diagonal sqrt(2)
        let sqrt2 = 2.0_f32.sqrt();
        let distances = vec![
            1.0, // d(1,0)
            1.0, sqrt2, // d(2,0), d(2,1)
            sqrt2, 1.0, 1.0, // d(3,0), d(3,1), d(3,2)
        ];
        let barcodes = ripser(&distances, 4, 1, f32::INFINITY);

        let h0: Vec<_> = filter_by_dim(&barcodes, 0);
        let h1: Vec<_> = filter_by_dim(&barcodes, 1);

        // Should have 4 H0 bars (one infinite)
        assert_eq!(h0.len(), 4);
        assert_eq!(infinite_bars(&h0).len(), 1);

        // Should have 1 H1 bar (the square loop)
        assert_eq!(h1.len(), 1);
    }

    #[test]
    fn test_full_to_lower() {
        let full = vec![0.0, 1.0, 2.0, 1.0, 0.0, 3.0, 2.0, 3.0, 0.0];
        let lower = full_to_lower_triangular(&full, 3);
        assert_eq!(lower, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_euclidean_distance() {
        // Three points: (0,0), (1,0), (0,1)
        let points = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let distances = euclidean_distance_matrix(&points, 3, 2);

        assert_eq!(distances.len(), 3);
        assert!((distances[0] - 1.0).abs() < 1e-6); // d(1,0) = 1
        assert!((distances[1] - 1.0).abs() < 1e-6); // d(2,0) = 1
        assert!((distances[2] - 2.0_f32.sqrt()).abs() < 1e-6); // d(2,1) = sqrt(2)
    }

    #[test]
    fn test_threshold() {
        // Two points at distance 2, with threshold 1
        let distances = vec![2.0];
        let barcodes = ripser(&distances, 2, 1, 1.0);

        // With threshold 1, the edge at distance 2 is never added
        // So we should have 2 independent H0 components
        let h0: Vec<_> = filter_by_dim(&barcodes, 0);
        let infinite = infinite_bars(&h0);

        // Both points remain as separate infinite components
        assert_eq!(infinite.len(), 2);
    }

    #[test]
    #[should_panic(expected = "Distance matrix has wrong size")]
    fn test_wrong_size() {
        let distances = vec![1.0, 2.0]; // Wrong size for n=3
        let _ = ripser(&distances, 3, 1, f32::INFINITY);
    }
}
