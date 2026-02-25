//! # ripser-rs
//!
//! Minimal Rust bindings for [Ripser](https://github.com/Ripser/ripser), a lean C++ code
//! for computation of Vietoris-Rips persistence barcodes.
//!
//! ## Example
//!
//! ```rust
//! use ripser::{ripser, Barcode};
//!
//! // Lower triangular distance matrix for 3 points (equilateral triangle)
//! // Layout: d(1,0), d(2,0), d(2,1)
//! let distances = vec![1.0, 1.0, 1.0];
//!
//! let barcodes = ripser(&distances, 3, 1, None);
//!
//! for bar in &barcodes {
//!     println!("H{}: [{}, {})", bar.dim, bar.birth, bar.death);
//! }
//! ```
//!
//! ## Distance Matrix Format
//!
//! Lower triangular, row-major order:
//!
//! ```text
//!      0   1   2   3
//! 0    -
//! 1  d10   -
//! 2  d20 d21   -
//! 3  d30 d31 d32   -
//! ```
//!
//! Flattened: `[d10, d20, d21, d30, d31, d32, ...]`
//!
//! Length must be exactly `n * (n - 1) / 2`.
//!
//! ## Platform Support
//!
//! The standard `ripser()` function is available on all platforms.
//!
//! The `ripser128()` function (128-bit simplex indices for large point clouds) is only
//! available on Unix platforms (Linux, macOS) because MSVC does not support `__int128`.

/// A persistence barcode (birth-death pair) with dimension.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Barcode {
    /// Homological dimension (0 = connected components, 1 = loops, 2 = voids, etc.)
    pub dim: i32,
    /// Birth time (filtration value where the feature appears)
    pub birth: f32,
    /// Death time (filtration value where the feature disappears)
    /// `f32::INFINITY` indicates an essential feature that never dies.
    pub death: f32,
}

impl Barcode {
    /// Returns true if this is an infinite bar (essential feature).
    #[inline]
    pub fn is_infinite(&self) -> bool {
        self.death.is_infinite()
    }

    /// Returns the persistence (death - birth).
    /// Returns `f32::INFINITY` for infinite bars.
    #[inline]
    pub fn persistence(&self) -> f32 {
        self.death - self.birth
    }
}

// FFI types - must match C++ definitions exactly
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

#[cfg(not(target_os = "windows"))]
extern "C" {
    fn ripser128_from_lower_distance_matrix(
        distances: *const f32,
        n: usize,
        max_dim: i32,
        threshold: f32,
    ) -> BarcodeResult;

    fn free_barcodes128(result: BarcodeResult);
}

/// RAII guard for panic safety - ensures FFI memory is always freed.
struct BarcodeGuard(BarcodeResult);

impl Drop for BarcodeGuard {
    fn drop(&mut self) {
        unsafe {
            free_barcodes(BarcodeResult {
                data: self.0.data,
                len: self.0.len,
            })
        }
    }
}

/// Compute persistent homology of a point cloud given as a distance matrix.
///
/// # Arguments
///
/// * `distances` - Lower triangular distance matrix in row-major order.
///   Length must be exactly `n * (n - 1) / 2`.
///   Layout: `[d(1,0), d(2,0), d(2,1), d(3,0), d(3,1), d(3,2), ...]`
///   Values should be finite and non-negative.
///
/// * `n` - Number of points.
///
/// * `max_dim` - Maximum homological dimension to compute (0, 1, 2, ...).
///   - Dimension 0: connected components
///   - Dimension 1: loops/holes
///   - Dimension 2: voids
///   - etc.
///
/// * `threshold` - Maximum filtration value. Edges with distance > threshold
///   are ignored. Use `None` to auto-compute the enclosing radius (recommended).
///
/// # Returns
///
/// A vector of `Barcode` structs representing the persistence diagram.
/// Barcodes are ordered by dimension, then by appearance in the filtration.
///
/// # Panics
///
/// - If `distances.len() != n * (n - 1) / 2`
/// - If `max_dim < 0`
///
/// # Example
///
/// ```rust
/// use ripser::ripser;
///
/// // Two points at distance 1
/// let distances = vec![1.0];
/// let barcodes = ripser(&distances, 2, 0, None);
///
/// // H0: two components merge at time 1
/// assert_eq!(barcodes.len(), 2);
/// ```
pub fn ripser(distances: &[f32], n: usize, max_dim: i32, threshold: Option<f32>) -> Vec<Barcode> {
    // Validate inputs
    let expected_len = n * (n.saturating_sub(1)) / 2;
    assert_eq!(
        distances.len(),
        expected_len,
        "Distance matrix has wrong size: expected {} for n={}, got {}",
        expected_len,
        n,
        distances.len()
    );
    assert!(
        max_dim >= 0,
        "max_dim must be non-negative, got {}",
        max_dim
    );

    // Handle degenerate cases
    if n == 0 {
        return Vec::new();
    }

    // Convert threshold: None or infinite → f32::MAX (triggers auto-compute in C++)
    let thresh = match threshold {
        None => f32::MAX,
        Some(t) if t.is_infinite() => f32::MAX,
        Some(t) => t,
    };

    // Call FFI and wrap result in guard for panic safety
    let guard = BarcodeGuard(unsafe {
        ripser_from_lower_distance_matrix(distances.as_ptr(), n, max_dim, thresh)
    });

    // Convert to Vec (guard ensures cleanup even if this panics)
    if guard.0.data.is_null() || guard.0.len == 0 {
        Vec::new()
    } else {
        let slice = unsafe { std::slice::from_raw_parts(guard.0.data, guard.0.len) };
        slice.to_vec()
    }
    // guard drops here, always freeing the C++ memory
}

/// RAII guard for panic safety - ensures FFI memory is always freed (128-bit version).
#[cfg(not(target_os = "windows"))]
struct BarcodeGuard128(BarcodeResult);

#[cfg(not(target_os = "windows"))]
impl Drop for BarcodeGuard128 {
    fn drop(&mut self) {
        unsafe {
            free_barcodes128(BarcodeResult {
                data: self.0.data,
                len: self.0.len,
            })
        }
    }
}

/// Compute persistent homology using the 128-bit index variant of Ripser.
///
/// This is identical to [`ripser()`] but uses `__int128_t` internally for simplex
/// indices, allowing it to handle larger simplicial complexes that would overflow
/// 64-bit indices.
///
/// See [`ripser()`] for full documentation of parameters and return values.
///
/// # Platform Availability
///
/// This function is only available on Unix platforms (Linux, macOS).
/// MSVC does not support `__int128`.
#[cfg(not(target_os = "windows"))]
pub fn ripser128(
    distances: &[f32],
    n: usize,
    max_dim: i32,
    threshold: Option<f32>,
) -> Vec<Barcode> {
    // Validate inputs
    let expected_len = n * (n.saturating_sub(1)) / 2;
    assert_eq!(
        distances.len(),
        expected_len,
        "Distance matrix has wrong size: expected {} for n={}, got {}",
        expected_len,
        n,
        distances.len()
    );
    assert!(
        max_dim >= 0,
        "max_dim must be non-negative, got {}",
        max_dim
    );

    // Handle degenerate cases
    if n == 0 {
        return Vec::new();
    }

    // Convert threshold: None or infinite → f32::MAX (triggers auto-compute in C++)
    let thresh = match threshold {
        None => f32::MAX,
        Some(t) if t.is_infinite() => f32::MAX,
        Some(t) => t,
    };

    // Call FFI and wrap result in guard for panic safety
    let guard = BarcodeGuard128(unsafe {
        ripser128_from_lower_distance_matrix(distances.as_ptr(), n, max_dim, thresh)
    });

    // Convert to Vec (guard ensures cleanup even if this panics)
    if guard.0.data.is_null() || guard.0.len == 0 {
        Vec::new()
    } else {
        let slice = unsafe { std::slice::from_raw_parts(guard.0.data, guard.0.len) };
        slice.to_vec()
    }
    // guard drops here, always freeing the C++ memory
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_point() {
        let distances: Vec<f32> = vec![];
        let barcodes = ripser(&distances, 1, 1, None);

        assert_eq!(barcodes.len(), 1);
        assert_eq!(barcodes[0].dim, 0);
        assert_eq!(barcodes[0].birth, 0.0);
        assert!(barcodes[0].is_infinite());
    }

    #[test]
    fn test_two_points() {
        let distances = vec![1.0];
        let barcodes = ripser(&distances, 2, 1, None);

        let h0: Vec<_> = barcodes.iter().filter(|b| b.dim == 0).collect();
        assert_eq!(h0.len(), 2);

        let finite: Vec<_> = h0.iter().filter(|b| !b.is_infinite()).collect();
        let infinite: Vec<_> = h0.iter().filter(|b| b.is_infinite()).collect();

        assert_eq!(finite.len(), 1);
        assert_eq!(infinite.len(), 1);
        assert_eq!(finite[0].death, 1.0);
    }

    #[test]
    fn test_triangle() {
        let distances = vec![1.0, 1.0, 1.0];
        let barcodes = ripser(&distances, 3, 1, None);

        let h0: Vec<_> = barcodes.iter().filter(|b| b.dim == 0).collect();
        let h1: Vec<_> = barcodes.iter().filter(|b| b.dim == 1).collect();

        assert_eq!(h0.len(), 3);
        assert_eq!(h0.iter().filter(|b| b.is_infinite()).count(), 1);
        // H1 loop has zero persistence (birth == death), filtered out
        assert_eq!(h1.len(), 0);
    }

    #[test]
    fn test_threshold() {
        let distances = vec![2.0];
        let barcodes = ripser(&distances, 2, 1, Some(1.0));

        let infinite: Vec<_> = barcodes.iter().filter(|b| b.is_infinite()).collect();
        assert_eq!(infinite.len(), 2); // Two separate components
    }

    #[test]
    fn test_empty() {
        let distances: Vec<f32> = vec![];
        let barcodes = ripser(&distances, 0, 1, None);
        assert!(barcodes.is_empty());
    }

    #[test]
    #[should_panic(expected = "Distance matrix has wrong size")]
    fn test_wrong_size() {
        let distances = vec![1.0, 2.0];
        let _ = ripser(&distances, 3, 1, None);
    }

    #[test]
    #[should_panic(expected = "max_dim must be non-negative")]
    fn test_negative_dim() {
        let distances = vec![1.0];
        let _ = ripser(&distances, 2, -1, None);
    }

    // ========================================================================
    // ripser128 tests (not available on Windows)
    // ========================================================================

    #[cfg(not(target_os = "windows"))]
    #[test]
    fn test_128_single_point() {
        let distances: Vec<f32> = vec![];
        let barcodes = ripser128(&distances, 1, 1, None);

        assert_eq!(barcodes.len(), 1);
        assert_eq!(barcodes[0].dim, 0);
        assert_eq!(barcodes[0].birth, 0.0);
        assert!(barcodes[0].is_infinite());
    }

    #[cfg(not(target_os = "windows"))]
    #[test]
    fn test_128_two_points() {
        let distances = vec![1.0];
        let barcodes = ripser128(&distances, 2, 1, None);

        let h0: Vec<_> = barcodes.iter().filter(|b| b.dim == 0).collect();
        assert_eq!(h0.len(), 2);

        let finite: Vec<_> = h0.iter().filter(|b| !b.is_infinite()).collect();
        let infinite: Vec<_> = h0.iter().filter(|b| b.is_infinite()).collect();

        assert_eq!(finite.len(), 1);
        assert_eq!(infinite.len(), 1);
        assert_eq!(finite[0].death, 1.0);
    }

    #[cfg(not(target_os = "windows"))]
    #[test]
    fn test_128_triangle() {
        let distances = vec![1.0, 1.0, 1.0];
        let barcodes = ripser128(&distances, 3, 1, None);

        let h0: Vec<_> = barcodes.iter().filter(|b| b.dim == 0).collect();
        let h1: Vec<_> = barcodes.iter().filter(|b| b.dim == 1).collect();

        assert_eq!(h0.len(), 3);
        assert_eq!(h0.iter().filter(|b| b.is_infinite()).count(), 1);
        assert_eq!(h1.len(), 0);
    }

    #[cfg(not(target_os = "windows"))]
    #[test]
    fn test_128_matches_64bit() {
        // Verify ripser() and ripser128() produce identical results
        let distances = vec![1.0, 1.41, 2.0, 1.0, 1.41, 1.0];
        let n = 4;

        let barcodes_64 = ripser(&distances, n, 2, None);
        let barcodes_128 = ripser128(&distances, n, 2, None);

        assert_eq!(
            barcodes_64.len(),
            barcodes_128.len(),
            "Different number of barcodes: 64-bit={}, 128-bit={}",
            barcodes_64.len(),
            barcodes_128.len()
        );

        for (b64, b128) in barcodes_64.iter().zip(barcodes_128.iter()) {
            assert_eq!(b64.dim, b128.dim);
            assert_eq!(b64.birth, b128.birth);
            assert_eq!(b64.death, b128.death);
        }
    }
}
