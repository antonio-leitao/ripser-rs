# ripser-rs

Rust bindings for [Ripser](https://github.com/Ripser/ripser), a lean C++ code for computation of Vietoris-Rips persistence barcodes.

## Features

- **Zero-copy input**: Distance matrix is borrowed, not copied
- **Minimal overhead**: Single allocation for output
- **Identical results**: Same algorithm as original Ripser CLI
- **Safe Rust API**: Automatic memory management

## Usage

```rust
use ripser::{ripser, Barcode, filter_by_dim};
// Lower triangular distance matrix for 4 points
// Layout: d(1,0), d(2,0), d(2,1), d(3,0), d(3,1), d(3,2)
let distances = vec![
    1.0,        // d(1,0)
    2.0, 1.5,   // d(2,0), d(2,1)
    3.0, 2.5, 1.0, // d(3,0), d(3,1), d(3,2)
];
// Compute persistent homology up to dimension 1
let barcodes = ripser(&distances, 4, 1, None);
for bar in &barcodes {
    println!("H{}: [{}, {})", bar.dim, bar.birth, bar.death);
}
// Filter by dimension
let h0 = filter_by_dim(&barcodes, 0);
let h1 = filter_by_dim(&barcodes, 1);
```

## 128-bit Index Mode

If you encounter this error:

```
libc++abi: terminating due to uncaught exception of type std::overflow_error:
simplex index 9239466805310005056 in filtration is larger than maximum index 36028797018963967
```

your dataset is too large for the default 64-bit index mode. Use `ripser128` instead:

```rust
use ripser::{ripser128, Barcode, filter_by_dim};
let barcodes = ripser128(&distances, n, max_dim, None);
```

`ripser128` has the same API as `ripser` but uses 128-bit simplex indices, supporting much larger point clouds at the cost of slower computation.

> **Platform note:** `ripser128` is only available on Linux and macOS. It relies on `__int128`, which is not supported by MSVC on Windows. The function is gated behind `#[cfg(not(target_os = "windows"))]` and will not be compiled or available on Windows targets.

## API

### Main functions

```rust
pub fn ripser(
    distances: &[f32],  // Lower triangular distance matrix
    n: usize,           // Number of points
    max_dim: i32,       // Maximum homology dimension
    threshold: Option<f32>, // Maximum filtration value (None for no limit)
) -> Vec<Barcode>

#[cfg(not(target_os = "windows"))]
pub fn ripser128(
    distances: &[f32],  // Lower triangular distance matrix
    n: usize,           // Number of points
    max_dim: i32,       // Maximum homology dimension
    threshold: Option<f32>, // Maximum filtration value (None for no limit)
) -> Vec<Barcode>
```

### Barcode struct

```rust
pub struct Barcode {
    pub dim: i32,    // Homological dimension (0, 1, 2, ...)
    pub birth: f32,  // Birth time (filtration value)
    pub death: f32,  // Death time (f32::INFINITY for essential features)
}
```

### Helper functions

- `filter_by_dim(barcodes, dim)` - Filter barcodes by dimension
- `infinite_bars(barcodes)` - Get only infinite (essential) bars
- `finite_bars(barcodes)` - Get only finite bars
- `full_to_lower_triangular(matrix, n)` - Convert full matrix to lower triangular
- `euclidean_distance_matrix(points, n, dim)` - Compute distances from point cloud

## Distance Matrix Format

The lower triangular distance matrix is stored in row-major order:

```
     0   1   2   3
0    -
1  d10   -
2  d20 d21   -
3  d30 d31 d32   -
```

Flattened: `[d10, d20, d21, d30, d31, d32, ...]`
Length: `n * (n - 1) / 2`

## Building

```bash
cargo build --release
```

## Running Tests

```bash
cargo test
```

## Running Examples

```bash
cargo run --example basic
```

## License

MIT License (same as original Ripser)

## Acknowledgments

This is a Rust wrapper for Ripser by Ulrich Bauer:
https://github.com/Ripser/ripser
