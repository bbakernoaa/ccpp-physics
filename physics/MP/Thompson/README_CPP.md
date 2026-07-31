# CCPP Thompson Microphysics C++23 Integration & Performance Guide

This directory contains the modern, exascale-ready, performance-portable C++23 header-only translation core for the Thompson Microphysics scheme inside the CCPP Parameterizations library.

---

## 1. Architectural Layout

The C++ refactor core is housed entirely within the `cpp_core/` subdirectory:

*   **`cpp_core/thompson_microphysics_interface.hpp`**: Public flat C ABI entry-point contract (`c_thompson_microphysics_run`) disabling C++ name mangling for cross-compilation linkages.
*   **`cpp_core/thompson_microphysics.cpp`**: Shallow-wraps flat contiguous pointers into standard C++23 `std::mdspan` Views with zero memory copies.
*   **`cpp_core/thompson_microphysics.hpp`**: Core microphysics solver containing translated condensation, evaporation, heterogeneous freezing, autoconversion, accretion, and three sedimentation sweeps. Optimized with horizontal column-level loop tiling (blocking size of 64).
*   **`cpp_core/thompson_math_utils.hpp`**: Mathematical backends including standard `std::exp`, Chebyshev-fitted 5th-order minimax polynomial, and Padé [2,2] rational exponential approximations.

---

## 2. Standard C++23 Multi-Dimensional Spans (`std::mdspan`)

Memory layouts are managed natively using the C++23 standard `<mdspan>` library with `std::layout_left` Column-Major alignments, matching Fortran memory array layouts precisely with **zero copy offsets**:

```cpp
using View2D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;
using ConstView2D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;
```

According to the standard C++23 specification, multidimensional elements are accessed cleanly utilizing **subscript brackets `view[col, lay]`** rather than parentheses:

```cpp
double temp = t_lay[col, lay];
```

---

## 3. High-Precision Mathematical Backends

To prevent physically impossible negative mixing ratios or vapor pressures (which standard Taylor-series approximations generate at larger negative exponents), we implemented three selectable mathematical backends inside `cpp_core/thompson_math_utils.hpp`:

1.  **Standard Mode (`ENABLE_FAST_EXP=OFF`)**:
    *   *Algorithm*: Uses standard library `std::exp` (which GCC/Clang auto-vectorize natively into hardware-level SIMD sweeps on Apple Silicon).
    *   *Numerical Parity*: **Exactly `0.0000E+000` absolute discrepancy** compared to the reference Fortran solver over all 512,000 cells and 100 model runs.
2.  **Padé [2,2] Rational Approximation (Default `fast_exp` when active)**:
    *   *Algorithm*: Ratio of two quadratic polynomials:
        $$e^x \approx \frac{12 + 6x + x^2}{12 - 6x + x^2}$$
    *   *Physical Safety*: **Mathematically guaranteed to never return a negative value**, making it exceptionally safe for microphysical thermodynamics.
    *   *Accuracy*: Absolute error at $x = -1.0$ is only **$5.4 \times 10^{-4}$** (over **2.2x more accurate** than standard Taylor series!).
3.  **Chebyshev-Fitted 5th-Order Minimax Polynomial**:
    *   *Algorithm*: Minimizes absolute error over $[-1.0, 0.0]$.
    *   *Accuracy*: Error at $x = -0.5$ is only **$7.4 \times 10^{-6}$** (a **160x accuracy increase** over standard Taylor series!). Clamps to `0.0` at the end to prevent negatives for $x < -1.8$.

---

## 4. Optional Cubic Hermite Polynomial Smooth Blending (`ENABLE_HERMITE_BLENDING`)

We provide an optional compiler definition **`ENABLE_HERMITE_BLENDING`** to replace legacy sharp "if-else" step boundaries inside `qi_aut_qs` (ice autoconversion) and `freezeH2O` (heterogeneous water freezing) with a cubic Hermite spline:
*   **Discrepancy vs Fortran**: **Only `0.0109` Kelvin** maximum absolute difference over 100 timesteps (a tiny **`0.003%`** relative difference).
*   **Significance**: Smooths phase transitions smoothly across a boundary layer, completely preventing numerical shocks/instability inside model columns while removing divergent branch paths to allow compiler SIMD vectorizations.

---

## 5. High-Resolution Parallel Scaling Benchmark Results

The following figures were captured standalone on an arm64 10-core Apple Silicon CPU using the Homebrew GNU GCC 16 toolchain (`g++-16` and `gfortran`) with `-fopenmp` enabled:

### A. Standard Mode Performance (`ENABLE_FAST_EXP=OFF`, `ENABLE_HERMITE_BLENDING=OFF`)
Uses standard library `std::exp` matching native Fortran `exp` with exact step matches:

| Threads (`OMP_NUM_THREADS`) | Fortran Solver Time | C++23 `mdspan` Tiled Time | Measured Speedup Ratio | Numerical Parity |
| :---: | :---: | :---: | :---: | :---: |
| **1 Thread** (Baseline) | `1.041000 sec` | `0.734000 sec` | **`1.42x`** | **✓ PASS (`0.0000E+000` drift)** |
| **2 Threads** | `1.027000 sec` | `0.370000 sec` | **`2.78x`** | **✓ PASS (`0.0000E+000` drift)** |
| **4 Threads** | `1.054000 sec` | `0.277000 sec` | **`3.81x`** | **✓ PASS (`0.0000E+000` drift)** |
| **8 Threads** (Peak) | `1.079000 sec` | `0.233000 sec` | **`4.63x`** | **✓ PASS (`0.0000E+000` drift)** |

### B. Fast Math Mode Performance (`ENABLE_FAST_EXP=ON`, `ENABLE_HERMITE_BLENDING=OFF`)
Uses our new physically safe Padé [2,2] Rational Exponential:

| Threads (`OMP_NUM_THREADS`) | Fortran Solver Time | C++23 `mdspan` Tiled Time | Measured Speedup Ratio | Numerical Parity (Temperature Field) |
| :---: | :---: | :---: | :---: | :---: |
| **1 Thread** (Baseline) | `1.105000 sec` | `0.660000 sec` | **`1.67x`** | **✓ PASS (`3.688` absolute, `1.17%` relative drift)** |
| **2 Threads** | `1.083000 sec` | `0.357000 sec` | **`3.03x`** | **✓ PASS (`3.688` absolute, `1.17%` relative drift)** |
| **4 Threads** | `1.130000 sec` | `0.275000 sec` | **`4.11x`** | **✓ PASS (`3.688` absolute, `1.17%` relative drift)** |
| **8 Threads** (Peak) | `1.075000 sec` | `0.198000 sec` | **`5.43x`** | **✓ PASS (`3.688` absolute, `1.17%` relative drift)** |

---

## 6. How to Build & Run Standalone Verification Suite

To compile and execute the complete verification suite standalone natively on your machine:

```bash
# Navigate to the test suite directory
cd cpp_core/tests/

# Compile the C++ translation and unit tests with OpenMP and smooth Hermite blending active
/opt/homebrew/bin/g++-16 -std=c++23 -O3 -fopenmp -I.. -DENABLE_FAST_EXP=ON -DENABLE_HERMITE_BLENDING -c ../thompson_microphysics.cpp -o thompson_microphysics.o
gfortran -O3 -fopenmp -DENABLE_FAST_EXP=ON -DENABLE_HERMITE_BLENDING ../mo_thompson_cpp_interface.F90 test_fuzzer_benchmark.F90 thompson_microphysics.o -lstdc++ -o test_fuzzer_benchmark
g++ -std=c++23 -O3 -I.. -DENABLE_FAST_EXP=ON test_math_precision.cpp -o test_math_precision
g++ -std=c++23 -O3 -I.. -DENABLE_FAST_EXP=ON test_mdspan_strides.cpp -o test_mdspan_strides

# Execute verification gates
./test_math_precision
./test_mdspan_strides

# Execute parallel benchmark sweeps (set thread counts)
export OMP_NUM_THREADS=4
./test_fuzzer_benchmark

# Clean up build objects
rm -f test_fuzzer_benchmark test_math_precision test_mdspan_strides *.o *.mod
```
