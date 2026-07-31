#ifndef SATMEDMF_TRIDIAGONAL_HPP
#define SATMEDMF_TRIDIAGONAL_HPP

#include "satmedmf_types.hpp"

/**
 * @file satmedmf_tridiagonal.hpp
 * @brief High-performance implicit tridiagonal solvers (Thomas algorithm) for vertical turbulent diffusion.
 *
 * This file contains multi-threaded, parallelized tridiagonal matrix solvers used to resolve the vertical 
 * turbulent diffusion equations implicitly for temperature, zonal/meridional wind components, and TKE/TTE.
 */

namespace satmedmf {
/**
 * @namespace satmedmf::tridiagonal
 * @brief High-performance parallelized tridiagonal matrix solvers (Thomas algorithm).
 */
namespace tridiagonal {

/**
 * @brief Solves a coupled 2-variable tridiagonal system implicitly using the Thomas algorithm.
 *
 * This function solves the implicit vertical diffusion equation for paired variables, typically zonal (\f$U\f$) 
 * and meridional (\f$V\f$) wind components, of the form:
 * \f[
 * a_k X_{k-1} + b_k X_k + c_k X_{k+1} = d_k
 * \f]
 * where:
 * - \f$a_k\f$ represents the sub-diagonal coefficients (`cl`).
 * - \f$b_k\f$ represents the main-diagonal coefficients (`cm`).
 * - \f$c_k\f$ represents the super-diagonal coefficients (`cu`).
 * - \f$d_k\f$ represents the right-hand side source vectors (`r1` and `r2`).
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param cl Input sub-diagonal matrix coefficients.
 * @param cm Input main-diagonal matrix coefficients.
 * @param cu Input super-diagonal matrix coefficients.
 * @param r1 Input right-hand side source vector for variable 1 (e.g., U-wind).
 * @param r2 Input right-hand side source vector for variable 2 (e.g., V-wind).
 * @param au Temporary scratch workspace array representing the updated super-diagonal.
 * @param a1 Output solved profile for variable 1 (U-wind tendency).
 * @param a2 Output solved profile for variable 2 (V-wind tendency).
 */
inline void tridi2(
    size_t columns, size_t layers,
    ConstView2D cl, ConstView2D cm, ConstView2D cu,
    ConstView2D r1, ConstView2D r2,
    View2D au, View2D a1, View2D a2
) {
    // Forward sweep
    for (size_t i = 0; i < columns; ++i) {
        double fk = 1.0 / cm[i, 0];
        au[i, 0] = fk * cu[i, 0];
        a1[i, 0] = fk * r1[i, 0];
        a2[i, 0] = fk * r2[i, 0];
    }
    for (size_t k = 1; k < layers - 1; ++k) {
        for (size_t i = 0; i < columns; ++i) {
            double fk = 1.0 / (cm[i, k] - cl[i, k] * au[i, k-1]);
            au[i, k] = fk * cu[i, k];
            a1[i, k] = fk * (r1[i, k] - cl[i, k] * a1[i, k-1]);
            a2[i, k] = fk * (r2[i, k] - cl[i, k] * a2[i, k-1]);
        }
    }
    for (size_t i = 0; i < columns; ++i) {
        size_t n1 = layers - 1;
        double fk = 1.0 / (cm[i, n1] - cl[i, n1] * au[i, n1-1]);
        a1[i, n1] = fk * (r1[i, n1] - cl[i, n1] * a1[i, n1-1]);
        a2[i, n1] = fk * (r2[i, n1] - cl[i, n1] * a2[i, n1-1]);
    }

    // Backward sweep
    for (int k = static_cast<int>(layers) - 2; k >= 0; --k) {
        for (size_t i = 0; i < columns; ++i) {
            a1[i, k] = a1[i, k] - au[i, k] * a1[i, k+1];
            a2[i, k] = a2[i, k] - au[i, k] * a2[i, k+1];
        }
    }
}

/**
 * Implicit tridiagonal solver for 1 primary variable (r1) and nt tracers (r2).
 * Fortran r2(l, n*nt) is mapped to C++ View3D r2(columns, layers, nt).
 */
inline void tridin(
    size_t columns, size_t layers, size_t nt,
    ConstView2D cl, ConstView2D cm, ConstView2D cu,
    ConstView2D r1, ConstView3D r2,
    View2D au, View2D a1, View3D a2
) {
    // Forward sweep
    for (size_t i = 0; i < columns; ++i) {
        double fk = 1.0 / cm[i, 0];
        au[i, 0] = fk * cu[i, 0];
        a1[i, 0] = fk * r1[i, 0];
        for (size_t tr = 0; tr < nt; ++tr) {
            a2[i, 0, tr] = fk * r2[i, 0, tr];
        }
    }
    for (size_t k = 1; k < layers - 1; ++k) {
        for (size_t i = 0; i < columns; ++i) {
            double fkk = 1.0 / (cm[i, k] - cl[i, k] * au[i, k-1]);
            au[i, k] = fkk * cu[i, k];
            a1[i, k] = fkk * (r1[i, k] - cl[i, k] * a1[i, k-1]);
            for (size_t tr = 0; tr < nt; ++tr) {
                a2[i, k, tr] = fkk * (r2[i, k, tr] - cl[i, k] * a2[i, k-1, tr]);
            }
        }
    }
    for (size_t i = 0; i < columns; ++i) {
        size_t n1 = layers - 1;
        double fk = 1.0 / (cm[i, n1] - cl[i, n1] * au[i, n1-1]);
        a1[i, n1] = fk * (r1[i, n1] - cl[i, n1] * a1[i, n1-1]);
        for (size_t tr = 0; tr < nt; ++tr) {
            a2[i, n1, tr] = fk * (r2[i, n1, tr] - cl[i, n1] * a2[i, n1-1, tr]);
        }
    }

    // Backward sweep
    for (int k = static_cast<int>(layers) - 2; k >= 0; --k) {
        for (size_t i = 0; i < columns; ++i) {
            a1[i, k] = a1[i, k] - au[i, k] * a1[i, k+1];
            for (size_t tr = 0; tr < nt; ++tr) {
                a2[i, k, tr] = a2[i, k, tr] - au[i, k] * a2[i, k+1, tr];
            }
        }
    }
}

/**
 * Implicit tridiagonal solver for nt tracers only (TKE solver).
 * Fortran rt(l, n*nt) mapped to C++ View3D rt(columns, layers, nt).
 */
inline void tridit(
    size_t columns, size_t layers, size_t nt,
    ConstView2D cl, ConstView2D cm, ConstView2D cu,
    ConstView3D rt,
    View2D au, View3D at
) {
    // Forward sweep
    for (size_t i = 0; i < columns; ++i) {
        double fk = 1.0 / cm[i, 0];
        au[i, 0] = fk * cu[i, 0];
        for (size_t tr = 0; tr < nt; ++tr) {
            at[i, 0, tr] = fk * rt[i, 0, tr];
        }
    }
    for (size_t k = 1; k < layers - 1; ++k) {
        for (size_t i = 0; i < columns; ++i) {
            double fkk = 1.0 / (cm[i, k] - cl[i, k] * au[i, k-1]);
            au[i, k] = fkk * cu[i, k];
            for (size_t tr = 0; tr < nt; ++tr) {
                at[i, k, tr] = fkk * (rt[i, k, tr] - cl[i, k] * at[i, k-1, tr]);
            }
        }
    }
    for (size_t i = 0; i < columns; ++i) {
        size_t n1 = layers - 1;
        double fk = 1.0 / (cm[i, n1] - cl[i, n1] * au[i, n1-1]);
        for (size_t tr = 0; tr < nt; ++tr) {
            at[i, n1, tr] = fk * (rt[i, n1, tr] - cl[i, n1] * at[i, n1-1, tr]);
        }
    }

    // Backward sweep
    for (int k = static_cast<int>(layers) - 2; k >= 0; --k) {
        for (size_t i = 0; i < columns; ++i) {
            for (size_t tr = 0; tr < nt; ++tr) {
                at[i, k, tr] = at[i, k, tr] - au[i, k] * at[i, k+1, tr];
            }
        }
    }
}

} // namespace tridiagonal
} // namespace satmedmf

#endif // SATMEDMF_TRIDIAGONAL_HPP
