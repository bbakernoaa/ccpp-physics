#ifndef RRTMG_TYPES_HPP
#define RRTMG_TYPES_HPP

#include <mdspan>
#include <stddef.h>
#include <cmath>

namespace rrtmg {

// Define precision alias used throughout the radiation math solvers (FR-010)
using Real = double;

// Highly-efficient 5th-order minimax polynomial exponential (Horner's Scheme)
#ifdef ENABLE_FAST_EXP
inline double fast_exp(double x) {
    if (x < -15.0) return 0.0;
    return 1.0 + x * (1.0 + x * (0.5 + x * (0.16666666666666667 + x * (0.041666666666666664 + x * 0.008333333333333333))));
}
#else
inline double fast_exp(double x) {
    return std::exp(x);
}
#endif

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 3D View (Read-Write)
using View3D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 3D)
using ConstView3D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Radiative Sounding Profile for SW and LW solvers
struct RadiativeSounding {
    ConstView2D t_lay;         // Temperature (K)
    ConstView2D q_vap;         // Specific humidity (kg/kg)
    ConstView2D o3_vap;        // Ozone concentration (kg/kg)
    ConstView2D cld_frac;      // Cloud fraction profile
    ConstView2D p_lay;         // Mean layer pressure (Pa)
    ConstView2D p_int;         // Pressure at layer interfaces (Pa)
    ConstView2D surface_param; // Surface albedo (SW) or emissivity (LW) per spectral band
};

} // namespace rrtmg

#endif // RRTMG_TYPES_HPP
