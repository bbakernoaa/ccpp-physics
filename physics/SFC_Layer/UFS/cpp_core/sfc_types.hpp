#ifndef SFC_TYPES_HPP
#define SFC_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace sfc {

// Define precision alias used throughout the surface layer math solvers (FR-010)
using Real = double;

// Standard View (Read-Write, 1D)
using View1D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 1D)
using ConstView1D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent>, std::layout_left>;

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Atmospheric Surface Sounding coordinate structures for SFC diffusion
struct SurfaceSounding {
    ConstView1D u1;           // Zonal wind component (m/s)
    ConstView1D v1;           // Meridional wind component (m/s)
    ConstView1D t1;           // Temperature (K)
    ConstView1D q1;           // Specific humidity (kg/kg)
    ConstView1D z1;           // Boundary height (m)
    ConstView1D ps;           // Surface pressure (Pa)
    ConstView1D tskin;        // Surface skin temperature (K)
};

} // namespace sfc

#endif // SFC_TYPES_HPP
