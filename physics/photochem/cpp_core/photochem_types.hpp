#ifndef PHOTOCHEM_TYPES_HPP
#define PHOTOCHEM_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace photochem {

// Define precision alias used throughout the photochemistry math solvers (FR-009)
using Real = double;

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 3D View (Read-Write)
using View3D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 3D)
using ConstView3D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Atmospheric Species Sounding Profile for Ozone and H2O solvers
struct PhotochemSounding {
    ConstView2D t_lay;         // Temperature (K)
    ConstView2D p_lay;         // Mean layer pressure (Pa)
    ConstView2D dp;            // Layer pressure thickness (Pa)
    ConstView2D species;       // Species concentration profile (Ozone or H2O, kg/kg)
};

} // namespace photochem

#endif // PHOTOCHEM_TYPES_HPP
