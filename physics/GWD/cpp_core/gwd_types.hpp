#ifndef GWD_TYPES_HPP
#define GWD_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace gwd {

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 3D View (Read-Write)
using View3D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 3D)
using ConstView3D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Atmospheric Sounding Profile for GWD
struct GwdSounding {
    ConstView2D ugrs;     // Zonal wind component (m/s)
    ConstView2D vgrs;     // Meridional wind component (m/s)
    ConstView2D tgrs;     // Air temperature (K)
    ConstView2D q1;       // Specific humidity profile (kg/kg)
    ConstView2D prsl;     // Mean layer pressure profile (Pa)
    ConstView2D prsi;     // Interface layer pressure profile (Pa)
    ConstView2D prslk;    // Exner function profile
    ConstView2D phil;     // Geopotential height at layer centers (m)
    ConstView2D phii;     // Geopotential height at layer interfaces (m)
    ConstView2D del;      // Layer pressure thickness profile (Pa)
};

} // namespace gwd

#endif // GWD_TYPES_HPP
