#ifndef SAMF_TYPES_HPP
#define SAMF_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace samf {

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 3D View (Read-Write)
using View3D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 3D)
using ConstView3D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 2D Integer View (Read-Write)
using IntView2D = std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant 2D Integer View (Read-Only)
using ConstIntView2D = std::mdspan<const int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Physical State for Convective Sounding Input Arrays
struct ConvectiveSounding {
    View2D t_lay;         // Temperature (K)
    View2D q_vap;         // Specific humidity / water vapor (kg/kg)
    View2D u_wind;        // Zonal wind component (m/s)
    View2D v_wind;        // Meridional wind component (m/s)
    ConstView2D p_lay;    // Mean layer pressure (Pa)
    ConstView2D p_int;    // Pressure at layer interfaces (Pa)
    ConstView2D z_lay;    // Geopotential height of layer centers (m)
    ConstView2D z_int;    // Geopotential height of layer interfaces (m)
};

} // namespace samf

#endif // SAMF_TYPES_HPP
