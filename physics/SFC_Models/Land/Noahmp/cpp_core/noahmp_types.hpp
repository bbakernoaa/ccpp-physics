#ifndef NOAHMP_TYPES_HPP
#define NOAHMP_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace noahmp {

// Define precision alias used throughout the Noah-MP solvers
using Real = double;

// Standard double-precision View types matching Column-Major layouts (using mdspan)
using View1D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent>, std::layout_left>;
using ConstView1D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent>, std::layout_left>;
using View2D = std::mdspan<Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;
using ConstView2D = std::mdspan<const Real, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Radiative Sounding coordinates for GFS Noah-MP (using mdspan)
struct LandSounding {
    ConstView2D stc;           // Soil temperature profile [columns, soil_layers]
    ConstView2D smc;           // Total soil volumetric moisture profile [columns, soil_layers]
    ConstView2D sh2o;          // Liquid soil volumetric moisture profile [columns, soil_layers]
    ConstView2D sldpst;        // Soil layer depth profile [columns, soil_layers]
    ConstView1D tg;            // Ground surface temperature [columns]
    ConstView1D tv;            // Canopy vegetation temperature [columns]
};

/**
 * @brief Configuration parameter options struct containing all 19 Noah-MP scientific options.
 */
struct NoahMP_Config {
    int idveg;      // Dynamic vegetation option
    int iopt_crs;   // Canopy stomatal resistance option
    int iopt_btr;   // Soil moisture factor for stomatal resistance option
    int iopt_run;   // Runoff and groundwater option
    int iopt_sfc;   // Surface layer drag coeff option
    int iopt_frz;   // Supercooled liquid water option
    int iopt_inf;   // Frozen soil permeability option
    int iopt_rad;   // Radiation transfer option
    int iopt_alb;   // Ground snow surface albedo option
    int iopt_snf;   // Partitioning precipitation into rainfall & snowfall option
    int iopt_tbot;  // Lower boundary condition of soil temperature option
    int iopt_stc;   // Snow/soil temperature time scheme option
    int iopt_trs;   // Thermal roughness scheme option
    int iopt_diag;  // Surface diagnose approach option
    int iopt_rsf;   // Surface resistance option
    int iopt_soil;  // Soil parameter treatment option
    int iopt_pedo;  // Pedotransfer function option
    int iopt_crop;  // Crop model option
    int iopt_gla;   // Glacier treatment option
    int iopt_z0m;   // z0m treatment option
};

} // namespace noahmp

#endif // NOAHMP_TYPES_HPP
