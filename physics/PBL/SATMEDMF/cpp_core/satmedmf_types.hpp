#ifndef SATMEDMF_TYPES_HPP
#define SATMEDMF_TYPES_HPP

#include <mdspan>
#include <stddef.h>

namespace satmedmf {

// Standard View (Read-Write, 2D)
using View2D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 2D)
using ConstView2D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 3D View (Read-Write for Tracers)
using View3D = std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant View (Read-Only, 3D for Tracers)
using ConstView3D = std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// 2D Integer View (Read-Write)
using IntView2D = std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Constant 2D Integer View (Read-Only)
using ConstIntView2D = std::mdspan<const int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left>;

// Comprehensive Physical Options & Constants Struct
struct PhysicsOptions {
    bool tte_edmf;        // Total Turbulent Energy flag (true = TTE, false = TKE)
    bool sa3dtke;         // Scale-aware 3D-TKE scheme flag (Deardorff 1980)
    bool dspheat;         // Dissipative heating temperature feedback flag
    bool use_oceanuv;     // Ocean surface currents active flag
    bool use_lpt;         // Liquid water potential temperature usage flag
    bool do_canopy;       // Sub-grid forest canopy active flag
    bool gen_tend;        // Add tendencies to state arrays flag
    bool ldiag3d;         // Save 3D physical diagnostics flag
};

// Comprehensive Physical State & Meteorological Fields
struct PhysicalState {
    View2D t_lay;         // Temperature (K)
    View2D u_wind;        // Zonal wind component (m/s)
    View2D v_wind;        // Meridional wind component (m/s)
    View2D q_vap;         // Specific humidity / water vapor (kg/kg)
    View2D te;            // TKE or TTE concentration (m2/s2)
    ConstView2D p_lay;    // Mean layer pressure (Pa)
    ConstView2D prsi;     // Pressure at layer interfaces (Pa)
    ConstView2D rho;      // Air density (kg/m3)
    
    // Boundary & Surface Fluxes
    const double* heat;   // Kinematic upward surface sensible heat flux [columns] (K m/s)
    const double* evap;   // Kinematic upward surface latent heat flux [columns] (kg/kg m/s)
    const double* stress; // Surface wind stress [columns] (m2/s2)
    const double* u10m;   // Zonal wind at 10m [columns] (m/s)
    const double* v10m;   // Meridional wind at 10m [columns] (m/s)
    const double* tsea;   // Sea surface temperature [columns] (K)
    const double* fm;     // Surface Monin-Obukhov stability function for momentum [columns]
    const double* fh;     // Surface Monin-Obukhov stability function for heat [columns]
    const double* rbsoil; // Surface Richardson number [columns]
    const double* cfch;   // Canopy forest height [columns] (m)
};

// Comprehensive Output Tendencies & Diagnostics
struct DiagnosticTendencies {
    View2D dt_temp;       // Temperature tendency (K/s)
    View2D dt_u;          // Zonal wind tendency (m/s2)
    View2D dt_v;          // Meridional wind tendency (m/s2)
    View2D dt_te;         // TKE/TTE tendency (m2/s3)
    
    // Surface tendencies and PBL heights diagnostics
    double* dusfc;        // Zonal wind surface tendency [columns] (m/s2)
    double* dvsfc;        // Meridional wind surface tendency [columns] (m/s2)
    double* dtsfc;        // Temperature surface tendency [columns] (K/s)
    double* dqsfc;        // Specific humidity surface tendency [columns] (kg/kg/s)
    double* hpbl;         // Diagnostic PBL height [columns] (m)
    double* kpbl;         // Diagnostic PBL top layer index [columns]
    double* tkeh;         // Diagnostic interface TKE array [columns * layers] (m2/s2)
    double* dkt;          // Diagnostic heat vertical diffusion coefficient [columns * layers] (m2/s)
    double* dku;          // Diagnostic momentum vertical diffusion coefficient [columns * layers] (m2/s)
};

} // namespace satmedmf

#endif // SATMEDMF_TYPES_HPP
