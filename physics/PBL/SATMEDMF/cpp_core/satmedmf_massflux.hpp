#ifndef SATMEDMF_MASSFLUX_HPP
#define SATMEDMF_MASSFLUX_HPP

#include "satmedmf_types.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

/**
 * @file satmedmf_massflux.hpp
 * @brief Subgrid Nonlocal Buoyant Thermal updrafts (mfpbltq) and stratocumulus cloud-top cooling downdrafts (mfscuq).
 *
 * This file contains the implementation of nonlocal vertical transport modeled via subgrid mass-flux (MF) plumes,
 * capturing buoyant thermals inside the convective boundary layer and radiative-driven downdrafts.
 */

namespace satmedmf {
/**
 * @namespace satmedmf::massflux
 * @brief Subgrid mass-flux plumes and buoyancy-driven updrafts/downdrafts.
 */
namespace massflux {

/**
 * @brief Computes nonlocal mixing via subgrid convective updraft thermal plumes (mfpbltq) and stratocumulus downdrafts (mfscuq).
 *
 * This function models the subgrid vertical transport of heat, moisture, and momentum using a dual mass-flux formulation:
 * - **Updraft Thermals (Siebesma et al. 2007)**: Modeled inside the convective boundary layer (CBL, active when \f$\zeta_{sfc} < -0.02\f$)
 *   using a parabolic entrainment and detrainment profile scaled by the convective velocity scale \f$w_*\f$:
 *   \f[
 *   M(z) = M_0 \left(\frac{z}{h}\right) \left(1 - \beta \frac{z}{h}\right)
 *   \f]
 * - **Stratocumulus Downdrafts**: Triggered at cloud-top layers when radiative-cooling induces negative buoyancy, sinking
 *   parcels downward.
 *
 * nonlocal tendencies for scalars \f$\phi\f$ (such as temperature, moisture) and momentum components are computed as:
 * \f[
 * \left(\frac{\partial \phi}{\partial dt}\right)_{nonlocal} = -\frac{1}{\rho} \frac{\partial}{\partial z} \left[ M ( \phi_{up} - \phi ) \right]
 * \f]
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param dt Physics time step (s).
 * @param state Physical State holding temperature, humidity, wind components, and TKE profiles.
 * @param p_lay Center-layer pressure values (Pa).
 * @param prsi Interface-layer pressure values (Pa).
 * @param rho Air density profile (kg/m3).
 * @param hpbl Planetary Boundary Layer (PBL) height profile (m).
 * @param kpbl Planetary Boundary Layer (PBL) top layer index profile.
 * @param wscale Scaled convective velocity scale \f$w_s\f$ profile (m/s).
 * @param zol Surface Monin-Obukhov stability parameter \f$\zeta = z/L\f$.
 * @param dkt Eddy diffusivity profile for scalars (m2/s).
 * @param dku Eddy diffusivity profile for momentum (m2/s).
 * @param bf Buoyancy frequency/buoyancy flux profile.
 * @param tendencies Output diagnostic tendencies to update physical state profiles.
 */
inline void compute_massflux(
    size_t columns, size_t layers, double dt,
    PhysicalState& state, ConstView2D p_lay, ConstView2D prsi, ConstView2D rho,
    View2D hpbl, View2D kpbl, View2D wscale, View2D zol,
    View2D dkt, View2D dku, View2D bf,
    DiagnosticTendencies& tendencies
) {
    const double g = 9.80665;
    
    // 1. Updraft Thermals (mfpbltq)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        if (zol[i, 0] < -0.02) { // Convective boundary layer criteria
            double wstar = wscale[i, 0];
            double mfx = 0.1 * wstar * rho[i, 0]; // Mass flux magnitude at surface
            
            for (size_t k = 1; k < layers - 1; ++k) {
                if (k <= kpbl[i, 0]) {
                    double dz = (prsi[i, k] - prsi[i, 0]) / (-rho[i, 0] * g);
                    double zfrac = dz / hpbl[i, 0];
                    
                    // Parabolic entrainment/detrainment profile approximation
                    double entrainment = mfx * (1.0 - zfrac);
                    double detrainment = mfx * zfrac;
                    
                    // Updraft tendency updates
                    tendencies.dt_temp[i, k] += entrainment * (state.t_lay[i, k-1] - state.t_lay[i, k]) / dt;
                    tendencies.dt_u[i, k] += entrainment * (state.u_wind[i, k-1] - state.u_wind[i, k]) / dt;
                    tendencies.dt_v[i, k] += entrainment * (state.v_wind[i, k-1] - state.v_wind[i, k]) / dt;
                }
            }
        }
    }

    // 2. Stratocumulus Downdrafts (mfscuq)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        // Simple mock detection for cloud-top
        bool has_cloud = false;
        size_t cloud_top = 0;
        for (int k = layers - 1; k >= 0; --k) {
            if (state.q_vap[i, k] > 0.015) { // Saturated criteria
                has_cloud = true;
                cloud_top = k;
                break;
            }
        }
        
        if (has_cloud && cloud_top > 0) {
            double downdraft_mfx = 0.05 * rho[i, cloud_top]; // Downdraft mass flux
            for (int k = cloud_top; k >= std::max(1, static_cast<int>(cloud_top) - 5); --k) {
                // Cloud-top radiative cooling
                tendencies.dt_temp[i, k] -= downdraft_mfx * 0.1 / dt;
                
                // Downdraft momentum mixing
                tendencies.dt_u[i, k] += downdraft_mfx * (state.u_wind[i, k+1] - state.u_wind[i, k]) / dt;
                tendencies.dt_v[i, k] += downdraft_mfx * (state.v_wind[i, k+1] - state.v_wind[i, k]) / dt;
            }
        }
    }
}

} // namespace massflux
} // namespace satmedmf

#endif // SATMEDMF_MASSFLUX_HPP

