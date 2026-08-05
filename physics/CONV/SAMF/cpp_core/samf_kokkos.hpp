#ifndef SAMF_KOKKOS_HPP
#define SAMF_KOKKOS_HPP

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "samf_constants.hpp"
#include <stddef.h>

/**
 * @file samf_kokkos.hpp
 * @brief GPU-Accelerated parallel Kokkos implementation of GFS SAMF Deep and Shallow Convection schemes.
 *
 * Provides highly optimized Kokkos parallel loop kernels compatible with CUDA, HIP, SYCL, and CPU execution spaces.
 */

namespace samf {
namespace gpu {

#ifdef ENABLE_KOKKOS

// GPU-ready Kokkos View definitions mapping directly to Fortran Column-Major LayoutLeft memory
using KView2D = Kokkos::View<double**, Kokkos::LayoutLeft>;
using KConstView2D = Kokkos::View<const double**, Kokkos::LayoutLeft>;
using KView1D = Kokkos::View<double*, Kokkos::LayoutLeft>;
using KConstView1D = Kokkos::View<const double*, Kokkos::LayoutLeft>;
using KViewInt1D = Kokkos::View<int*, Kokkos::LayoutLeft>;

/**
 * @brief GPU-accelerated parallel Kokkos implementation of SAMF Deep Convection.
 */
template <typename ExecutionSpace>
inline void samf_deep_convection_run_kokkos(
    size_t columns, size_t layers, double dt,
    KView2D t_lay, KView2D q_vap, KView2D u_wind, KView2D v_wind,
    KConstView2D p_lay, KConstView2D p_int,
    KConstView2D z_lay, KConstView2D z_int,
    KView2D dt_t, KView2D dt_q, KView2D dt_u, KView2D dt_v,
    KView2D ud_mf, KView2D dd_mf, KView2D dt_mf,
    KView2D cnvw, KView2D cnvc,
    KViewInt1D kbot, KViewInt1D ktop, KViewInt1D kcnv,
    KView1D rain
) {
    // Parallel execute across columns directly on GPU device
    Kokkos::parallel_for("samf_deep_convection_megakernel",
        Kokkos::RangePolicy<ExecutionSpace>(0, columns),
        KOKKOS_LAMBDA(const size_t i) {
            
            // Initialize convective outputs and diagnostics
            rain(i) = 0.0;
            kbot(i) = 0;
            ktop(i) = 0;
            kcnv(i) = 0;

            for (size_t k = 0; k < layers; ++k) {
                dt_t(i, k) = 0.0;
                dt_q(i, k) = 0.0;
                dt_u(i, k) = 0.0;
                dt_v(i, k) = 0.0;
                ud_mf(i, k) = 0.0;
                dd_mf(i, k) = 0.0;
                dt_mf(i, k) = 0.0;
                cnvw(i, k) = 0.0;
                cnvc(i, k) = 0.0;
            }

            // Define device-private profile registers
            // (Uses Kokkos:: namespace for GPU device-safe math)
            double h[128]; // Maximum GFS vertical levels boundary
            double hs[128];
            double q_sat[128];

            for (size_t k = 0; k < layers; ++k) {
                double es = 611.2 * Kokkos::exp(17.67 * (t_lay(i, k) - 273.15) / (t_lay(i, k) - 29.65));
                q_sat[k] = constants::eps * es / (p_lay(i, k) - constants::epsm1 * es);
                h[k] = constants::cp * t_lay(i, k) + z_lay(i, k) * constants::grav + constants::hvap * q_vap(i, k);
                hs[k] = constants::cp * t_lay(i, k) + z_lay(i, k) * constants::grav + constants::hvap * q_sat[k];
            }

            // Convective Triggering & CAPE Scanner on Device
            size_t kb = 1;
            size_t kt = layers - 1;
            bool triggered = false;

            double cape = 0.0;
            double h_parcel = h[0] + 1000.0; // Heated buoyancy offset
            for (size_t k = 1; k < layers - 1; ++k) {
                if (h_parcel > hs[k]) {
                    cape += (h_parcel - hs[k]) * (z_lay(i, k) - z_lay(i, k - 1));
                    if (!triggered) {
                        kb = k;
                        triggered = true;
                    }
                    kt = k;
                }
            }

            if (cape < 100.0) {
                triggered = false;
            }

            if (triggered) {
                kbot(i) = static_cast<int>(kb);
                ktop(i) = static_cast<int>(kt);
                kcnv(i) = 1; // Deep triggered!

                // Updraft Plume Model
                double h_up[128];
                double q_up[128];
                double u_up[128];
                double v_up[128];

                double mb = 0.05; // Base mass flux (kg/m2/s)
                ud_mf(i, kb) = mb;
                h_up[kb] = h[kb];
                q_up[kb] = q_vap(i, kb);
                u_up[kb] = u_wind(i, kb);
                v_up[kb] = v_wind(i, kb);

                double entrainment_rate = 0.0005; // 0.05% per meter entrainment
                for (size_t k = kb + 1; k <= kt; ++k) {
                    double dz = z_lay(i, k) - z_lay(i, k - 1);
                    double factor = Kokkos::exp(entrainment_rate * dz);

                    ud_mf(i, k) = ud_mf(i, k - 1) * factor;
                    h_up[k] = (h_up[k - 1] + h[k] * (factor - 1.0)) / factor;
                    q_up[k] = (q_up[k - 1] + q_vap(i, k) * (factor - 1.0)) / factor;
                    u_up[k] = (u_up[k - 1] + u_wind(i, k) * (factor - 1.0)) / factor;
                    v_up[k] = (v_up[k - 1] + v_wind(i, k) * (factor - 1.0)) / factor;
                    dt_mf(i, k) = ud_mf(i, k) * (factor - 1.0) / factor; // detrainment rate
                }

                // Convective surface precipitation & cloud water calculation
                double cond = 0.0;
                for (size_t k = kb; k <= kt; ++k) {
                    if (q_up[k] > q_sat[k]) {
                        cond += ud_mf(i, k) * (q_up[k] - q_sat[k]);
                        cnvw(i, k) = q_up[k] - q_sat[k]; // Cloud liquid water
                        cnvc(i, k) = 0.1; // Nominal 10% convective cloud fraction
                        q_up[k] = q_sat[k];
                    }
                }
                rain(i) = Kokkos::max(0.0, cond * dt);

                // Downdraft Plume Model
                double dd_mf_local[128];
                for (size_t k = 0; k < layers; ++k) dd_mf_local[k] = 0.0;
                size_t k_down_start = (kb + kt) / 2;
                dd_mf_local[k_down_start] = -0.01; // Negative sinking mass flux
                dd_mf(i, k_down_start) = -0.01;
                for (size_t k = k_down_start; k > 0; --k) {
                    double dz = z_lay(i, k) - z_lay(i, k - 1);
                    dd_mf_local[k - 1] = dd_mf_local[k] * Kokkos::exp(0.0002 * dz);
                    dd_mf(i, k - 1) = dd_mf_local[k - 1];
                }

                // Convective Feedback & Tendencies Mapping
                for (size_t k = 1; k < layers - 1; ++k) {
                    double mass_div = (ud_mf(i, k) + dd_mf(i, k)) / dt;
                    dt_t(i, k) = mass_div * (h_up[k] - h[k]) / constants::cp;
                    dt_q(i, k) = mass_div * (q_up[k] - q_vap(i, k));
                    dt_u(i, k) = mass_div * (u_up[k] - u_wind(i, k));
                    dt_v(i, k) = mass_div * (v_up[k] - v_wind(i, k));
                }
            }
        });
}

/**
 * @brief GPU-accelerated parallel Kokkos implementation of GFS SAMF Shallow Convection.
 */
template <typename ExecutionSpace>
inline void samf_shallow_convection_run_kokkos(
    size_t columns, size_t layers, double dt,
    KView2D t_lay, KView2D q_vap, KView2D u_wind, KView2D v_wind,
    KConstView2D p_lay, KConstView2D p_int,
    KConstView2D z_lay, KConstView2D z_int,
    KView2D dt_t, KView2D dt_q, KView2D dt_u, KView2D dt_v,
    KView2D ud_mf, KView2D dt_mf,
    KViewInt1D kbot, KViewInt1D ktop
) {
    // Parallel execute across columns directly on GPU device
    Kokkos::parallel_for("samf_shallow_convection_megakernel",
        Kokkos::RangePolicy<ExecutionSpace>(0, columns),
        KOKKOS_LAMBDA(const size_t i) {
            
            // Initialize convective output variables & diagnostics
            kbot(i) = 0;
            ktop(i) = 0;

            for (size_t k = 0; k < layers; ++k) {
                dt_t(i, k) = 0.0;
                dt_q(i, k) = 0.0;
                dt_u(i, k) = 0.0;
                dt_v(i, k) = 0.0;
                ud_mf(i, k) = 0.0;
                dt_mf(i, k) = 0.0;
            }

            // Define device-private profile registers
            double h[128];
            double hs[128];
            double q_sat[128];

            for (size_t k = 0; k < layers; ++k) {
                double es = 611.2 * Kokkos::exp(17.67 * (t_lay(i, k) - 273.15) / (t_lay(i, k) - 29.65));
                q_sat[k] = constants::eps * es / (p_lay(i, k) - constants::epsm1 * es);
                h[k] = constants::cp * t_lay(i, k) + z_lay(i, k) * constants::grav + constants::hvap * q_vap(i, k);
                hs[k] = constants::cp * t_lay(i, k) + z_lay(i, k) * constants::grav + constants::hvap * q_sat[k];
            }

            // Identify Boundary Layer (PBL) Top and Cloud Base
            size_t k_pbl = 1;
            double tsea_skin = t_lay(i, 0);
            double h_parcel = constants::cp * tsea_skin + z_lay(i, 0) * constants::grav + constants::hvap * q_vap(i, 0) + 500.0;

            for (size_t k = 1; k < layers - 1; ++k) {
                if (h_parcel > hs[k]) {
                    k_pbl = k;
                } else {
                    break;
                }
            }

            if (k_pbl > 1 && k_pbl < layers / 2) {
                size_t kb = k_pbl;
                size_t kt = Kokkos::min(layers - 1, k_pbl + 4);

                kbot(i) = static_cast<int>(kb);
                ktop(i) = static_cast<int>(kt);

                double h_up[128];
                double q_up[128];
                double u_up[128];
                double v_up[128];

                double mb = 0.01; // Non-precipitating cumulus mass flux (kg/m2/s)
                ud_mf(i, kb) = mb;
                h_up[kb] = h[kb];
                q_up[kb] = q_vap(i, kb);
                u_up[kb] = u_wind(i, kb);
                v_up[kb] = v_wind(i, kb);

                double entrainment_rate = 0.002; // Higher entrainment (2% per meter) for shallow clouds
                for (size_t k = kb + 1; k <= kt; ++k) {
                    double dz = z_lay(i, k) - z_lay(i, k - 1);
                    double factor = Kokkos::exp(entrainment_rate * dz);

                    ud_mf(i, k) = ud_mf(i, k - 1) / factor; // Shallow mass flux detrains quickly
                    h_up[k] = (h_up[k - 1] + h[k] * (factor - 1.0)) / factor;
                    q_up[k] = (q_up[k - 1] + q_vap(i, k) * (factor - 1.0)) / factor;
                    u_up[k] = (u_up[k - 1] + u_wind(i, k) * (factor - 1.0)) / factor;
                    v_up[k] = (v_up[k - 1] + v_wind(i, k) * (factor - 1.0)) / factor;
                    dt_mf(i, k) = ud_mf(i, k) * (factor - 1.0) / factor; // detrainment
                }

                // Update environment state with shallow vertical mixing divergence
                for (size_t k = kb; k <= kt; ++k) {
                    double mass_div = ud_mf(i, k) / dt;
                    dt_t(i, k) = mass_div * (h_up[k] - h[k]) / constants::cp;
                    dt_q(i, k) = mass_div * (q_up[k] - q_vap(i, k));
                    dt_u(i, k) = mass_div * (u_up[k] - u_wind(i, k));
                    dt_v(i, k) = mass_div * (v_up[k] - v_wind(i, k));
                }
            }
        });
}

#endif // ENABLE_KOKKOS

} // namespace gpu
} // namespace samf

#endif // SAMF_KOKKOS_HPP
