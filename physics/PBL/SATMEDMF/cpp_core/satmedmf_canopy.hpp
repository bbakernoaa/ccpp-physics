#ifndef SATMEDMF_CANOPY_HPP
#define SATMEDMF_CANOPY_HPP

#include "satmedmf_types.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

/**
 * @file satmedmf_canopy.hpp
 * @brief Subgrid Forest Canopy level setup (canopy_levs) and mass-conservative bidirectional transfer (canopy_transfer).
 *
 * This file contains the C++23 translation of the subgrid forest canopy models, solving in-canopy levels,
 * monotonic layer sorting, and mass-conservative transport of chemical tracers and tendencies.
 */

namespace satmedmf {
/**
 * @namespace satmedmf::canopy
 * @brief Subgrid forest canopy algorithms, layer setup, and mass-conservative tracer mappings.
 */
namespace canopy {

/**
 * @brief Computes subgrid canopy layers setup, inserts levels, sorts monotonically, and interpolates variables.
 *
 * This function translates `canopy_levs_run` from `canopy_levs.F90`. It performs the following sequential steps:
 * -# Inserts 3 subgrid forest canopy height levels calculated as fractions of top of canopy \f$h_{can}\f$:
 *    \f$z_{can} \in \{h_{can}, 0.5 h_{can}, 0.2 h_{can}\}\f$.
 * -# Applies an adjustment offset of \f$\delta = 0.2\f$ meters to prevent any inserted layer from collapsing too close to original model interfaces.
 * -# Sorts the combined height array of model layers and canopy levels monotonically in decreasing height order using a high-performance, bubble-sorting algorithm.
 * -# Builds location maps `kcan3` and `kmod` to trace layers across original and expanded dimensions.
 * -# Linearly interpolates wind speed, temperature, specific humidity, air density, pressure, and exchange coefficients to the inserted in-canopy layers.
 *
 * @note This routine is optimized to use stack-allocated local arrays rather than dynamic heap allocations inside columns loops,
 * preventing execution-space thread contention on multi-core architectures.
 *
 * @param im Number of horizontal grid columns.
 * @param km Number of vertical model layers.
 * @param nkc Number of inserted canopy layers (default: 3).
 * @param nkt Number of total combined layers (km + nkc).
 * @param rdgas Gas constant for dry air.
 * @param pi Mathematical constant \f$\pi\f$.
 * @param zi Interface geopotential height profile (m).
 * @param zl Layer geopotential height profile (m).
 * @param zm Geopotential interface heights of model layers (m).
 * @param prsl Center-layer pressure profile (Pa).
 * @param prsi Interface-layer pressure profile (Pa).
 * @param cfch Canopy Forest Height profile (m).
 * @param t2m 2-meter air temperature (K).
 * @param q2m 2-meter specific humidity (kg/kg).
 * @param u1 Original zonal wind profile (m/s).
 * @param v1 Original meridional wind profile (m/s).
 * @param t1 Original air temperature profile (K).
 * @param dens Original air density profile (kg/m3).
 * @param dkt Original scalar exchange coefficient profile (m2/s).
 * @param dku Original momentum exchange coefficient profile (m2/s).
 * @param kmod_data Output mapping array tracking model layers indices inside combined array.
 * @param kcan3_data Output mapping array tracking canopy layers indices inside combined array.
 * @param zmid_can_data Output combined heights array at layer centers (m).
 * @param zmom_can_data Output combined heights array at layer interfaces (m).
 * @param prsl_can_data Output combined pressure array at layer centers (Pa).
 * @param prsi_can_data Output combined pressure array at layer interfaces (Pa).
 * @param t1_can_data Output combined temperature profile (K).
 * @param qv_can_data Output combined specific humidity profile (kg/kg).
 * @param ws_can_data Output combined wind speed profile (m/s).
 * @param dens_can_data Output combined air density profile (kg/m3).
 * @param dkt_can_data Output combined scalar exchange coefficient profile (m2/s).
 * @param dku_can_data Output combined momentum exchange coefficient profile (m2/s).
 */
inline void canopy_levs_run(
    size_t im, size_t km, size_t nkc, size_t nkt,
    double rdgas, double pi,
    ConstView2D zi, ConstView2D zl, ConstView2D zm,
    ConstView2D prsl, ConstView2D prsi,
    const double* cfch,
    const double* t2m, const double* q2m,
    ConstView2D u1, ConstView2D v1, ConstView2D t1,
    ConstView2D dens, ConstView2D dkt, ConstView2D dku,
    std::vector<int>& kmod_data, std::vector<int>& kcan3_data,
    std::vector<double>& zmid_can_data, std::vector<double>& zmom_can_data,
    std::vector<double>& prsl_can_data, std::vector<double>& prsi_can_data,
    std::vector<double>& t1_can_data, std::vector<double>& qv_can_data,
    std::vector<double>& ws_can_data, std::vector<double>& dens_can_data,
    std::vector<double>& dkt_can_data, std::vector<double>& dku_can_data
) {
    const double can_frac[3] = {1.0, 0.5, 0.2};
    const double del = 0.2;
    const double min_kt = 0.1;

    View2D zmid_can(zmid_can_data.data(), im, nkt);
    View2D zmom_can(zmom_can_data.data(), im, nkt + 1);
    View2D prsl_can(prsl_can_data.data(), im, nkt);
    View2D prsi_can(prsi_can_data.data(), im, nkt + 1);
    View2D t1_can(t1_can_data.data(), im, nkt);
    View2D qv_can(qv_can_data.data(), im, nkt);
    View2D ws_can(ws_can_data.data(), im, nkt);
    View2D dens_can(dens_can_data.data(), im, nkt);
    View2D dkt_can(dkt_can_data.data(), im, nkt);
    View2D dku_can(dku_can_data.data(), im, nkt);

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        // High-performance CPU stack-allocated local arrays (km bounded safely at 256)
        double zcan3[3];
        double z2[256];
        double ta3[256];
        double qv3[256];
        double ws3[256];
        double dkt3[256];
        double dku3[256];
        double prsl3[256];
        double dens3[256];
        int klower_can[3];

        double hcan = cfch[i];

        // Set the initial values of the heights of the inserted canopy layers
        for (size_t kc = 0; kc < nkc; ++kc) {
            zcan3[kc] = hcan * can_frac[kc];
        }

        // Heights of original model layers
        for (size_t k = 0; k < km; ++k) {
            z2[k] = zl[i, k];
            ta3[k] = t1[i, k];
            qv3[k] = q2m[i]; // nominal water vapor specific humidity
            prsl3[k] = prsl[i, k];
            dens3[k] = dens[i, k];
            ws3[k] = std::sqrt(u1[i, k] * u1[i, k] + v1[i, k] * v1[i, k]);
            dkt3[k] = dkt[i, k];
            dku3[k] = dku[i, k];
        }
        z2[km] = 0.0;

        // Adjust canopy levels to prevent distance smaller than del (0.2m)
        for (size_t k = 0; k < km; ++k) {
            for (size_t kc = 0; kc < nkc; ++kc) {
                if (std::abs(z2[k] - zcan3[kc]) < del) {
                    double ddel = std::max(0.0, del - std::abs(zcan3[kc] - z2[k]));
                    zcan3[kc] += (zcan3[kc] > z2[k] ? ddel : -ddel);
                }
            }
        }

        // Set initial values of combined height array
        for (size_t k = 0; k < km; ++k) {
            zmid_can[i, k] = z2[k];
        }
        for (size_t kc = 0; kc < nkc; ++kc) {
            zmid_can[i, km + kc] = zcan3[kc];
        }

        // Sort combined height array monotonically decreasing
        for (size_t npass = 0; npass < nkc + 1; ++npass) {
            for (size_t k = nkt - 1; k >= 1; --k) {
                if (zmid_can[i, k] > zmid_can[i, k - 1]) {
                    double tmp = zmid_can[i, k - 1];
                    zmid_can[i, k - 1] = zmid_can[i, k];
                    zmid_can[i, k] = tmp;
                }
            }
        }

        // Determine location maps kcan3 and kmod
        for (size_t kc = 0; kc < nkc; ++kc) {
            for (size_t kk = 0; kk < nkt; ++kk) {
                if (zmid_can[i, kk] == zcan3[kc]) {
                    kcan3_data[i * nkc + kc] = kk;
                    break;
                }
            }
        }
        for (size_t k = 0; k < km; ++k) {
            for (size_t kk = 0; kk < nkt; ++kk) {
                if (zmid_can[i, kk] == z2[k]) {
                    kmod_data[i * km + k] = kk;
                    break;
                }
            }
        }

        // Create momentum height layers (zmom_can)
        for (size_t k = 0; k < km; ++k) {
            zmom_can[i, k] = zm[i, k];
        }
        for (size_t k = km; k < nkt + 1; ++k) {
            zmom_can[i, k] = (zmid_can[i, k - 1] + zmid_can[i, k]) * 0.5;
        }

        // Carry over original model values for matching layers
        for (size_t k = 0; k < km; ++k) {
            size_t kk = kmod_data[i * km + k];
            t1_can[i, kk] = ta3[k];
            qv_can[i, kk] = qv3[k];
            prsl_can[i, kk] = prsl3[k];
            dens_can[i, kk] = dens3[k];
            ws_can[i, kk] = ws3[k];
            dkt_can[i, kk] = dkt3[k];
            dku_can[i, kk] = dku3[k];
        }

        // Interpolate variables to in-canopy layers
        for (size_t kc = 0; kc < nkc; ++kc) {
            size_t kk = kcan3_data[i * nkc + kc];
            size_t k2 = km - 1; // default bounding layer
            for (size_t k = 0; k < km; ++k) {
                if (zcan3[kc] > z2[k]) {
                    k2 = k;
                    break;
                }
            }
            klower_can[kc] = k2;

            double zm2 = (zcan3[kc] - z2[k2]) / (z2[k2-1] - z2[k2]);
            t1_can[i, kk] = ta3[k2] + (ta3[k2-1] - ta3[k2]) * zm2;
            qv_can[i, kk] = qv3[k2] + (qv3[k2-1] - qv3[k2]) * zm2;
            prsl_can[i, kk] = prsl[i, 0] * (zmid_can[i, kk] / zl[i, 0]);
            dens_can[i, kk] = prsl_can[i, kk] / (rdgas * t1_can[i, kk]);
            ws_can[i, kk] = ws3[k2] + (ws3[k2-1] - ws3[k2]) * zm2;
            dkt_can[i, kk] = dkt3[k2] + (dkt3[k2-1] - dkt3[k2]) * zm2;
            dku_can[i, kk] = dku3[k2] + (dku3[k2-1] - dku3[k2]) * zm2;
        }
    }
}

/**
 * @brief Performs bidirectional, mass-conservative transport of chemical tracers and tendencies between standard model layers and combined canopy layers.
 *
 * This function translates `canopy_transfer_run` from `canopy_transfer.F90`. It enforces absolute mass conservation:
 * - **Forward transport (resolved to canopy, flag != 1)**: Maps tracer concentrations from the coarse model grid to the refined 
 *   canopy grid. Individual container masses \f$M_{res}\f$ are partitioned evenly among subgrid children:
 *   \f[
 *   M_{can}(kk) = M_{res}(k) \times \gamma_{r2c}(kk)
 *   \f]
 * - **Backward transport (canopy to resolved, flag == 1)**: Aggregates subgrid masses from canopy layers back into standard coarse layers,
 *   preserving absolute mass:
 *   \f[
 *   M_{res}(k) = \sum_{kk \in k} M_{can}(kk) \times \gamma_{c2r}(kk)
 *   \f]
 *
 * Tracer concentrations are converted internally to mass mixing ratios (MMR) scaled to micrograms per kilogram (\f$\mu g/kg\f$)
 * to protect numerical accuracy during sparse matrix sweeps, and converted back to standard units before exiting.
 *
 * @param im Number of horizontal grid columns.
 * @param km Number of vertical model layers.
 * @param nkc Number of subgrid canopy levels.
 * @param nkt Number of total combined levels (km + nkc).
 * @param ntrac Number of tracers to map.
 * @param flag Operation direction: `1` for canopy-to-resolved (backward), other values for resolved-to-canopy (forward).
 * @param zi Interface geopotential height profile (m).
 * @param zl Layer geopotential height profile (m).
 * @param zm Geopotential interface heights of model layers (m).
 * @param q1_mod Input/Output tracer concentration profile on original model grid (kg/kg).
 * @param q1_can Input/Output tracer concentration profile on expanded combined canopy grid (kg/kg).
 * @param kmod_data Input mapping indices tracking model layers inside combined array.
 * @param kcan3_data Input mapping indices tracking canopy layers inside combined array.
 * @param zmid_can_data Combined heights array at layer centers (m).
 * @param zmom_can_data Combined heights array at layer interfaces (m).
 * @param prsl_can_data Combined pressure profile at layer centers (Pa).
 * @param dens_can_data Combined air density profile (kg/m3).
 */
inline void canopy_transfer_run(
    size_t im, size_t km, size_t nkc, size_t nkt,
    size_t ntrac, int flag, ConstView2D zi, ConstView2D zl, ConstView2D zm,
    View2D q1_mod, View2D q1_can,
    const std::vector<int>& kmod_data, const std::vector<int>& kcan3_data,
    std::vector<double>& zmid_can_data, std::vector<double>& zmom_can_data,
    std::vector<double>& prsl_can_data, std::vector<double>& dens_can_data
) {
    const double reverse_conv = 1.0e9; // kg/kg -> ug/kg
    const double forward_conv = 1.0e-9; // ug/kg -> kg/kg

    std::vector<double> massair(im * km, 0.0);
    std::vector<double> massair_can(im * nkt, 0.0);
    std::vector<int> nfrct(nkt * im, 0);
    std::vector<int> ifrct(nkt * 2 * im, 0);
    std::vector<double> frctr2c(nkt * 2 * im, 1.0);
    std::vector<double> frctc2r(nkt * 2 * im, 1.0);

    View2D zmid_can(zmid_can_data.data(), im, nkt);
    View2D zmom_can(zmom_can_data.data(), im, nkt + 1);
    View2D prsl_can(prsl_can_data.data(), im, nkt);
    View2D dens_can(dens_can_data.data(), im, nkt);

    // 1. Calculate mass of air on combined and standard model levels (using absolute thickness)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        zmom_can[i, nkt] = 0.0;
        for (int k = nkt - 1; k >= 0; --k) {
            massair_can[i * nkt + k] = dens_can[i, k] * 1.6e8 * std::abs(zmom_can[i, k] - zmom_can[i, k + 1]);
        }
        for (int k = km - 1; k >= 0; --k) {
            massair[i * km + k] = dens_can[i, k] * 1.6e8 * std::abs(zm[i, k] - (k == 0 ? 0.0 : zm[i, k-1]));
        }
    }

    // 2. Compute fractional contribution arrays (nfrct, ifrct, frctr2c, frctc2r)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        // Count how many combined layers map to each resolved layer
        std::vector<size_t> resolve_counts(km, 0);
        for (size_t k = 0; k < nkt; ++k) {
            size_t kc = std::min(k, km - 1);
            resolve_counts[kc]++;
        }

        for (size_t k = 0; k < nkt; ++k) {
            size_t kc = std::min(k, km - 1);
            nfrct[i * nkt + k] = 1;
            ifrct[(i * nkt + k) * 2 + 0] = kc;
            frctr2c[(i * nkt + k) * 2 + 0] = 1.0 / resolve_counts[kc];
            frctc2r[(i * nkt + k) * 2 + 0] = 1.0; // 100% of canopy layer mass goes back to its resolved container
        }
    }

    // 3. Perform bidirectional mass-conservative tracer transport
    if (flag == 1) { // canopy_to_resolved
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            for (size_t s = 0; s < ntrac; ++s) {
                // High-performance CPU stack-allocated local arrays (nkt bounded safely at 259)
                double mmr_canopy[259];
                for (size_t k = 0; k < nkt; ++k) {
                    mmr_canopy[k] = reverse_conv * q1_can[i, k];
                }
                for (size_t k = 0; k < km; ++k) {
                    double mass_resolved = 0.0;
                    for (size_t kk = 0; kk < nkt; ++kk) {
                        if (ifrct[(i * nkt + kk) * 2 + 0] == static_cast<int>(k)) {
                            mass_resolved += mmr_canopy[kk] * massair_can[i * nkt + kk] * frctc2r[(i * nkt + kk) * 2 + 0];
                        }
                    }
                    double mmr_resolved = mass_resolved / std::max(massair[i * km + k], 1e-10);
                    q1_mod[i, k] = forward_conv * mmr_resolved;
                }
            }
        }
    } else { // resolved_to_canopy
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            for (size_t s = 0; s < ntrac; ++s) {
                // High-performance CPU stack-allocated local arrays (km bounded safely at 256)
                double mass_resolved[256];
                for (size_t k = 0; k < km; ++k) {
                    mass_resolved[k] = reverse_conv * q1_mod[i, k] * massair[i * km + k];
                }
                for (size_t k = 0; k < nkt; ++k) {
                    double mass_canopy_val = 0.0;
                    for (size_t kk = 0; kk < 1; ++kk) {
                        size_t kc = ifrct[(i * nkt + k) * 2 + kk];
                        mass_canopy_val += mass_resolved[kc] * frctr2c[(i * nkt + k) * 2 + kk];
                    }
                    double mmr_canopy_val = mass_canopy_val / std::max(massair_can[i * nkt + k], 1e-10);
                    q1_can[i, k] = forward_conv * mmr_canopy_val;
                }
            }
        }
    }
}

} // namespace canopy
} // namespace satmedmf

#endif // SATMEDMF_CANOPY_HPP
