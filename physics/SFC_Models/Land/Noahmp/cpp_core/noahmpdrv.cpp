#include "noahmpdrv.hpp"
#include "noahmp_science/common.hpp"
#include "noahmp_science/atmosphere.hpp"
#include "noahmp_science/stomata.hpp"
#include "noahmp_science/energy.hpp"
#include "noahmp_science/water.hpp"
#include "noahmp_science/radiation.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

/**
 * @file noahmpdrv.cpp
 * @brief High-performance C++23 std::mdspan implementation of the GFS Noah-MP LSM CCPP Driver core (noahmpdrv.F90).
 */

namespace noahmp {

/**
 * @brief High-performance parallel physical solver for GFS Noah-MP Land Surface Model driver.
 */
void noahmp_sflx_run(
    size_t columns, size_t soil_layers,
    double dt,
    LandSounding sounding,
    const double* sfctmp,
    const double* sfcprs,
    const double* q2,
    const double* soldn,
    const double* lwdn,
    const double* wind,
    const double* ch,
    const int* is_glacier,
    const NoahMP_Config& config,
    View1D sheat,
    View1D eta,
    View1D gflux,
    View1D runoff
) {
    // 1. Validate forcing boundaries before starting parallel execution loops (FR-006, FR-007)
    for (size_t i = 0; i < columns; ++i) {
        if (sfctmp[i] <= 0.0) {
            throw std::runtime_error("Fatal: Surface air temperature (sfctmp) is below or equal to absolute zero: " + std::to_string(sfctmp[i]) + " K.");
        }
        if (sfcprs[i] <= 0.0) {
            throw std::runtime_error("Fatal: Surface atmospheric pressure (sfcprs) is below or equal to zero: " + std::to_string(sfcprs[i]) + " Pa.");
        }
    }

    // Thread-safe parallel execution across horizontal columns
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        // Initialize outputs
        sheat[i] = 0.0;
        eta[i] = 0.0;
        gflux[i] = 0.0;
        runoff[i] = 0.0;

        // Apply automatic physical clamping constraints for grid interpolation math errors (FR-007, EC-002)
        double q2_clamped = science::clamp_to_physical_minimum(q2[i]);
        double soldn_clamped = science::clamp_to_physical_minimum(soldn[i]);
        double lwdn_clamped = science::clamp_to_physical_minimum(lwdn[i]);

        // Compute dry air density at the surface using atmosphere preprocessor
        double t_air = sfctmp[i];
        double p_air = sfcprs[i];
        double rho_air = science::compute_dry_air_density(p_air, t_air);

        // Fetch inputs from our interface (removing the hardcoded dummy values)
        double wind_speed = wind[i];
        double exchange_coeff = ch[i];

        // 2. Compute albedos using radiation sub-solver with Sellers canopy-soil blending and radiation options
        double albedo_short = 0.0;
        double albedo_long = 0.0;
        double cosz = 0.8; // simulated zenith angle
        science::compute_radiation_fluxes(config, cosz, 1.5, 0.5, albedo_short, albedo_long);

        // 3. Dynamic surface type branching / glacier bypassing (EC-003, FR-008)
        bool is_glacier_tile = (is_glacier[i] != 0); // Correctly utilize landuse classifications instead of T_ground proxy

        double r_stomata = 0.0;
        if (!is_glacier_tile) {
            // Vegetated surface: compute stomatal resistance via stomata sub-solver (with soil moisture stress control!)
            double q_sat_canopy = 3.8e-3 * std::exp(17.67 * (sounding.tv[i] - 273.15) / (sounding.tv[i] - 29.65));
            double soil_moisture_layer1 = sounding.sh2o[i, 0];
            r_stomata = science::compute_stomatal_resistance(config, t_air, q_sat_canopy, q2_clamped, soil_moisture_layer1);
        } else {
            // Non-vegetated glacier tile: completely bypass stomatal resistance/transpiration calculations
            r_stomata = 1e6; // infinite resistance representing zero stomatal transpiration
        }

        // 4. Compute Sensible Heat Flux (sheat) using energy solver
        double t_canopy = sounding.tv[i];
        sheat[i] = science::compute_sensible_heat_flux(rho_air, exchange_coeff, wind_speed, t_canopy, t_air);

        // 5. Compute Latent Heat Flux (eta) using water solver
        eta[i] = science::compute_latent_heat_flux(rho_air, exchange_coeff, wind_speed, t_canopy, q2_clamped, r_stomata);

        // 6. Compute Ground Heat Flux (gflux) using energy conduction solver
        double t_ground = sounding.tg[i];
        double t_soil1 = sounding.stc[i, 0];
        double dz1 = sounding.sldpst[i, 0];
        gflux[i] = science::compute_ground_heat_flux(t_ground, t_soil1, dz1);

        // 7. Compute surface runoff based on soil saturation excess using water solver (with dynamic options)
        double soil_saturation = sounding.smc[i, 0]; // 0-based index for first soil layer
        double precip_flux = 0.01 * soldn_clamped / 1000.0; // Simulated precip proxy scaled from radiation
        runoff[i] = science::compute_runoff(config, soil_saturation, precip_flux);
    }
}

} // namespace noahmp

extern "C" {

/**
 * @brief Flat C ABI entry point for the C++23 GFS Noah-MP LSM CCPP Driver.
 * Returns an integer error status code to safely handle exceptions across the Fortran boundary.
 */
int c_noahmp_sflx_run(
    size_t columns, size_t soil_layers,
    double dt,
    const double* stc, const double* smc, const double* sh2o, const double* sldpst,
    const double* tg, const double* tv,
    const double* sfctmp, const double* sfcprs, const double* q2,
    const double* soldn, const double* lwdn,
    const double* wind, const double* ch,
    const int* is_glacier,
    // 19 configuration parameters passed individually to keep C ABI flat
    int idveg, int iopt_crs, int iopt_btr, int iopt_run, int iopt_sfc,
    int iopt_frz, int iopt_inf, int iopt_rad, int iopt_alb, int iopt_snf,
    int iopt_tbot, int iopt_stc, int iopt_trs, int iopt_diag, int iopt_rsf,
    int iopt_soil, int iopt_pedo, int iopt_crop, int iopt_gla, int iopt_z0m,
    double* sheat, double* eta, double* gflux, double* runoff
) {
    try {
        noahmp::ConstView2D stc_view(stc, columns, soil_layers);
        noahmp::ConstView2D smc_view(smc, columns, soil_layers);
        noahmp::ConstView2D sh2o_view(sh2o, columns, soil_layers);
        noahmp::ConstView2D sldpst_view(sldpst, columns, soil_layers);
        noahmp::ConstView1D tg_view(tg, columns);
        noahmp::ConstView1D tv_view(tv, columns);

        noahmp::View1D sheat_view(sheat, columns);
        noahmp::View1D eta_view(eta, columns);
        noahmp::View1D gflux_view(gflux, columns);
        noahmp::View1D runoff_view(runoff, columns);

        noahmp::LandSounding sounding{
            stc_view,
            smc_view,
            sh2o_view,
            sldpst_view,
            tg_view,
            tv_view
        };

        // Package individual configuration flags into our config struct
        noahmp::NoahMP_Config config{
            idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc,
            iopt_frz, iopt_inf, iopt_rad, iopt_alb, iopt_snf,
            iopt_tbot, iopt_stc, iopt_trs, iopt_diag, iopt_rsf,
            iopt_soil, iopt_pedo, iopt_crop, iopt_gla, iopt_z0m
        };

        noahmp::noahmp_sflx_run(
            columns, soil_layers,
            dt,
            sounding,
            sfctmp,
            sfcprs,
            q2,
            soldn,
            lwdn,
            wind,
            ch,
            is_glacier,
            config,
            sheat_view,
            eta_view,
            gflux_view,
            runoff_view
        );
        return 0; // Success
    } catch (const std::exception& e) {
        std::cerr << "\n[C++ Boundary Error] " << e.what() << "\n" << std::endl;
        return 1; // Validation or runtime error status
    } catch (...) {
        std::cerr << "\n[C++ Boundary Error] An unknown fatal error occurred inside the C++ physical solver.\n" << std::endl;
        return 2; // Unknown fatal error status
    }
}

}
