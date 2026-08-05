#ifndef NOAHMP_SCIENCE_WATER_HPP
#define NOAHMP_SCIENCE_WATER_HPP

#include "../noahmp_types.hpp"
#include "../noahmp_constants.hpp"
#include "common.hpp"
#include <cmath>
#include <algorithm>

namespace noahmp {
namespace science {

/**
 * @brief Computes latent heat flux/evapotranspiration and Clausius-Clapeyron saturation deficit.
 */
inline Real compute_latent_heat_flux(
    Real rho_air, Real ch, Real wind_speed, Real t_canopy, Real q2_clamped, Real r_stomata
) {
    Real l_vapor = 2.501e6; // J/kg
    Real q_sat_canopy = 3.8e-3 * std::exp(17.67 * (t_canopy - 273.15) / (t_canopy - 29.65));
    Real evap_rate = rho_air * ch * wind_speed * std::max(0.0, q_sat_canopy - q2_clamped) / (1.0 + ch * wind_speed * r_stomata);
    return l_vapor * evap_rate;
}

/**
 * @brief Computes surface water runoff rate based on soil saturation excess and runoff options (iopt_run).
 */
inline Real compute_runoff(const NoahMP_Config& config, Real soil_saturation, Real precip_flux) {
    if (config.iopt_run == 1) {
        // Option 1: Simple TOPMODEL (SIMGM) runoff: exponential saturation drainage
        return precip_flux * std::exp(-2.5 * (1.0 - soil_saturation));
    } else if (config.iopt_run == 2) {
        // Option 2: Schaake-type saturation excess runoff
        return precip_flux * std::pow(soil_saturation, 4);
    } else if (config.iopt_run == 3) {
        // Option 3: BATS-type runoff
        Real factor = std::max(0.0, soil_saturation - 0.15) / 0.85;
        return precip_flux * std::pow(factor, 3);
    } else {
        // Default Option: Baseline power-law runoff
        return precip_flux * std::pow(soil_saturation, 4);
    }
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCE_WATER_HPP
