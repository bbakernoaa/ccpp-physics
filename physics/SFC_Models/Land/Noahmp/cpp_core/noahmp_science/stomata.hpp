#ifndef NOAHMP_SCIENCESTOMATA_HPP
#define NOAHMP_SCIENCESTOMATA_HPP

#include "../noahmp_types.hpp"
#include "../noahmp_constants.hpp"
#include "common.hpp"
#include <cmath>
#include <algorithm>

namespace noahmp {
namespace science {

/**
 * @brief Computes stomatal resistance adjusted for vapor pressure deficit, temperature, and soil water stress.
 * Dynamically branches based on iopt_crs (resistance scheme) and iopt_btr (soil moisture stress scheme).
 */
inline Real compute_stomatal_resistance(
    const NoahMP_Config& config,
    Real t_air, Real q_sat, Real q2, Real soil_moisture_layer1, Real r_min = 40.0
) {
    // 1. Establish soil moisture stress coefficient (beta) based on iopt_btr
    Real beta = 1.0;
    const Real theta_wilt = 0.10;
    const Real theta_fc   = 0.35;

    if (config.iopt_btr == 2) {
        // Option 2: Logarithmic CLM-type soil moisture stress
        Real fraction = (soil_moisture_layer1 - theta_wilt) / (theta_fc - theta_wilt);
        fraction = std::max(0.01, std::min(1.0, fraction));
        beta = std::log(1.0 + 9.0 * fraction) / std::log(10.0);
    } else if (config.iopt_btr == 3) {
        // Option 3: Quadratic SSiB-type soil moisture stress
        Real fraction = (soil_moisture_layer1 - theta_wilt) / (theta_fc - theta_wilt);
        fraction = std::max(0.01, std::min(1.0, fraction));
        beta = std::pow(fraction, 2);
    } else {
        // Option 1 (Default): Linear Noah-type soil moisture stress
        beta = (soil_moisture_layer1 - theta_wilt) / (theta_fc - theta_wilt);
    }
    beta = std::max(0.01, std::min(1.0, beta)); // Clamp stress coefficient to [0.01, 1.0]

    // 2. Establish basic environmental scaling factors
    Real deficit = std::max(0.0, q_sat - q2);
    Real f_vpd = 1.0 / (1.0 + 4000.0 * deficit); // vapor pressure deficit factor
    Real f_temp = 1.0 - 0.0016 * std::pow(298.15 - t_air, 2); // temperature factor
    f_temp = std::max(0.01, f_temp);

    // 3. Select stomatal resistance scheme based on iopt_crs
    if (config.iopt_crs == 1) {
        // Option 1: Ball-Berry photosynthesis-driven stomatal conductance
        // Conducting resistance is inversely proportional to net carbon assimilation rate (simulated)
        Real net_assimilation_proxy = 15.0 * beta * f_temp; // mu_mol / m2 / s
        Real bb_conductance = 0.01 + 9.0 * net_assimilation_proxy / std::max(0.01, 1.0 + 1000.0 * deficit);
        return 1.0 / std::max(0.001, bb_conductance);
    } else {
        // Option 2 (Default): Jarvis-type empirical stomatal resistance
        return r_min / (beta * f_vpd * f_temp);
    }
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCESTOMATA_HPP
