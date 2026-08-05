#ifndef NOAHMP_SCIENCE_RADIATION_HPP
#define NOAHMP_SCIENCE_RADIATION_HPP

#include "../noahmp_types.hpp"
#include "../noahmp_constants.hpp"
#include "common.hpp"
#include <cmath>
#include <algorithm>

namespace noahmp {
namespace science {

/**
 * @brief Computes dynamic canopy-soil blended albedos based on leaf area index (LAI) and selected options.
 * Matches CCPP option options (iopt_rad and iopt_alb).
 */
inline void compute_radiation_fluxes(
    const NoahMP_Config& config,
    Real cosz, Real elai, Real esai,
    Real& albedo_short, Real& albedo_long
) {
    // 1. Establish albedo profiles based on snow surface albedo option (iopt_alb)
    Real soil_albedo_short = 0.15;
    Real soil_albedo_long  = 0.20;
    Real canopy_albedo_short = 0.22;
    Real canopy_albedo_long  = 0.26;

    if (config.iopt_alb == 2) {
        // CLASS-type snow albedo parameters
        soil_albedo_short = 0.17;
        soil_albedo_long  = 0.22;
    }

    // 2. Branch calculations based on radiation transfer option (iopt_rad)
    if (config.iopt_rad == 2) {
        // Option 2: Simplified, linear BATS-type shielding fraction
        Real total_lai = elai + esai;
        Real shielding = std::min(1.0, total_lai / 2.0); // simple linear fraction proxy
        albedo_short = (1.0 - shielding) * soil_albedo_short + shielding * canopy_albedo_short;
        albedo_long  = (1.0 - shielding) * soil_albedo_long  + shielding * canopy_albedo_long;
    } else {
        // Option 1 (Default): Two-stream Sellers exponential shielding
        Real total_lai = elai + esai;
        Real shielding = 1.0 - std::exp(-0.5 * total_lai / std::max(0.01, cosz));
        albedo_short = (1.0 - shielding) * soil_albedo_short + shielding * canopy_albedo_short;
        albedo_long  = (1.0 - shielding) * soil_albedo_long  + shielding * canopy_albedo_long;
    }
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCE_RADIATION_HPP
