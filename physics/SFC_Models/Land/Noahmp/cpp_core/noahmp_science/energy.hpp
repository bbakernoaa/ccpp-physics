#ifndef NOAHMP_SCIENCE_ENERGY_HPP
#define NOAHMP_SCIENCE_ENERGY_HPP

#include "../noahmp_types.hpp"
#include "../noahmp_constants.hpp"
#include "common.hpp"
#include <cmath>

namespace noahmp {
namespace science {

/**
 * @brief Computes sensible heat flux using temperature gradients between surface canopy and ground layers.
 */
inline Real compute_sensible_heat_flux(
    Real rho_air, Real ch, Real wind_speed, Real t_canopy, Real t_air
) {
    return rho_air * constants::cp * ch * wind_speed * (t_canopy - t_air);
}

/**
 * @brief Computes ground conduction heat flux into first soil layer.
 */
inline Real compute_ground_heat_flux(
    Real t_ground, Real t_soil1, Real dz1, Real k_soil = 1.2
) {
    return k_soil * (t_ground - t_soil1) / std::max(0.01, dz1);
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCE_ENERGY_HPP
