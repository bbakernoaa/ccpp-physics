#ifndef NOAHMP_SCIENCE_ATMOSPHERE_HPP
#define NOAHMP_SCIENCE_ATMOSPHERE_HPP

#include "../noahmp_types.hpp"
#include "../noahmp_constants.hpp"
#include "common.hpp"
#include <cmath>

namespace noahmp {
namespace science {

/**
 * @brief Computes surface boundary-layer scaling height and air density at reference level.
 */
inline Real compute_dry_air_density(Real p_air, Real t_air) {
    return p_air / (constants::rd * t_air);
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCE_ATMOSPHERE_HPP
