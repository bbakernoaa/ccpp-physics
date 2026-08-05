#ifndef NOAHMP_SCIENCE_COMMON_HPP
#define NOAHMP_SCIENCE_COMMON_HPP

#include "../noahmp_types.hpp"
#include <algorithm>

namespace noahmp {
namespace science {

/**
 * @brief Clamps a near-zero negative atmospheric forcing parameter to non-negative limits (FR-007, EC-002).
 */
inline Real clamp_to_physical_minimum(Real val, Real min_val = 0.0) {
    return std::max(min_val, val);
}

} // namespace science
} // namespace noahmp

#endif // NOAHMP_SCIENCE_COMMON_HPP
