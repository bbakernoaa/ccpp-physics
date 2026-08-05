#ifndef PHOTOCHEM_CONSTANTS_HPP
#define PHOTOCHEM_CONSTANTS_HPP

#include "photochem_types.hpp"

namespace photochem {
namespace constants {

// Standard GFS photochemistry and thermodynamic constants
constexpr Real grav = 9.80665;             // Acceleration due to gravity (m/s2)
constexpr Real r_methane_co2 = 2.0;        // Methane oxidation H2O source ratio multiplier

} // namespace constants
} // namespace photochem

#endif // PHOTOCHEM_CONSTANTS_HPP
