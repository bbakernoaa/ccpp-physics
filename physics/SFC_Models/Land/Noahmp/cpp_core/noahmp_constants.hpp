#ifndef NOAHMP_CONSTANTS_HPP
#define NOAHMP_CONSTANTS_HPP

#include "noahmp_types.hpp"

namespace noahmp {
namespace constants {

// Standard GFS/Noah-MP physical constants
constexpr Real grav = 9.80665;             // Acceleration due to gravity (m/s2)
constexpr Real cp = 1004.6;                // Specific heat of dry air (J/kg/K)
constexpr Real rd = 287.05;                // Gas constant for dry air (J/kg/K)
constexpr Real sbc = 5.670367e-8;          // Stefan-Boltzmann constant (W/m2/K4)

} // namespace constants
} // namespace noahmp

#endif // NOAHMP_CONSTANTS_HPP
