#ifndef RRTMG_CONSTANTS_HPP
#define RRTMG_CONSTANTS_HPP

#include "rrtmg_types.hpp"

namespace rrtmg {
namespace constants {

// Standard GFS/RRTMG radiative and thermodynamic constants
constexpr Real grav = 9.80665;             // Acceleration due to gravity (m/s2)
constexpr Real cp = 1004.6;                // Specific heat of dry air (J/kg/K)
constexpr Real sbc = 5.670367e-8;          // Stefan-Boltzmann constant (W/m2/K4)
constexpr Real solcon = 1361.0;            // Nominal solar constant (W/m2)
constexpr Real pi = 3.141592653589793;     // Pi

} // namespace constants
} // namespace rrtmg

#endif // RRTMG_CONSTANTS_HPP
