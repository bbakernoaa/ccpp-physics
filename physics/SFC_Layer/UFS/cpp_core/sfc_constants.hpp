#ifndef SFC_CONSTANTS_HPP
#define SFC_CONSTANTS_HPP

#include "sfc_types.hpp"

namespace sfc {
namespace constants {

// Standard GFS physical constants for surface layer calculations
constexpr Real grav = 9.80665;             // Acceleration due to gravity (m/s2)
constexpr Real cp = 1004.6;                // Specific heat of dry air (J/kg/K)
constexpr Real rd = 287.05;                // Gas constant for dry air (J/kg/K)
constexpr Real karman = 0.4;               // Von Karman constant

} // namespace constants
} // namespace sfc

#endif // SFC_CONSTANTS_HPP
