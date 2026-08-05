#ifndef GWD_CONSTANTS_HPP
#define GWD_CONSTANTS_HPP

namespace gwd {
namespace constants {

// Standard GFS/UGWP-v1 physical constants
constexpr double grav = 9.80665;             // Acceleration due to gravity (m/s2)
constexpr double pi = 3.141592653589793;     // Mathematical constant Pi
constexpr double rd = 287.05;                // Gas constant for dry air (J/kg/K)
constexpr double cp = 1004.6;                // Specific heat of dry air (J/kg/K)
constexpr double rv = 461.5;                 // Gas constant for water vapor (J/kg/K)
constexpr double eps = rd / rv;              // Ratio of gas constants (0.622)
constexpr double epsm1 = eps - 1.0;          // eps - 1.0 (-0.378)
constexpr double rerth = 6.3712e6;           // Mean radius of Earth (m)
constexpr double omega = 7.2921e-5;          // Angular velocity of Earth rotation (rad/s)

} // namespace constants
} // namespace gwd

#endif // GWD_CONSTANTS_HPP
