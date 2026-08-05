#ifndef SAMF_CONSTANTS_HPP
#define SAMF_CONSTANTS_HPP

namespace samf {
namespace constants {

// Standard GFS/CCPP physical constants
constexpr double grav = 9.80665;          // Acceleration due to gravity (m/s2)
constexpr double pi = 3.141592653589793;  // Mathematical constant Pi
constexpr double rd = 287.05;             // Gas constant for dry air (J/kg/K)
constexpr double cp = 1004.6;             // Specific heat of dry air at constant pressure (J/kg/K)
constexpr double rv = 461.5;              // Gas constant for water vapor (J/kg/K)
constexpr double hvap = 2.5e6;            // Latent heat of vaporization (J/kg)
constexpr double hfus = 3.33e5;           // Latent heat of fusion (J/kg)
constexpr double eps = rd / rv;           // Ratio of gas constants (0.622)
constexpr double epsm1 = eps - 1.0;       // eps - 1.0 (-0.378)
constexpr double rhowater = 1000.0;       // Density of liquid water (kg/m3)

} // namespace constants
} // namespace samf

#endif // SAMF_CONSTANTS_HPP
