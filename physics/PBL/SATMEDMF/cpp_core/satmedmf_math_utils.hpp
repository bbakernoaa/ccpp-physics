#ifndef SATMEDMF_MATH_UTILS_HPP
#define SATMEDMF_MATH_UTILS_HPP

#include <cmath>

namespace satmedmf {

// Physical Constants
constexpr double P0 = 100000.0;     // Reference pressure (Pa)
constexpr double RD = 287.05;       // Gas constant for dry air (J/kg/K)
constexpr double CP = 1004.6;       // Specific heat at constant pressure (J/kg/K)
constexpr double RD_OVER_CP = RD / CP; // Poisson constant (approx 0.2857)

/**
 * Computes the dimensionless Exner function value: PI = (P / P0)^(R_d / C_p)
 *
 * @param pressure Pressure (Pa)
 * @return Exner function value
 */
inline double compute_exner(double pressure) {
    if (pressure <= 0.0) return 0.0;
    return std::pow(pressure / P0, RD_OVER_CP);
}

/**
 * Computes the saturation vapor pressure using standard Tetens formula (over liquid water).
 *
 * @param temp Temperature (K)
 * @return Saturation vapor pressure (Pa)
 */
inline double compute_sat_vapor_pressure(double temp) {
    const double tc = temp - 273.15; // Convert to Celsius
    if (tc > 0.0) {
        // Over liquid water
        return 611.2 * std::exp((17.67 * tc) / (tc + 243.5));
    } else {
        // Over ice
        return 611.2 * std::exp((21.874 * tc) / (tc + 265.5));
    }
}

} // namespace satmedmf

#endif // SATMEDMF_MATH_UTILS_HPP
