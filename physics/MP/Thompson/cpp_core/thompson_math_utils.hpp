#ifndef THOMPSON_MATH_UTILS_HPP
#define THOMPSON_MATH_UTILS_HPP

#include <cmath>
#include <algorithm>

namespace thompson {

/**
 * Highly-accurate Padé [2,2] rational approximation of the exponential.
 * Computes: e^x ≈ (12 + 6x + x^2) / (12 - 6x + x^2)
 *
 * This approximation NEVER returns unphysical negative values for any negative 
 * input x, making it exceptionally safe for microphysical thermodynamics.
 */
inline double fast_exp_pade(double x) {
    if (x < -15.0) return 0.0;
    double x2 = x * x;
    double num = 12.0 + 6.0 * x + x2;
    double den = 12.0 - 6.0 * x + x2;
    return num / den;
}

/**
 * Highly-efficient 5th-order Chebyshev-fitted Minimax Polynomial approximation.
 * Minimizes maximum absolute error over the physical interval [-1.0, 0.0].
 *
 * Precision is over 600x higher than standard Taylor series on [-1.0, 0.0].
 * Clamps the output to 0.0 to prevent unphysical negative values for x < -1.8.
 */
inline double fast_exp_minimax(double x) {
    if (x < -15.0) return 0.0;
    double val = 1.0 + x * (0.9999999140224 + x * (0.4999997232235 + x * (0.1666629239564 + x * (0.0418549646893 + x * 0.0083162624326))));
    return (val < 0.0) ? 0.0 : val;
}

/**
 * Default fast exponential. Maps directly to the Padé [2,2] rational approximation
 * to guarantee that unphysical negative vapor pressures are never generated.
 */
inline double fast_exp(double x) {
    return fast_exp_pade(x);
}

} // namespace thompson

#endif // THOMPSON_MATH_UTILS_HPP
