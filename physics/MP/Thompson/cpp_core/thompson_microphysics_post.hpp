#ifndef THOMPSON_MICROPHYSICS_POST_HPP
#define THOMPSON_MICROPHYSICS_POST_HPP

#include <stddef.h>
#include <mdspan>
#include <algorithm>
#include <string>

/**
 * @file thompson_microphysics_post.hpp
 * @brief C++23 translation of the Thompson microphysics postprocessing scheme (mp_thompson_post.F90).
 *
 * Implements potential temperature tendency limiters protecting thermodynamic matrices 
 * from numerical stiffness or non-physical microphysical changes at high timesteps.
 */

namespace thompson {
namespace post {

// Translation-unit initialization states
static bool is_initialized = false;
static bool apply_limiter = false;

/**
 * @brief Initializes the Thompson postprocessing limiter states.
 *
 * Matches the logic of mp_thompson_post_init.
 *
 * @param ttendlim Maximum permitted potential temperature tendency change threshold. If negative, the limiter is bypassed.
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void thompson_post_init(double ttendlim, std::string& errmsg, int& errflg) {
    errmsg = "";
    errflg = 0;

    if (is_initialized) return;

    if (ttendlim < 0.0) {
        apply_limiter = false;
    } else {
        apply_limiter = true;
    }

    is_initialized = true;
}

/**
 * @brief Executes potential temperature tendency limiting for Thompson microphysics.
 *
 * This routine translates `mp_thompson_post_run`. If `apply_limiter` is active, it enforces:
 * \f[
 * \left(\frac{\partial \theta}{\partial dt}\right)_{limited} = \max\left(-T_{lim}, \min\left(T_{lim}, \frac{\partial T}{\partial dt} \times \Pi^{-1}\right)\right)
 * \f]
 * and maps the limited temperature tendency back to physical space:
 * \f[
 * \left(\frac{\partial T}{\partial dt}\right)_{limited} = \left(\frac{\partial \theta}{\partial dt}\right)_{limited} \times \Pi
 * \f]
 * where \f$\Pi\f$ represents the Exner function profile (`prslk`).
 *
 * @param ncol Number of horizontal grid columns.
 * @param nlev Number of vertical model layers.
 * @param dtgrs Input/Output temperature tendency profile (K/s) [columns, layers].
 * @param tgrs Input current temperature profile (K) [columns, layers].
 * @param prslk Input Exner function profile [columns, layers].
 * @param dtp Physics timestep (s).
 * @param ttendlim Maximum potential temperature tendency threshold (K/s).
 * @param kdt Physics timestep counter.
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void thompson_post_run(
    size_t ncol, size_t nlev,
    std::mdspan<double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left> dtgrs,
    std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left> tgrs,
    std::mdspan<const double, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_left> prslk,
    double dtp, double ttendlim, int kdt,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;

    if (!is_initialized) {
        errmsg = "thompson_post_run called before thompson_post_init";
        errflg = 1;
        return;
    }

    if (!apply_limiter) return;

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < nlev; ++k) {
        for (size_t i = 0; i < ncol; ++i) {
            double mp_tend = dtgrs[i, k] / prslk[i, k];
            mp_tend = std::max(-ttendlim, std::min(ttendlim, mp_tend));
            dtgrs[i, k] = mp_tend * prslk[i, k];
        }
    }
}

/**
 * @brief Finalizes the Thompson postprocessing state.
 *
 * Matches the logic of mp_thompson_post_finalize.
 *
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void thompson_post_finalize(std::string& errmsg, int& errflg) {
    errmsg = "";
    errflg = 0;

    if (!is_initialized) return;

    is_initialized = false;
    apply_limiter = false;
}

} // namespace post
} // namespace thompson

#endif // THOMPSON_MICROPHYSICS_POST_HPP
