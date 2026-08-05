#ifndef UGWPV1_GSLDRAG_HPP
#define UGWPV1_GSLDRAG_HPP

#include "gwd_types.hpp"
#include "gwd_constants.hpp"
#include <stddef.h>

namespace gwd {
namespace solver {

/**
 * @brief Public C++23 interface for UGWP v1 GSL gravity wave drag solver (ugwpv1_gsldrag).
 *
 * Runs the physical parameterizations of the Unified Gravity Wave Physics (UGWP) v1 solvers,
 * including subgrid orographic gravity wave drag (OGW), turbulent orographic form drag (TOFD),
 * and non-orographic gravity wave drag (NGW).
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param dtp Physics time step (s).
 * @param sounding Active physical atmospheric profiles.
 * @param do_gsl_drag_ls_bl Active GSL large-scale boundary-layer drag flag.
 * @param do_gsl_drag_ss Active GSL small-scale gravity wave drag flag.
 * @param do_gsl_drag_tofd Active Turbulent Orographic Form Drag (TOFD) flag.
 * @param do_ngw_ec Active Non-Orographic Gravity Wave (NGW) drag flag.
 * @param hprime topographical standard deviation profile.
 * @param oc Topographical convexity profile.
 * @param theta Topographical mountain orientation profile.
 * @param sigma Topographical mountain slope profile.
 * @param gamma Topographical mountain anisotropy profile.
 * @param elvmax Maximum elevation profile.
 * @param clx Fraction of mountain subgrid cover.
 * @param oa4 Topographical multi-directional orientation profiles.
 * @param varss Subgrid mountain elevation variance.
 * @param dx Grid horizontal spacing (m).
 * @param xlat Grid columns latitudes (rad).
 * @param area Grid columns area values (m2).
 * @param dudt_ogw Output computed orographic wind zonal tendency (m/s/s).
 * @param dvdt_ogw Output computed orographic wind meridional tendency (m/s/s).
 * @param dudt_ngw Output computed non-orographic wind zonal tendency (m/s/s).
 * @param dvdt_ngw Output computed non-orographic wind meridional tendency (m/s/s).
 * @param dtdt_ngw Output computed non-orographic temperature tendency (K/s).
 * @param dudt_ofd Output computed TOFD wind zonal tendency (m/s/s).
 * @param dvdt_ofd Output computed TOFD wind meridional tendency (m/s/s).
 * @param tau_ogw Output orographic gravity wave stress profile (Pa).
 * @param tau_ngw Output non-orographic gravity wave stress profile (Pa).
 */
void ugwpv1_gsldrag_run(
    size_t columns, size_t layers, double dtp,
    GwdSounding sounding,
    bool do_gsl_drag_ls_bl, bool do_gsl_drag_ss, bool do_gsl_drag_tofd, bool do_ngw_ec,
    const double* hprime, const double* oc, const double* theta, const double* sigma, const double* gamma,
    const double* elvmax, const double* clx, const double* oa4, const double* varss,
    const double* dx, const double* xlat, const double* area,
    View2D dudt_ogw, View2D dvdt_ogw,
    View2D dudt_ngw, View2D dvdt_ngw, View2D dtdt_ngw,
    View2D dudt_ofd, View2D dvdt_ofd,
    View2D tau_ogw, View2D tau_ngw
);

} // namespace solver
} // namespace gwd

#endif // UGWPV1_GSLDRAG_HPP
