#ifndef UGWPV1_GSLDRAG_POST_HPP
#define UGWPV1_GSLDRAG_POST_HPP

#include "gwd_types.hpp"
#include <cmath>
#include <algorithm>
#include <string>

/**
 * @file ugwpv1_gsldrag_post.hpp
 * @brief High-performance C++23 std::mdspan implementation of the GFS GWD postprocessing diagnostics (ugwpv1_gsldrag_post.F90).
 */

namespace gwd {
namespace post {

/**
 * @brief Computes time-averaged GWD diagnostics (CCPP-compliant ugwpv1_gsldrag_post_run).
 *
 * Accumulates zonal mountain stress, gravity wave stress, and wind tendencies over physics steps.
 *
 * @param im Number of horizontal grid columns.
 * @param levs Number of vertical model layers.
 * @param ldiag_ugwp Active GWD diagnostics flag.
 * @param dtf Fractional timestep scaling factor (s).
 * @param zobl Subgrid boundary-layer diagnostic profile [columns].
 * @param zlwb Surface wave-drag breaking diagnostics [columns].
 * @param zogw Orographic surface gravity wave-drag [columns].
 * @param tau_ogw Orographic surface stress values [columns].
 * @param tau_ngw Non-orographic surface stress values [columns].
 * @param du_ofdcol-du_oblcol Diagnostic wind drag components [columns].
 * @param tot_mtb-tot_ngw Accumulated output surface diagnostics.
 * @param tot_zmtb-tot_zogw Accumulated output geopotential/PBL height diagnostics.
 * @param dudt_gw-dvdt_gw Convective wind tendencies [columns, layers].
 * @param dudt_obl-dudt_ogw Orographic and boundary-layer diagnostic wind tendencies.
 * @param du3dt_mtb-dv3dt_ngw Output accumulated 3D diagnostic wind tendencies.
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void ugwpv1_gsldrag_post_run(
    size_t im, size_t levs, bool ldiag_ugwp, double dtf,
    const double* zobl, const double* zlwb, const double* zogw,
    const double* tau_ogw, const double* tau_ngw,
    const double* du_ofdcol, const double* du_oblcol,
    double* tot_mtb, double* tot_ogw, double* tot_tofd, double* tot_ngw,
    double* tot_zmtb, double* tot_zlwb, double* tot_zogw,
    ConstView2D dudt_gw, ConstView2D dvdt_gw,
    ConstView2D dudt_obl, ConstView2D dudt_ofd, ConstView2D dudt_ogw,
    View2D du3dt_mtb, View2D du3dt_tms, View2D du3dt_ogw,
    View2D du3dt_ngw, View2D dv3dt_ngw,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;

    if (!ldiag_ugwp) return;

    // Standard CPU parallelization over independent columns
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        tot_zmtb[i] += dtf * zobl[i];
        tot_zlwb[i] += dtf * zlwb[i];
        tot_zogw[i] += dtf * zogw[i];

        tot_tofd[i] += dtf * du_ofdcol[i];
        tot_mtb[i]  += dtf * du_oblcol[i];
        tot_ogw[i]  += dtf * tau_ogw[i];
        tot_ngw[i]  += dtf * tau_ngw[i];
    }

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < levs; ++k) {
        for (size_t i = 0; i < im; ++i) {
            du3dt_mtb[i, k] += dtf * dudt_obl[i, k];
            du3dt_tms[i, k] += dtf * dudt_ofd[i, k];
            du3dt_ogw[i, k] += dtf * dudt_ogw[i, k];
            du3dt_ngw[i, k] += dtf * dudt_gw[i, k];
            dv3dt_ngw[i, k] += dtf * dvdt_gw[i, k];
        }
    }
}

} // namespace post
} // namespace gwd

extern "C" {

/**
 * @brief Flat C ABI entry point for GFS UGWP v1 Gravity Wave Drag Postprocessing diagnostics.
 */
void c_ugwpv1_gsldrag_post_run(
    size_t im, size_t levs, int ldiag_ugwp, double dtf,
    const double* zobl, const double* zlwb, const double* zogw,
    const double* tau_ogw, const double* tau_ngw,
    const double* du_ofdcol, const double* du_oblcol,
    double* tot_mtb, double* tot_ogw, double* tot_tofd, double* tot_ngw,
    double* tot_zmtb, double* tot_zlwb, double* tot_zogw,
    const double* dudt_gw, const double* dvdt_gw,
    const double* dudt_obl, const double* dudt_ofd, const double* dudt_ogw,
    double* du3dt_mtb, double* du3dt_tms, double* du3dt_ogw,
    double* du3dt_ngw, double* dv3dt_ngw
) {
    gwd::ConstView2D dudt_gw_view(dudt_gw, im, levs);
    gwd::ConstView2D dvdt_gw_view(dvdt_gw, im, levs);
    gwd::ConstView2D dudt_obl_view(dudt_obl, im, levs);
    gwd::ConstView2D dudt_ofd_view(dudt_ofd, im, levs);
    gwd::ConstView2D dudt_ogw_view(dudt_ogw, im, levs);

    gwd::View2D du3dt_mtb_view(du3dt_mtb, im, levs);
    gwd::View2D du3dt_tms_view(du3dt_tms, im, levs);
    gwd::View2D du3dt_ogw_view(du3dt_ogw, im, levs);
    gwd::View2D du3dt_ngw_view(du3dt_ngw, im, levs);
    gwd::View2D dv3dt_ngw_view(dv3dt_ngw, im, levs);

    std::string errmsg;
    int errflg = 0;

    gwd::post::ugwpv1_gsldrag_post_run(
        im, levs, ldiag_ugwp != 0, dtf,
        zobl, zlwb, zogw, tau_ogw, tau_ngw,
        du_ofdcol, du_oblcol,
        tot_mtb, tot_ogw, tot_tofd, tot_ngw,
        tot_zmtb, tot_zlwb, tot_zogw,
        dudt_gw_view, dvdt_gw_view,
        dudt_obl_view, dudt_ofd_view, dudt_ogw_view,
        du3dt_mtb_view, du3dt_tms_view, du3dt_ogw_view,
        du3dt_ngw_view, dv3dt_ngw_view,
        errmsg, errflg
    );
}

}
#endif // UGWPV1_GSLDRAG_POST_HPP
