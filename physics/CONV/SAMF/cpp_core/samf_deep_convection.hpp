#ifndef SAMF_DEEP_CONVECTION_HPP
#define SAMF_DEEP_CONVECTION_HPP

#include "samf_types.hpp"
#include "samf_constants.hpp"
#include <stddef.h>

namespace samf {
namespace deep {

/**
 * @brief Public interface for GFS SAMF Deep Convection Solver core (samfdeepcnv).
 *
 * Runs physical parameterizations of the scale-aware deep mass-flux convection scheme.
 * Handles zero heap allocations internally within threads.
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param dt Physics time step (s).
 * @param sounding Active physical vertical sounding profiles.
 * @param dt_t Output computed temperature tendency (K/s) [columns, layers].
 * @param dt_q Output computed specific humidity tendency (kg/kg/s) [columns, layers].
 * @param dt_u Output computed zonal wind tendency (m/s/s) [columns, layers].
 * @param dt_v Output computed meridional wind tendency (m/s/s) [columns, layers].
 * @param ud_mf Output convective updraft mass flux profile (kg/m2/s) [columns, layers].
 * @param dd_mf Output convective downdraft mass flux profile (kg/m2/s) [columns, layers].
 * @param dt_mf Output convective detrainment mass flux profile (kg/m2/s) [columns, layers].
 * @param cnvw Output convective cloud liquid water profile (kg/kg) [columns, layers].
 * @param cnvc Output convective cloud fraction profile (0-1) [columns, layers].
 * @param kbot Output convective cloud base level indices [columns].
 * @param ktop Output convective cloud top level indices [columns].
 * @param kcnv Output convective deep trigger flag indices [columns].
 * @param rain Output convective surface rain rate (mm/s) [columns].
 */
void samf_deep_convection_run(
    size_t columns, size_t layers, double dt,
    ConvectiveSounding sounding,
    View2D dt_t, View2D dt_q, View2D dt_u, View2D dt_v,
    View2D ud_mf, View2D dd_mf, View2D dt_mf,
    View2D cnvw, View2D cnvc,
    int* kbot, int* ktop, int* kcnv,
    double* rain
);

} // namespace deep
} // namespace samf

#endif // SAMF_DEEP_CONVECTION_HPP
