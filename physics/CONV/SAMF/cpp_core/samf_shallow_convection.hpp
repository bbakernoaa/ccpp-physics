#ifndef SAMF_SHALLOW_CONVECTION_HPP
#define SAMF_SHALLOW_CONVECTION_HPP

#include "samf_types.hpp"
#include "samf_constants.hpp"
#include <stddef.h>

namespace samf {
namespace shallow {

/**
 * @brief Public interface for GFS SAMF Shallow Convection Solver core (samfshalcnv).
 *
 * Runs physical parameterizations of the scale-aware shallow mass-flux convection scheme.
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
 * @param dt_mf Output convective detrainment mass flux profile (kg/m2/s) [columns, layers].
 * @param kbot Output convective cloud base level indices [columns].
 * @param ktop Output convective cloud top level indices [columns].
 */
void samf_shallow_convection_run(
    size_t columns, size_t layers, double dt,
    ConvectiveSounding sounding,
    View2D dt_t, View2D dt_q, View2D dt_u, View2D dt_v,
    View2D ud_mf, View2D dt_mf,
    int* kbot, int* ktop
);

} // namespace shallow
} // namespace samf

#endif // SAMF_SHALLOW_CONVECTION_HPP
