#ifndef NOAHMPDRV_HPP
#define NOAHMPDRV_HPP

#include "noahmp_types.hpp"
#include "noahmp_constants.hpp"
#include <stddef.h>

namespace noahmp {

/**
 * @brief Public C++23 interface for GFS Noah-MP LSM CCPP Driver.
 *
 * @param columns Number of horizontal grid columns.
 * @param soil_layers Number of vertical soil layers (standard: 4 layers).
 * @param dt Physics timestep (s).
 * @param sounding Input vertical land/soil profiles state.
 * @param sfctmp Input surface air temperature profile [columns].
 * @param sfcprs Input surface air pressure profile [columns].
 * @param q2 Input specific humidity at 2m height profile [columns].
 * @param soldn Input downward solar radiation flux profile [columns].
 * @param lwdn Input downward longwave thermal radiation flux profile [columns].
 * @param wind Input wind speed profile [columns].
 * @param ch Input boundary-layer turbulent heat exchange coefficient profile [columns].
 * @param is_glacier Input surface classification identifying glacier/ice tiles [columns].
 * @param config Input physical science configuration options (all 19 options).
 * @param sheat Output computed sensible heat flux profile [columns].
 * @param eta Output computed latent heat flux profile [columns].
 * @param gflux Output computed ground heat flux profile [columns].
 * @param runoff Output computed surface water runoff profile [columns].
 */
void noahmp_sflx_run(
    size_t columns, size_t soil_layers,
    double dt,
    LandSounding sounding,
    const double* sfctmp,
    const double* sfcprs,
    const double* q2,
    const double* soldn,
    const double* lwdn,
    const double* wind,
    const double* ch,
    const int* is_glacier,
    const NoahMP_Config& config,
    View1D sheat,
    View1D eta,
    View1D gflux,
    View1D runoff
);

} // namespace noahmp

#endif // NOAHMPDRV_HPP
