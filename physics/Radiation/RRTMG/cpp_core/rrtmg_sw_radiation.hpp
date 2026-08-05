#ifndef RRTMG_SW_RADIATION_HPP
#define RRTMG_SW_RADIATION_HPP

#include "rrtmg_types.hpp"
#include "rrtmg_constants.hpp"
#include <stddef.h>

namespace rrtmg {
namespace sw {

/**
 * @brief Public interface for GFS RRTMG Shortwave Radiation Solver (rrtmg_sw).
 *
 * Runs the physical parameterizations of the shortwave solar radiation absorption,
 * scattering, and reflection loops.
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param bands Number of shortwave solar spectral bands (standard: 14 bands).
 * @param sounding Input vertical radiative atmospheric profile state.
 * @param cos_solar_zenith Input cosine of solar zenith angle values [columns].
 * @param sw_heating_rate Output computed shortwave warming rates profile (K/day) [columns, layers].
 * @param sw_flux_down Output computed downward shortwave flux profile (W/m2) [columns, layers + 1].
 * @param sw_flux_up Output computed upward shortwave flux profile (W/m2) [columns, layers + 1].
 */
void rrtmg_sw_radiation_run(
    size_t columns, size_t layers, size_t bands,
    RadiativeSounding sounding,
    const Real* cos_solar_zenith,
    View2D sw_heating_rate,
    View2D sw_flux_down,
    View2D sw_flux_up
);

} // namespace sw
} // namespace rrtmg

#endif // RRTMG_SW_RADIATION_HPP
