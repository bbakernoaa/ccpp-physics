#ifndef RRTMG_LW_RADIATION_HPP
#define RRTMG_LW_RADIATION_HPP

#include "rrtmg_types.hpp"
#include "rrtmg_constants.hpp"
#include <stddef.h>

namespace rrtmg {
namespace lw {

/**
 * @brief Public interface for GFS RRTMG Longwave Radiation Solver (rrtmg_lw).
 *
 * Runs the physical parameterizations of the terrestrial longwave thermal emission,
 * atmospheric absorption, and radiative cooling loops.
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param bands Number of longwave spectral bands (standard: 16 bands).
 * @param sounding Input vertical radiative atmospheric profile state.
 * @param lw_heating_rate Output computed longwave cooling rates profile (K/day) [columns, layers].
 * @param lw_flux_down Output computed downward longwave flux profile (W/m2) [columns, layers + 1].
 * @param lw_flux_up Output computed upward longwave flux profile (W/m2) [columns, layers + 1].
 */
void rrtmg_lw_radiation_run(
    size_t columns, size_t layers, size_t bands,
    RadiativeSounding sounding,
    View2D lw_heating_rate,
    View2D lw_flux_down,
    View2D lw_flux_up
);

} // namespace lw
} // namespace rrtmg

#endif // RRTMG_LW_RADIATION_HPP
