#ifndef SFC_NST_HPP
#define SFC_NST_HPP

#include "sfc_types.hpp"
#include "sfc_constants.hpp"
#include <stddef.h>

namespace sfc {
namespace nst {

/**
 * @brief Public interface for GFS Near-Surface Sea Temperature Ocean Solver (sfc_nst_run).
 *
 * Runs the ocean diurnal warming skin layers and thermal profiles.
 *
 * @param columns Number of horizontal grid columns.
 * @param sol_flux Input shortwave solar warming flux (W/m2) [columns].
 * @param wind_stress Input surface wind stress magnitude (N/m2) [columns].
 * @param tskin_wat Output computed sea skin temperatures (K) [columns].
 * @param cool_skin Output computed cool skin amplitude (K) [columns].
 * @param warm_layer Output computed warm layer amplitude (K) [columns].
 */
void sfc_nst_run(
    size_t columns,
    const Real* sol_flux,
    const Real* wind_stress,
    View1D tskin_wat,
    View1D cool_skin,
    View1D warm_layer
);

} // namespace nst
} // namespace sfc

#endif // SFC_NST_HPP
