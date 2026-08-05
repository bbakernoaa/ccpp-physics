#ifndef SFC_DIFF_HPP
#define SFC_DIFF_HPP

#include "sfc_types.hpp"
#include "sfc_constants.hpp"
#include <stddef.h>

namespace sfc {
namespace diff {

/**
 * @brief Public interface for GFS Surface Diffusion Exchange Coefficients Solver (sfc_diff_run).
 *
 * Runs stability profile iterations and computes turbulent exchange coefficients
 * for momentum and heat.
 *
 * @param columns Number of horizontal grid columns.
 * @param sounding Input boundary sounding profiles state.
 * @param z0 Input surface roughness length values [columns].
 * @param cm Output momentum exchange coefficient profile [columns].
 * @param ch Output heat/moisture exchange coefficient profile [columns].
 * @param ustar Output friction velocity profile (m/s) [columns].
 * @param stress Output surface wind stress profile (N/m2) [columns].
 */
void sfc_diff_run(
    size_t columns,
    SurfaceSounding sounding,
    const Real* z0,
    int sfc_z0_type,
    View1D cm,
    View1D ch,
    View1D ustar,
    View1D stress
);

} // namespace diff
} // namespace sfc

#endif // SFC_DIFF_HPP
