#ifndef H2OPHYS_HPP
#define H2OPHYS_HPP

#include "photochem_types.hpp"
#include "photochem_constants.hpp"
#include <stddef.h>

namespace photochem {
namespace h2o {

/**
 * @brief Public C++23 interface for GFS Stratospheric Water Vapor Chemistry Solver.
 *
 * Runs stratospheric methane oxidation source reactions and water vapor photolysis sink profiles.
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param dt Physics timestep (s).
 * @param sounding Input vertical species profile state.
 * @param h2opltc Input water vapor photolysis loss coefficients [columns, layers].
 * @param dqv_dt_prd Output computed stratospheric water vapor production rate (kg/kg/s) [columns, layers].
 * @param dqv_dt_qv Output computed stratospheric water vapor photolysis sink tendency (kg/kg/s) [columns, layers].
 */
void run_h2ophys(
    size_t columns, size_t layers,
    double dt,
    PhotochemSounding sounding,
    ConstView2D h2opltc,
    View2D dqv_dt_prd,
    View2D dqv_dt_qv
);

} // namespace h2o
} // namespace photochem

#endif // H2OPHYS_HPP
