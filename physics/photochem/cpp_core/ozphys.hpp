#ifndef OZPHYS_HPP
#define OZPHYS_HPP

#include "photochem_types.hpp"
#include "photochem_constants.hpp"
#include <stddef.h>

namespace photochem {
namespace ozone {

/**
 * @brief Public interface for GFS Cariolle Parameterized Ozone Chemistry Solver (run_o3prog_2015).
 *
 * Runs local Cariolle parameterized ozone production, loss, and temperature feedback equations.
 *
 * @param columns Number of horizontal grid columns.
 * @param layers Number of vertical model layers.
 * @param con_1ovg Reciprocal of gravity scaling constant (s2/m).
 * @param dt Physics timestep (s).
 * @param sounding Input vertical radiative species profile state.
 * @param ozpl Input diagnostic ozone production and loss climatology [columns, layers].
 * @param do3_dt_prd Output computed ozone production rate profile (kg/kg/s) [columns, layers].
 * @param do3_dt_temp Output computed temperature feedback tendency profile (kg/kg/s) [columns, layers].
 */
void run_o3prog_2015(
    size_t columns, size_t layers,
    double con_1ovg, double dt,
    PhotochemSounding sounding,
    ConstView2D ozpl,
    View2D do3_dt_prd,
    View2D do3_dt_temp
);

} // namespace ozone
} // namespace photochem

#endif // OZPHYS_HPP
