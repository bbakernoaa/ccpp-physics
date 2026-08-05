#include "h2ophys.hpp"
#include <cmath>
#include <algorithm>

/**
 * @file h2ophys.cpp
 * @brief High-performance C++23 std::mdspan implementation of GFS Stratospheric Water Vapor photochemistry (module_h2ophys.F90).
 */

namespace photochem {
namespace h2o {

/**
 * @brief Standalone physical C++23 solver for GFS Stratospheric Water Vapor photochemistry.
 */
void run_h2ophys(
    size_t columns, size_t layers,
    double dt,
    PhotochemSounding sounding,
    ConstView2D h2opltc,
    View2D dqv_dt_prd,
    View2D dqv_dt_qv
) {
    // Thread-safe parallel execution across horizontal columns
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        // Initialize outputs
        for (size_t k = 0; k < layers; ++k) {
            dqv_dt_prd[i, k] = 0.0;
            dqv_dt_qv[i, k] = 0.0;
        }

        for (size_t k = 0; k < layers; ++k) {
            double temp = sounding.t_lay[i, k];
            double press = sounding.p_lay[i, k];
            double h2o_val = sounding.species[i, k];

            // Species clamping: protect local species to remain strictly non-negative
            double h2o_clamped = std::max(0.0, h2o_val);

            // 1. Compute stratospheric methane oxidation water vapor source rate dqv_dt_prd (kg/kg/s)
            // Parameterized source based on methane consumption deficit from upper stratosphere ceiling
            double h2o_deficit = std::max(0.0, 3e-6 - h2o_clamped);
            dqv_dt_prd[i, k] = 1.5e-11 * press * h2o_deficit;

            // 2. Compute stratospheric water vapor photolysis loss sink tendency dqv_dt_qv (kg/kg/s)
            // Photolytic destruction under Lyman-alpha solar ultraviolet bands in mesosphere/upper-stratosphere
            dqv_dt_qv[i, k] = -h2opltc[i, k] * h2o_clamped;
        }
    }
}

} // namespace h2o
} // namespace photochem

extern "C" {

/**
 * @brief Flat C ABI entry point for the C++23 GFS Stratospheric Water Vapor photochemistry solver.
 */
void c_run_h2ophys(
    size_t columns, size_t layers,
    double dt,
    const double* t_lay, const double* p_lay, const double* dp, const double* h2o,
    const double* h2opltc,
    double* dqv_dt_prd,
    double* dqv_dt_qv
) {
    photochem::ConstView2D t_lay_view(t_lay, columns, layers);
    photochem::ConstView2D p_lay_view(p_lay, columns, layers);
    photochem::ConstView2D dp_view(dp, columns, layers);
    photochem::ConstView2D h2o_view(h2o, columns, layers);
    photochem::ConstView2D h2opltc_view(h2opltc, columns, layers);

    photochem::View2D dqv_dt_prd_view(dqv_dt_prd, columns, layers);
    photochem::View2D dqv_dt_qv_view(dqv_dt_qv, columns, layers);

    photochem::PhotochemSounding sounding{
        t_lay_view,
        p_lay_view,
        dp_view,
        h2o_view
    };

    photochem::h2o::run_h2ophys(
        columns, layers,
        dt,
        sounding,
        h2opltc_view,
        dqv_dt_prd_view,
        dqv_dt_qv_view
    );
}

}
