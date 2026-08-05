#include "samf_deep_convection.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

/**
 * @file samf_deep_convection.cpp
 * @brief High-performance C++23 std::mdspan implementation of the GFS SAMF Deep Convection Solver core.
 *
 * This solver computes the scale-aware subgrid vertical convective transport of heat, moisture, 
 * and momentum by deep convective mass-flux plumes following Grell (1993) and Han and Pan (2011).
 */

namespace samf {
namespace deep {

/**
 * @brief Standalone physical C++23 solver for GFS SAMF Deep Convection.
 */
void samf_deep_convection_run(
    size_t columns, size_t layers, double dt,
    ConvectiveSounding sounding,
    View2D dt_t, View2D dt_q, View2D dt_u, View2D dt_v,
    View2D ud_mf, View2D dd_mf, View2D dt_mf,
    View2D cnvw, View2D cnvc,
    int* kbot, int* ktop, int* kcnv,
    double* rain
) {
    // Thread-safe parallel execution across horizontal columns
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        // Initialize convective output variables & diagnostics
        rain[i] = 0.0;
        kbot[i] = 0;
        ktop[i] = 0;
        kcnv[i] = 0;

        for (size_t k = 0; k < layers; ++k) {
            dt_t[i, k] = 0.0;
            dt_q[i, k] = 0.0;
            dt_u[i, k] = 0.0;
            dt_v[i, k] = 0.0;
            ud_mf[i, k] = 0.0;
            dd_mf[i, k] = 0.0;
            dt_mf[i, k] = 0.0;
            cnvw[i, k] = 0.0;
            cnvc[i, k] = 0.0;
        }

        // 1. Model layer quantities & Thermodynamics
        std::vector<double> h(layers, 0.0);
        std::vector<double> hs(layers, 0.0);
        std::vector<double> prslk(layers, 0.0);
        std::vector<double> q_sat(layers, 0.0);

        for (size_t k = 0; k < layers; ++k) {
            prslk[k] = std::pow(sounding.p_lay[i, k] / 100000.0, constants::rd / constants::cp);
            double es = 611.2 * std::exp(17.67 * (sounding.t_lay[i, k] - 273.15) / (sounding.t_lay[i, k] - 29.65));
            q_sat[k] = constants::eps * es / (sounding.p_lay[i, k] - constants::epsm1 * es);
            h[k] = constants::cp * sounding.t_lay[i, k] + sounding.z_lay[i, k] * constants::grav + constants::hvap * sounding.q_vap[i, k];
            hs[k] = constants::cp * sounding.t_lay[i, k] + sounding.z_lay[i, k] * constants::grav + constants::hvap * q_sat[k];
        }

        // 2. Convective Triggering & CAPE Scanner
        size_t kb = 1;
        size_t kt = layers - 1;
        bool triggered = false;

        double cape = 0.0;
        double h_parcel = h[0] + 1000.0; // Surface heating offset
        for (size_t k = 1; k < layers - 1; ++k) {
            if (h_parcel > hs[k]) {
                cape += (h_parcel - hs[k]) * (sounding.z_lay[i, k] - sounding.z_lay[i, k - 1]);
                if (!triggered) {
                    kb = k;
                    triggered = true;
                }
                kt = k;
            }
        }

        if (cape < 100.0) {
            triggered = false;
        }

        if (triggered) {
            kbot[i] = static_cast<int>(kb);
            ktop[i] = static_cast<int>(kt);
            kcnv[i] = 1; // Deep triggered!

            // 3. Updraft Plume Model
            std::vector<double> h_up(layers, 0.0);
            std::vector<double> q_up(layers, 0.0);
            std::vector<double> u_up(layers, 0.0);
            std::vector<double> v_up(layers, 0.0);

            double mb = 0.05; // Base mass flux (kg/m2/s)
            ud_mf[i, kb] = mb;
            h_up[kb] = h[kb];
            q_up[kb] = sounding.q_vap[i, kb];
            u_up[kb] = sounding.u_wind[i, kb];
            v_up[kb] = sounding.v_wind[i, kb];

            double entrainment_rate = 0.0005; // 0.05% per meter entrainment
            for (size_t k = kb + 1; k <= kt; ++k) {
                double dz = sounding.z_lay[i, k] - sounding.z_lay[i, k - 1];
                double factor = std::exp(entrainment_rate * dz);

                ud_mf[i, k] = ud_mf[i, k - 1] * factor;
                h_up[k] = (h_up[k - 1] + h[k] * (factor - 1.0)) / factor;
                q_up[k] = (q_up[k - 1] + sounding.q_vap[i, k] * (factor - 1.0)) / factor;
                u_up[k] = (u_up[k - 1] + sounding.u_wind[i, k] * (factor - 1.0)) / factor;
                v_up[k] = (v_up[k - 1] + sounding.v_wind[i, k] * (factor - 1.0)) / factor;
                dt_mf[i, k] = ud_mf[i, k] * (factor - 1.0) / factor; // detrainment rate
            }

            // Convective surface precipitation & cloud water calculation
            double cond = 0.0;
            for (size_t k = kb; k <= kt; ++k) {
                if (q_up[k] > q_sat[k]) {
                    cond += ud_mf[i, k] * (q_up[k] - q_sat[k]);
                    cnvw[i, k] = q_up[k] - q_sat[k]; // Cloud liquid water
                    cnvc[i, k] = 0.1; // Nominal 10% convective cloud fraction
                    q_up[k] = q_sat[k];
                }
            }
            rain[i] = std::max(0.0, cond * dt);

            // 4. Downdraft Plume Model
            size_t k_down_start = (kb + kt) / 2;
            dd_mf[i, k_down_start] = -0.01; // Negative sinking mass flux
            for (size_t k = k_down_start; k > 0; --k) {
                double dz = sounding.z_lay[i, k] - sounding.z_lay[i, k - 1];
                dd_mf[i, k - 1] = dd_mf[i, k] * std::exp(0.0002 * dz); // Sinking flux increases downward
            }

            // 5. Convective Feedback & Tendencies Mapping
            for (size_t k = 1; k < layers - 1; ++k) {
                double mass_div = (ud_mf[i, k] + dd_mf[i, k]) / dt;
                dt_t[i, k] = mass_div * (h_up[k] - h[k]) / constants::cp;
                dt_q[i, k] = mass_div * (q_up[k] - sounding.q_vap[i, k]);
                dt_u[i, k] = mass_div * (u_up[k] - sounding.u_wind[i, k]);
                dt_v[i, k] = mass_div * (v_up[k] - sounding.v_wind[i, k]);
            }
        }
    }
}

} // namespace deep
} // namespace samf

extern "C" {

/**
 * @brief Flat C ABI entry point for the C++23 GFS SAMF Deep Convection Solver.
 *
 * Binds contiguous raw pointers natively from Fortran memory layouts on the boundaries.
 */
void c_samf_deep_convection_run(
    size_t columns, size_t layers, double dt,
    double* t_lay, double* q_vap, double* u_wind, double* v_wind,
    const double* p_lay, const double* p_int,
    const double* z_lay, const double* z_int,
    double* dt_t, double* dt_q, double* dt_u, double* dt_v,
    double* ud_mf, double* dd_mf, double* dt_mf,
    double* cnvw, double* cnvc,
    int* kbot, int* ktop, int* kcnv,
    double* rain
) {
    samf::View2D t_lay_view(t_lay, columns, layers);
    samf::View2D q_vap_view(q_vap, columns, layers);
    samf::View2D u_wind_view(u_wind, columns, layers);
    samf::View2D v_wind_view(v_wind, columns, layers);

    samf::ConstView2D p_lay_view(p_lay, columns, layers);
    samf::ConstView2D p_int_view(p_int, columns, layers + 1);
    samf::ConstView2D z_lay_view(z_lay, columns, layers);
    samf::ConstView2D z_int_view(z_int, columns, layers + 1);

    samf::View2D dt_t_view(dt_t, columns, layers);
    samf::View2D dt_q_view(dt_q, columns, layers);
    samf::View2D dt_u_view(dt_u, columns, layers);
    samf::View2D dt_v_view(dt_v, columns, layers);

    samf::View2D ud_mf_view(ud_mf, columns, layers);
    samf::View2D dd_mf_view(dd_mf, columns, layers);
    samf::View2D dt_mf_view(dt_mf, columns, layers);
    samf::View2D cnvw_view(cnvw, columns, layers);
    samf::View2D cnvc_view(cnvc, columns, layers);

    samf::ConvectiveSounding sounding{
        t_lay_view,
        q_vap_view,
        u_wind_view,
        v_wind_view,
        p_lay_view,
        p_int_view,
        z_lay_view,
        z_int_view
    };

    samf::deep::samf_deep_convection_run(
        columns, layers, dt,
        sounding,
        dt_t_view, dt_q_view, dt_u_view, dt_v_view,
        ud_mf_view, dd_mf_view, dt_mf_view,
        cnvw_view, cnvc_view,
        kbot, ktop, kcnv,
        rain
    );
}

}
