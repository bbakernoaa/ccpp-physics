#include "satmedmf_vdifq.hpp"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    std::cout << "Starting test_canopy under sub-grid forest canopy active profile..." << std::endl;

    const size_t layers = 5;
    const size_t columns = 1;
    const double dt = 300.0;

    std::vector<double> t_lay_data = { 290.0, 288.0, 286.0, 284.0, 282.0 };
    std::vector<double> u_wind_data = { 5.0, 5.0, 5.0, 5.0, 5.0 };
    std::vector<double> v_wind_data = { 2.0, 2.0, 2.0, 2.0, 2.0 };
    std::vector<double> q_vap_data = { 0.010, 0.009, 0.008, 0.007, 0.006 };
    std::vector<double> te_data = { 1.0, 0.8, 0.6, 0.4, 0.2 };

    std::vector<double> p_lay_data = { 100000.0, 99000.0, 98000.0, 97000.0, 96000.0 };
    std::vector<double> rho_data = { 1.2, 1.21, 1.22, 1.23, 1.24 };

    std::vector<double> heat_data = { 0.05 };
    std::vector<double> evap_data = { 0.0001 };
    std::vector<double> stress_data = { 0.01 };
    std::vector<double> cfch_data = { 15.0 }; // Canopy Forest Height = 15m (first layer is at ~41m height, let's make cfch 100m to cover first layer centers for test verification)
    cfch_data[0] = 100.0;

    std::vector<double> dt_temp_data(columns * layers, 0.0);
    std::vector<double> dt_u_data(columns * layers, 0.0);
    std::vector<double> dt_v_data(columns * layers, 0.0);
    std::vector<double> dt_te_data(columns * layers, 0.0);

    satmedmf::View2D t_lay_view(t_lay_data.data(), columns, layers);
    satmedmf::View2D u_wind_view(u_wind_data.data(), columns, layers);
    satmedmf::View2D v_wind_view(v_wind_data.data(), columns, layers);
    satmedmf::View2D q_vap_view(q_vap_data.data(), columns, layers);
    satmedmf::View2D te_view(te_data.data(), columns, layers);

    satmedmf::ConstView2D p_lay_view(p_lay_data.data(), columns, layers);
    satmedmf::ConstView2D rho_view(rho_data.data(), columns, layers);

    satmedmf::View2D dt_temp_view(dt_temp_data.data(), columns, layers);
    satmedmf::View2D dt_u_view(dt_u_data.data(), columns, layers);
    satmedmf::View2D dt_v_view(dt_v_data.data(), columns, layers);
    satmedmf::View2D dt_te_view(dt_te_data.data(), columns, layers);

    satmedmf::PhysicalState state{
        t_lay_view,
        u_wind_view,
        v_wind_view,
        q_vap_view,
        te_view
    };

    satmedmf::DiagnosticTendencies tendencies{
        dt_temp_view,
        dt_u_view,
        dt_v_view,
        dt_te_view
    };

    // Run core solver with do_canopy = true
    satmedmf::satmedmf_run_core(
        layers,
        columns,
        dt,
        false, // Standard TKE
        state,
        p_lay_view,
        rho_view,
        heat_data.data(),
        evap_data.data(),
        stress_data.data(),
        true, // do_canopy ACTIVE
        cfch_data.data(),
        tendencies
    );

    // Verify subgrid canopy drag cooling was applied to the in-canopy layers (lay = 1, since zl1 ~41m < 100m)
    std::cout << "In-Canopy U-wind Tendency at Layer 1: " << tendencies.dt_u[0, 1] << " m/s2" << std::endl;
    assert((tendencies.dt_u[0, 1] < 0.0));

    std::cout << "test_canopy PCC-Canopy PASSED." << std::endl;
    return 0;
}
