#include "satmedmf_vdifq.hpp"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    std::cout << "Starting test_updrafts under Convective Boundary Layer (CBL) profile..." << std::endl;

    const size_t layers = 5;
    const size_t columns = 1;
    const double dt = 300.0;

    // Unstable atmosphere: temperature decreases with height (buoyant plume updrafts)
    std::vector<double> t_lay_data = { 295.0, 292.0, 289.0, 286.0, 283.0 };
    std::vector<double> u_wind_data = { 2.0, 2.0, 2.0, 2.0, 2.0 };
    std::vector<double> v_wind_data = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector<double> q_vap_data = { 0.015, 0.012, 0.009, 0.006, 0.003 };
    std::vector<double> te_data = { 1.5, 1.2, 0.9, 0.6, 0.3 };

    std::vector<double> p_lay_data = { 100000.0, 99000.0, 98000.0, 97000.0, 96000.0 };
    std::vector<double> rho_data = { 1.15, 1.16, 1.17, 1.18, 1.19 };

    std::vector<double> heat_data = { 0.15 }; // Positive heat flux at surface (buoyant heating)
    std::vector<double> evap_data = { 0.0001 };
    std::vector<double> stress_data = { 0.01 };
    std::vector<double> cfch_data = { 0.0 };

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

    // Run core solver
    satmedmf::satmedmf_run_core(
        layers,
        columns,
        dt,
        true, // tte_edmf (TTE)
        state,
        p_lay_view,
        rho_view,
        heat_data.data(),
        evap_data.data(),
        stress_data.data(),
        false, // do_canopy
        cfch_data.data(),
        tendencies
    );

    // Verify unstable surface boundary temperature tendency is positive (heating)
    std::cout << "Surface Temp Tendency: " << tendencies.dt_temp[0, 0] << " K/s" << std::endl;
    assert((tendencies.dt_temp[0, 0] > 0.0));

    std::cout << "test_updrafts CBL PASSED." << std::endl;
    return 0;
}
