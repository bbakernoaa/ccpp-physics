#include "satmedmf_vdifq.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Starting SATMEDMF PBL Lifecycle Test..." << std::endl;

    // Grid Dimensions
    const size_t layers = 10;
    const size_t columns = 5;
    const double dt = 300.0;

    // Allocate flat memories
    std::vector<double> t_lay_data(columns * layers, 280.0);
    std::vector<double> u_wind_data(columns * layers, 5.0);
    std::vector<double> v_wind_data(columns * layers, 2.0);
    std::vector<double> q_vap_data(columns * layers, 0.001);
    std::vector<double> te_data(columns * layers, 0.1);

    std::vector<double> p_lay_data(columns * layers, 90000.0);
    std::vector<double> rho_data(columns * layers, 1.1);

    std::vector<double> heat_data(columns, 0.05);
    std::vector<double> evap_data(columns, 0.0001);
    std::vector<double> stress_data(columns, 0.02);
    std::vector<double> cfch_data(columns, 15.0);

    std::vector<double> dt_temp_data(columns * layers, 0.0);
    std::vector<double> dt_u_data(columns * layers, 0.0);
    std::vector<double> dt_v_data(columns * layers, 0.0);
    std::vector<double> dt_te_data(columns * layers, 0.0);

    // Build mdspan Views
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

    // Invoke solver
    satmedmf::satmedmf_run_core(
        layers,
        columns,
        dt,
        true, // tte_edmf
        state,
        p_lay_view,
        rho_view,
        heat_data.data(),
        evap_data.data(),
        stress_data.data(),
        true, // do_canopy
        cfch_data.data(),
        tendencies
    );

    // Verify output initialized to 0
    for (size_t col = 0; col < columns; ++col) {
        for (size_t lay = 0; lay < layers; ++lay) {
            if (tendencies.dt_temp[col, lay] != 0.0) {
                std::cerr << "Verification failed: dt_temp not initialized to 0 at ["
                          << col << ", " << lay << "]" << std::endl;
                return 1;
            }
        }
    }

    std::cout << "SATMEDMF PBL Lifecycle Test PASSED." << std::endl;
    return 0;
}
