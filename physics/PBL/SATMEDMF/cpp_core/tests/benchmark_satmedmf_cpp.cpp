#include "satmedmf_vdifq.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>

int main() {
    std::cout << "=============================================================" << std::endl;
    std::cout << "   CCPP SATMEDMF C++ 200,000-CELL HIGH-RESOLUTION BENCHMARK  " << std::endl;
    std::cout << "           Grid Dimensions: 2000 Columns x 100 Layers        " << std::endl;
    std::cout << "=============================================================" << std::endl;

    const size_t columns = 2000;
    const size_t layers = 100;
    const double dt = 300.0;

    // Allocate data
    std::vector<double> t_lay_data(columns * layers);
    std::vector<double> u_wind_data(columns * layers);
    std::vector<double> v_wind_data(columns * layers);
    std::vector<double> q_vap_data(columns * layers);
    std::vector<double> te_data(columns * layers);
    std::vector<double> p_lay_data(columns * layers);
    std::vector<double> rho_data(columns * layers);

    std::vector<double> heat_data(columns, 0.05);
    std::vector<double> evap_data(columns, 0.0001);
    std::vector<double> stress_data(columns, 0.02);
    std::vector<double> cfch_data(columns, 100.0);

    std::vector<double> dt_temp_data(columns * layers, 0.0);
    std::vector<double> dt_u_data(columns * layers, 0.0);
    std::vector<double> dt_v_data(columns * layers, 0.0);
    std::vector<double> dt_te_data(columns * layers, 0.0);

    // Initialize inputs matching Fortran benchmark
    for (size_t col = 0; col < columns; ++col) {
        for (size_t lay = 0; lay < layers; ++lay) {
            size_t idx = col + columns * lay;
            t_lay_data[idx] = 288.0 - 0.1 * (lay + 1);
            te_data[idx] = 1.0 - 0.001 * (lay + 1);
            p_lay_data[idx] = 100000.0 - 100.0 * (lay + 1);
            rho_data[idx] = 1.2 - 0.001 * (lay + 1);
            u_wind_data[idx] = 5.0 + 0.02 * (lay + 1);
            v_wind_data[idx] = 2.0;
            q_vap_data[idx] = 0.005;
        }
    }

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

    // Warm-up to ensure CPU caches and static vector caches are initialized
    {
        omp_set_num_threads(1);
        satmedmf::PhysicalState state{ t_lay_view, u_wind_view, v_wind_view, q_vap_view, te_view };
        satmedmf::DiagnosticTendencies tendencies{ dt_temp_view, dt_u_view, dt_v_view, dt_te_view };
        for (int i = 0; i < 3; ++i) {
            satmedmf::satmedmf_run_core(layers, columns, dt, true, state, p_lay_view, rho_view,
                                       heat_data.data(), evap_data.data(), stress_data.data(),
                                       true, cfch_data.data(), tendencies);
        }
    }

    // Benchmark across various thread counts
    std::vector<int> thread_counts = {1, 2, 4, 8, 12};
    for (int t : thread_counts) {
        omp_set_num_threads(t);

        // Reset outputs
        std::fill(dt_temp_data.begin(), dt_temp_data.end(), 0.0);
        std::fill(dt_u_data.begin(), dt_u_data.end(), 0.0);
        std::fill(dt_v_data.begin(), dt_v_data.end(), 0.0);
        std::fill(dt_te_data.begin(), dt_te_data.end(), 0.0);

        satmedmf::PhysicalState state{ t_lay_view, u_wind_view, v_wind_view, q_vap_view, te_view };
        satmedmf::DiagnosticTendencies tendencies{ dt_temp_view, dt_u_view, dt_v_view, dt_te_view };

        auto start = std::chrono::high_resolution_clock::now();
        const int iterations = 10;
        for (int iter = 0; iter < iterations; ++iter) {
            satmedmf::satmedmf_run_core(layers, columns, dt, true, state, p_lay_view, rho_view,
                                       heat_data.data(), evap_data.data(), stress_data.data(),
                                       true, cfch_data.data(), tendencies);
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        double avg_time = elapsed.count() / iterations;

        std::cout << "  C++ Threads: " << t << " | Average Execution Time: " << avg_time << "s" << std::endl;
    }

    std::cout << "=============================================================" << std::endl;
    return 0;
}
