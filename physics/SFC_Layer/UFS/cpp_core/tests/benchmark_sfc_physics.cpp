#include "sfc_diff.hpp"
#include "sfc_nst.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_sfc_physics..." << std::endl;

    const size_t columns = 50000;

    std::vector<double> u1_data(columns, 5.0);
    std::vector<double> v1_data(columns, 3.0);
    std::vector<double> t1_data(columns, 285.15);
    std::vector<double> q1_data(columns, 1.0e-3);
    std::vector<double> z1_data(columns, 10.0);
    std::vector<double> ps_data(columns, 101325.0);
    std::vector<double> tskin_data(columns, 288.15);
    std::vector<double> z0_data(columns, 0.05);

    std::vector<double> cm_data(columns, 0.0);
    std::vector<double> ch_data(columns, 0.0);
    std::vector<double> ustar_data(columns, 0.0);
    std::vector<double> stress_data(columns, 0.0);

    std::vector<double> sol_flux_data(columns, 500.0);
    std::vector<double> tskin_wat_data(columns, 295.15);
    std::vector<double> cool_skin_data(columns, 0.0);
    std::vector<double> warm_layer_data(columns, 0.0);

    sfc::ConstView1D u1(u1_data.data(), columns);
    sfc::ConstView1D v1(v1_data.data(), columns);
    sfc::ConstView1D t1(t1_data.data(), columns);
    sfc::ConstView1D q1(q1_data.data(), columns);
    sfc::ConstView1D z1(z1_data.data(), columns);
    sfc::ConstView1D ps(ps_data.data(), columns);
    sfc::ConstView1D tskin(tskin_data.data(), columns);

    sfc::View1D cm(cm_data.data(), columns);
    sfc::View1D ch(ch_data.data(), columns);
    sfc::View1D ustar(ustar_data.data(), columns);
    sfc::View1D stress(stress_data.data(), columns);

    sfc::View1D tskin_wat(tskin_wat_data.data(), columns);
    sfc::View1D cool_skin(cool_skin_data.data(), columns);
    sfc::View1D warm_layer(warm_layer_data.data(), columns);

    sfc::SurfaceSounding sounding{
        u1, v1, t1, q1, z1, ps, tskin
    };

    // Warm-up sweep
    sfc::diff::sfc_diff_run(columns, sounding, z0_data.data(), 6, cm, ch, ustar, stress);
    sfc::nst::sfc_nst_run(columns, sol_flux_data.data(), stress_data.data(), tskin_wat, cool_skin, warm_layer);

    std::cout << "  - Grid dimension: " << columns << " columns" << std::endl;
    std::cout << "  - Running 50 execution sweeps..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 50; ++step) {
        sfc::diff::sfc_diff_run(columns, sounding, z0_data.data(), 6, cm, ch, ustar, stress);
        sfc::nst::sfc_nst_run(columns, sol_flux_data.data(), stress_data.data(), tskin_wat, cool_skin, warm_layer);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count() / 50.0;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << elapsed << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << columns / elapsed << " columns / sec" << std::endl;
    std::cout << "benchmark_sfc_physics PASS!" << std::endl;
    return 0;
}
