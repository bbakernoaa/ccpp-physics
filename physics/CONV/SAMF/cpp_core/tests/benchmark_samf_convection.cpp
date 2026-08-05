#include "samf_deep_convection.hpp"
#include "samf_shallow_convection.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_samf_convection..." << std::endl;

    const size_t columns = 50000;
    const size_t layers = 127;
    const double dt = 300.0;

    std::vector<double> t_lay_data(columns * layers, 280.0);
    std::vector<double> q_vap_data(columns * layers, 0.005);
    std::vector<double> u_wind_data(columns * layers, 10.0);
    std::vector<double> v_wind_data(columns * layers, 5.0);
    std::vector<double> p_lay_data(columns * layers, 90000.0);
    std::vector<double> p_int_data(columns * (layers + 1), 91000.0);
    std::vector<double> z_lay_data(columns * layers, 500.0);
    std::vector<double> z_int_data(columns * (layers + 1), 100.0);

    std::vector<double> dt_t_data(columns * layers, 0.0);
    std::vector<double> dt_q_data(columns * layers, 0.0);
    std::vector<double> dt_u_data(columns * layers, 0.0);
    std::vector<double> dt_v_data(columns * layers, 0.0);

    // Deep diagnostics allocations
    std::vector<double> ud_mf_data(columns * layers, 0.0);
    std::vector<double> dd_mf_data(columns * layers, 0.0);
    std::vector<double> dt_mf_data(columns * layers, 0.0);
    std::vector<double> cnvw_data(columns * layers, 0.0);
    std::vector<double> cnvc_data(columns * layers, 0.0);
    std::vector<int> kbot_data(columns, 0);
    std::vector<int> ktop_data(columns, 0);
    std::vector<int> kcnv_data(columns, 0);
    std::vector<double> rain_data(columns, 0.0);

    // Shallow diagnostics allocations
    std::vector<double> ud_mf_sh_data(columns * layers, 0.0);
    std::vector<double> dt_mf_sh_data(columns * layers, 0.0);
    std::vector<int> kbot_sh_data(columns, 0);
    std::vector<int> ktop_sh_data(columns, 0);

    samf::View2D t_lay(t_lay_data.data(), columns, layers);
    samf::View2D q_vap(q_vap_data.data(), columns, layers);
    samf::View2D u_wind(u_wind_data.data(), columns, layers);
    samf::View2D v_wind(v_wind_data.data(), columns, layers);

    samf::ConstView2D p_lay(p_lay_data.data(), columns, layers);
    samf::ConstView2D p_int(p_int_data.data(), columns, layers + 1);
    samf::ConstView2D z_lay(z_lay_data.data(), columns, layers);
    samf::ConstView2D z_int(z_int_data.data(), columns, layers + 1);

    samf::View2D dt_t(dt_t_data.data(), columns, layers);
    samf::View2D dt_q(dt_q_data.data(), columns, layers);
    samf::View2D dt_u(dt_u_data.data(), columns, layers);
    samf::View2D dt_v(dt_v_data.data(), columns, layers);

    samf::View2D ud_mf(ud_mf_data.data(), columns, layers);
    samf::View2D dd_mf(dd_mf_data.data(), columns, layers);
    samf::View2D dt_mf(dt_mf_data.data(), columns, layers);
    samf::View2D cnvw(cnvw_data.data(), columns, layers);
    samf::View2D cnvc(cnvc_data.data(), columns, layers);

    samf::View2D ud_mf_sh(ud_mf_sh_data.data(), columns, layers);
    samf::View2D dt_mf_sh(dt_mf_sh_data.data(), columns, layers);

    samf::ConvectiveSounding sounding{
        t_lay,
        q_vap,
        u_wind,
        v_wind,
        p_lay,
        p_int,
        z_lay,
        z_int
    };

    // Warm-up sweep
    samf::deep::samf_deep_convection_run(
        columns, layers, dt, sounding, dt_t, dt_q, dt_u, dt_v,
        ud_mf, dd_mf, dt_mf, cnvw, cnvc,
        kbot_data.data(), ktop_data.data(), kcnv_data.data(),
        rain_data.data()
    );
    samf::shallow::samf_shallow_convection_run(
        columns, layers, dt, sounding, dt_t, dt_q, dt_u, dt_v,
        ud_mf_sh, dt_mf_sh,
        kbot_sh_data.data(), ktop_sh_data.data()
    );

    std::cout << "  - Grid dimension: " << columns << " columns, " << layers << " layers (" << columns * layers << " total cells)" << std::endl;
    std::cout << "  - Running 10 execution sweeps..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 10; ++step) {
        samf::deep::samf_deep_convection_run(
            columns, layers, dt, sounding, dt_t, dt_q, dt_u, dt_v,
            ud_mf, dd_mf, dt_mf, cnvw, cnvc,
            kbot_data.data(), ktop_data.data(), kcnv_data.data(),
            rain_data.data()
        );
        samf::shallow::samf_shallow_convection_run(
            columns, layers, dt, sounding, dt_t, dt_q, dt_u, dt_v,
            ud_mf_sh, dt_mf_sh,
            kbot_sh_data.data(), ktop_sh_data.data()
        );
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count() / 10.0;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << elapsed << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << (columns * layers) / elapsed << " grid cells / sec" << std::endl;
    std::cout << "benchmark_samf_convection PASS!" << std::endl;
    return 0;
}
