#include "ozphys.hpp"
#include "h2ophys.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_photochem..." << std::endl;

    const size_t columns = 50000;
    const size_t layers = 127;

    std::vector<double> t_lay_data(columns * layers, 240.0);
    std::vector<double> p_lay_data(columns * layers, 500.0);
    std::vector<double> dp_data(columns * layers, 100.0);
    std::vector<double> species_oz_data(columns * layers, 4.0e-6);
    std::vector<double> species_h2o_data(columns * layers, 2.0e-6);

    std::vector<double> ozpl_data(columns * layers, 2.0e-11);
    std::vector<double> h2opltc_data(columns * layers, 1.0e-7);

    std::vector<double> do3_dt_prd_data(columns * layers, 0.0);
    std::vector<double> do3_dt_temp_data(columns * layers, 0.0);

    std::vector<double> dqv_dt_prd_data(columns * layers, 0.0);
    std::vector<double> dqv_dt_qv_data(columns * layers, 0.0);

    photochem::ConstView2D t_lay(t_lay_data.data(), columns, layers);
    photochem::ConstView2D p_lay(p_lay_data.data(), columns, layers);
    photochem::ConstView2D dp(dp_data.data(), columns, layers);
    photochem::ConstView2D species_oz(species_oz_data.data(), columns, layers);
    photochem::ConstView2D species_h2o(species_h2o_data.data(), columns, layers);

    photochem::ConstView2D ozpl(ozpl_data.data(), columns, layers);
    photochem::ConstView2D h2opltc(h2opltc_data.data(), columns, layers);

    photochem::View2D do3_dt_prd(do3_dt_prd_data.data(), columns, layers);
    photochem::View2D do3_dt_temp(do3_dt_temp_data.data(), columns, layers);

    photochem::View2D dqv_dt_prd(dqv_dt_prd_data.data(), columns, layers);
    photochem::View2D dqv_dt_qv(dqv_dt_qv_data.data(), columns, layers);

    photochem::PhotochemSounding sounding_oz{
        t_lay, p_lay, dp, species_oz
    };

    photochem::PhotochemSounding sounding_h2o{
        t_lay, p_lay, dp, species_h2o
    };

    // Warm-up sweep
    photochem::ozone::run_o3prog_2015(columns, layers, 1.0/9.80665, 300.0, sounding_oz, ozpl, do3_dt_prd, do3_dt_temp);
    photochem::h2o::run_h2ophys(columns, layers, 300.0, sounding_h2o, h2opltc, dqv_dt_prd, dqv_dt_qv);

    std::cout << "  - Grid dimension: " << columns << " columns, " << layers << " layers (" << columns * layers << " total cells)" << std::endl;
    std::cout << "  - Running 20 execution sweeps..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 20; ++step) {
        photochem::ozone::run_o3prog_2015(columns, layers, 1.0/9.80665, 300.0, sounding_oz, ozpl, do3_dt_prd, do3_dt_temp);
        photochem::h2o::run_h2ophys(columns, layers, 300.0, sounding_h2o, h2opltc, dqv_dt_prd, dqv_dt_qv);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count() / 20.0;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << elapsed << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << (columns * layers) / elapsed << " cells / sec" << std::endl;
    std::cout << "benchmark_photochem PASS!" << std::endl;
    return 0;
}
