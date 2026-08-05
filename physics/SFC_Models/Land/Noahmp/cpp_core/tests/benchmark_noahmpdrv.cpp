#include "noahmpdrv.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_noahmpdrv..." << std::endl;

    const size_t columns = 50000;
    const size_t soil_layers = 4;

    std::vector<double> stc_data(columns * soil_layers, 280.0);
    std::vector<double> smc_data(columns * soil_layers, 0.35);
    std::vector<double> sh2o_data(columns * soil_layers, 0.35);
    std::vector<double> sldpst_data(columns * soil_layers, 0.1);
    std::vector<double> tg_data(columns, 288.15);
    std::vector<double> tv_data(columns, 287.15);

    std::vector<double> sfctmp_data(columns, 285.15);
    std::vector<double> sfcprs_data(columns, 101325.0);
    std::vector<double> q2_data(columns, 4.0e-3);
    std::vector<double> soldn_data(columns, 600.0);
    std::vector<double> lwdn_data(columns, 300.0);
    
    std::vector<double> wind_data(columns, 5.0);
    std::vector<double> ch_data(columns, 0.005);
    std::vector<int> is_glacier_data(columns, 0);

    std::vector<double> sheat_data(columns, 0.0);
    std::vector<double> eta_data(columns, 0.0);
    std::vector<double> gflux_data(columns, 0.0);
    std::vector<double> runoff_data(columns, 0.0);

    noahmp::ConstView2D stc(stc_data.data(), columns, soil_layers);
    noahmp::ConstView2D smc(smc_data.data(), columns, soil_layers);
    noahmp::ConstView2D sh2o(sh2o_data.data(), columns, soil_layers);
    noahmp::ConstView2D sldpst(sldpst_data.data(), columns, soil_layers);
    noahmp::ConstView1D tg(tg_data.data(), columns);
    noahmp::ConstView1D tv(tv_data.data(), columns);

    noahmp::View1D sheat(sheat_data.data(), columns);
    noahmp::View1D eta(eta_data.data(), columns);
    noahmp::View1D gflux(gflux_data.data(), columns);
    noahmp::View1D runoff(runoff_data.data(), columns);

    noahmp::LandSounding sounding{
        stc, smc, sh2o, sldpst, tg, tv
    };

    // Configure all 19 scientific parameter options
    noahmp::NoahMP_Config config{
        1, 1, 1, 1, 1,  // idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc
        1, 1, 1, 1, 1,  // iopt_frz, iopt_inf, iopt_rad, iopt_alb, iopt_snf
        1, 1, 1, 1, 4,  // iopt_tbot, iopt_stc, iopt_trs, iopt_diag, iopt_rsf
        1, 1, 0, 2, 1   // iopt_soil, iopt_pedo, iopt_crop, iopt_gla, iopt_z0m
    };

    // Warm-up sweep
    noahmp::noahmp_sflx_run(columns, soil_layers, 1800.0, sounding,
        sfctmp_data.data(), sfcprs_data.data(), q2_data.data(), soldn_data.data(), lwdn_data.data(),
        wind_data.data(), ch_data.data(), is_glacier_data.data(),
        config,
        sheat, eta, gflux, runoff);

    std::cout << "  - Grid dimension: " << columns << " columns (" << columns * soil_layers << " total soil cells)" << std::endl;
    std::cout << "  - Running 50 execution sweeps..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 50; ++step) {
        noahmp::noahmp_sflx_run(columns, soil_layers, 1800.0, sounding,
            sfctmp_data.data(), sfcprs_data.data(), q2_data.data(), soldn_data.data(), lwdn_data.data(),
            wind_data.data(), ch_data.data(), is_glacier_data.data(),
            config,
            sheat, eta, gflux, runoff);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count() / 50.0;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << elapsed << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << columns / elapsed << " columns / sec" << std::endl;
    std::cout << "benchmark_noahmpdrv PASS!" << std::endl;
    return 0;
}
