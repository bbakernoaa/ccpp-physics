#include "rrtmg_sw_radiation.hpp"
#include "rrtmg_lw_radiation.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_rrtmg_radiation..." << std::endl;

    const size_t columns = 50000;
    const size_t layers = 127;
    const size_t bands_sw = 14;
    const size_t bands_lw = 16;

    std::vector<double> t_lay_data(columns * layers, 250.0);
    std::vector<double> q_vap_data(columns * layers, 1.0e-4);
    std::vector<double> o3_vap_data(columns * layers, 1.0e-6);
    std::vector<double> cld_frac_data(columns * layers, 0.0);
    std::vector<double> p_lay_data(columns * layers, 50000.0);
    std::vector<double> p_int_data(columns * (layers + 1), 51000.0);
    std::vector<double> albedo_data(columns * bands_sw, 0.15);
    std::vector<double> emissivity_data(columns * bands_lw, 0.98);
    std::vector<double> cos_solar_zenith_data(columns, 0.8);

    std::vector<double> sw_heating_rate_data(columns * layers, 0.0);
    std::vector<double> sw_flux_down_data(columns * (layers + 1), 0.0);
    std::vector<double> sw_flux_up_data(columns * (layers + 1), 0.0);

    std::vector<double> lw_heating_rate_data(columns * layers, 0.0);
    std::vector<double> lw_flux_down_data(columns * (layers + 1), 0.0);
    std::vector<double> lw_flux_up_data(columns * (layers + 1), 0.0);

    rrtmg::ConstView2D t_lay(t_lay_data.data(), columns, layers);
    rrtmg::ConstView2D q_vap(q_vap_data.data(), columns, layers);
    rrtmg::ConstView2D o3_vap(o3_vap_data.data(), columns, layers);
    rrtmg::ConstView2D cld_frac(cld_frac_data.data(), columns, layers);
    rrtmg::ConstView2D p_lay(p_lay_data.data(), columns, layers);
    rrtmg::ConstView2D p_int(p_int_data.data(), columns, layers + 1);
    
    rrtmg::ConstView2D albedo(albedo_data.data(), columns, bands_sw);
    rrtmg::ConstView2D emissivity(emissivity_data.data(), columns, bands_lw);

    rrtmg::View2D sw_heating_rate(sw_heating_rate_data.data(), columns, layers);
    rrtmg::View2D sw_flux_down(sw_flux_down_data.data(), columns, layers + 1);
    rrtmg::View2D sw_flux_up(sw_flux_up_data.data(), columns, layers + 1);

    rrtmg::View2D lw_heating_rate(lw_heating_rate_data.data(), columns, layers);
    rrtmg::View2D lw_flux_down(lw_flux_down_data.data(), columns, layers + 1);
    rrtmg::View2D lw_flux_up(lw_flux_up_data.data(), columns, layers + 1);

    rrtmg::RadiativeSounding sounding_sw{
        t_lay, q_vap, o3_vap, cld_frac, p_lay, p_int, albedo
    };

    rrtmg::RadiativeSounding sounding_lw{
        t_lay, q_vap, o3_vap, cld_frac, p_lay, p_int, emissivity
    };

    // Warm-up sweep
    rrtmg::sw::rrtmg_sw_radiation_run(columns, layers, bands_sw, sounding_sw, cos_solar_zenith_data.data(), sw_heating_rate, sw_flux_down, sw_flux_up);
    rrtmg::lw::rrtmg_lw_radiation_run(columns, layers, bands_lw, sounding_lw, lw_heating_rate, lw_flux_down, lw_flux_up);

    std::cout << "  - Grid dimension: " << columns << " columns, " << layers << " layers (" << columns * layers << " total cells)" << std::endl;
    std::cout << "  - Running 10 execution sweeps..." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 10; ++step) {
        rrtmg::sw::rrtmg_sw_radiation_run(columns, layers, bands_sw, sounding_sw, cos_solar_zenith_data.data(), sw_heating_rate, sw_flux_down, sw_flux_up);
        rrtmg::lw::rrtmg_lw_radiation_run(columns, layers, bands_lw, sounding_lw, lw_heating_rate, lw_flux_down, lw_flux_up);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count() / 10.0;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << elapsed << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << (columns * layers) / elapsed << " grid cells / sec" << std::endl;
    std::cout << "benchmark_rrtmg_radiation PASS!" << std::endl;
    return 0;
}
