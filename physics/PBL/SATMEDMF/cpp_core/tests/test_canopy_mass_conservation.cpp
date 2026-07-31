#include "satmedmf_canopy.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>

int main() {
    std::cout << "Starting test_canopy_mass_conservation..." << std::endl;

    const size_t layers = 5;
    const size_t columns = 1;
    const size_t nkc = 3;
    const size_t nkt = layers + nkc;
    const size_t ntrac = 1;

    std::vector<double> zi_data(columns * (layers + 1), 0.0);
    std::vector<double> zl_data(columns * layers, 0.0);
    std::vector<double> zm_data(columns * layers, 0.0);
    std::vector<double> dens_data(columns * layers, 1.2);
    
    // Grid coordinate assignments (reversing layers so index 0 is top of the atmosphere)
    for (size_t lay = 0; lay < layers; ++lay) {
        zl_data[lay] = (layers - lay) * 100.0 - 50.0;
        zm_data[lay] = (layers - lay) * 100.0;
        zi_data[lay] = (layers - lay + 1) * 100.0;
    }
    zi_data[layers] = 0.0; // surface interface

    std::vector<double> q1_mod_data(columns * layers, 0.005); // Initial resolved tracer concentrations (5 g/kg)
    std::vector<double> q1_can_data(columns * nkt, 0.0);       // Canopy concentration arrays

    std::vector<int> kmod_data(columns * layers, 0);
    std::vector<int> kcan3_data(columns * nkc, 0);
    std::vector<double> zmid_can_data(columns * nkt, 0.0);
    std::vector<double> zmom_can_data(columns * (nkt + 1), 0.0);
    std::vector<double> prsl_can_data(columns * nkt, 100000.0);
    std::vector<double> dens_can_data(columns * nkt, 1.2);

    // Initialize canopy momentum heights with monotonically decreasing values (top-to-bottom)
    for (size_t k = 0; k <= nkt; ++k) {
        zmom_can_data[k] = (nkt - k) * 100.0;
    }

    satmedmf::ConstView2D zi(zi_data.data(), columns, layers + 1);
    satmedmf::ConstView2D zl(zl_data.data(), columns, layers);
    satmedmf::ConstView2D zm(zm_data.data(), columns, layers);
    satmedmf::View2D q1_mod(q1_mod_data.data(), columns, layers);
    satmedmf::View2D q1_can(q1_can_data.data(), columns, nkt);

    // Initial sum of resolved masses
    double initial_sum = std::accumulate(q1_mod_data.begin(), q1_mod_data.end(), 0.0);

    std::cout << "Initial values in q1_mod:" << std::endl;
    for (size_t k = 0; k < layers; ++k) {
        std::cout << "  q1_mod[0, " << k << "] = " << q1_mod[0, k] << std::endl;
    }

    // 1. Map Resolved -> Canopy (flag = 0)
    satmedmf::canopy::canopy_transfer_run(
        columns, layers, nkc, nkt,
        ntrac, 0, // resolved_to_canopy
        zi, zl, zm,
        q1_mod, q1_can,
        kmod_data, kcan3_data,
        zmid_can_data, zmom_can_data,
        prsl_can_data, dens_can_data
    );

    std::cout << "After resolved_to_canopy, q1_can values:" << std::endl;
    for (size_t k = 0; k < nkt; ++k) {
        std::cout << "  q1_can[0, " << k << "] = " << q1_can[0, k] << std::endl;
    }

    // Zero out resolved matrix
    std::fill(q1_mod_data.begin(), q1_mod_data.end(), 0.0);

    // 2. Map Canopy -> Resolved (flag = 1)
    satmedmf::canopy::canopy_transfer_run(
        columns, layers, nkc, nkt,
        ntrac, 1, // canopy_to_resolved
        zi, zl, zm,
        q1_mod, q1_can,
        kmod_data, kcan3_data,
        zmid_can_data, zmom_can_data,
        prsl_can_data, dens_can_data
    );

    std::cout << "After canopy_to_resolved, q1_mod values:" << std::endl;
    for (size_t k = 0; k < layers; ++k) {
        std::cout << "  q1_mod[0, " << k << "] = " << q1_mod[0, k] << std::endl;
    }

    // Final sum of resolved masses
    double final_sum = std::accumulate(q1_mod_data.begin(), q1_mod_data.end(), 0.0);

    std::cout << "Initial Sum: " << initial_sum << " | Final Sum: " << final_sum << std::endl;

    // Verify perfect mass conservation (0.0% loss)
    double tolerance = 1.0e-12;
    assert((std::abs(final_sum - initial_sum) < tolerance));

    std::cout << "test_canopy_mass_conservation PASSED." << std::endl;
    return 0;
}
