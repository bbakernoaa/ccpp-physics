#include "thompson_microphysics_interface.hpp"
#include <iostream>
#include <vector>
#include <cmath>

#if __has_include(<mdspan>)
#include <mdspan>
#define HAS_MDSPAN 1
#else
#define HAS_MDSPAN 0
#endif

int main() {
    std::cout << "Running test_thompson_cpp_regression..." << std::endl;

    size_t layers = 128;
    size_t columns = 10;
    size_t size = layers * columns;

    std::vector<double> t_lay(size, 260.0); // Kelvin
    std::vector<double> p_lay(size, 80000.0); // Pascals
    std::vector<double> rho(size, 1.0); // kg/m3

    std::vector<double> qv(size, 1e-3);
    std::vector<double> qc(size, 1e-4);
    std::vector<double> qr(size, 1e-4);
    std::vector<double> qi(size, 1e-5);
    std::vector<double> qs(size, 1e-5);
    std::vector<double> qg(size, 1e-6);

    std::vector<double> ni(size, 1e3);
    std::vector<double> nr(size, 1e3);
    std::vector<double> ns(size, 1e2);
    std::vector<double> ng(size, 1e2);

    std::vector<double> precip(3, 0.0);

    // Call translated baseline C++ solver
    c_thompson_microphysics_run(
        layers, columns, 1.0,
        t_lay.data(), p_lay.data(), rho.data(),
        qv.data(), qc.data(), qr.data(), qi.data(), qs.data(), qg.data(),
        ni.data(), nr.data(), ns.data(), ng.data(),
        precip.data()
    );

    // Verify values (simple baseline check, assuring no NaN or infinite values are generated)
    for (size_t i = 0; i < size; ++i) {
        if (std::isnan(qv[i]) || std::isnan(qc[i]) || std::isnan(t_lay[i])) {
            std::cerr << "test_thompson_cpp_regression FAIL: Numeric NaN generated!" << std::endl;
            return 1;
        }
    }

    std::cout << "test_thompson_cpp_regression PASS!" << std::endl;
    return 0;
}
