#include "thompson_microphysics_interface.hpp"
#include <iostream>
#include <vector>

#if __has_include(<mdspan>)
#include <mdspan>
#define HAS_MDSPAN 1
#else
#define HAS_MDSPAN 0
#endif

int main() {
    std::cout << "Running test_lifecycle..." << std::endl;

    // 1. Verify standard C++23 mdspan header exists
    #if HAS_MDSPAN
        std::cout << "  - Standard C++23 <mdspan> header is AVAILABLE." << std::endl;
    #else
        std::cout << "  - Standard C++23 <mdspan> header is NOT AVAILABLE." << std::endl;
        return 1;
    #endif

    // 2. Dry-run call to verify C++ entry point resolution
    size_t layers = 128;
    size_t columns = 1000;
    std::vector<double> dummy(layers * columns, 290.0);
    std::vector<double> precip(3, 0.0);

    c_thompson_microphysics_run(
        layers, columns, 1.0,
        dummy.data(), dummy.data(), dummy.data(),
        dummy.data(), dummy.data(), dummy.data(),
        dummy.data(), dummy.data(), dummy.data(),
        dummy.data(), dummy.data(), dummy.data(),
        dummy.data(),
        precip.data()
    );

    std::cout << "test_lifecycle PASS!" << std::endl;
    return 0;
}
