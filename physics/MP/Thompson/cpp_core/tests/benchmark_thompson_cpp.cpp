#include "thompson_microphysics_interface.hpp"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    std::cout << "Starting benchmark_thompson_cpp..." << std::endl;

    size_t layers = 128;
    size_t columns = 50000; // Huge grid size to profile exascale performance
    size_t size = layers * columns;

    std::vector<double> t_lay(size, 260.0);
    std::vector<double> p_lay(size, 80000.0);
    std::vector<double> rho(size, 1.0);

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

    std::cout << "  - Grid dimension: " << columns << " columns, " << layers << " layers (" << size << " total cells)" << std::endl;
    std::cout << "  - Running 10 execution sweeps..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    for (int run = 0; run < 10; ++run) {
        c_thompson_microphysics_run(
            layers, columns, 1.0,
            t_lay.data(), p_lay.data(), rho.data(),
            qv.data(), qc.data(), qr.data(), qi.data(), qs.data(), qg.data(),
            ni.data(), nr.data(), ns.data(), ng.data(),
            precip.data()
        );
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    double avg_time = diff.count() / 10.0;
    double cells_per_sec = static_cast<double>(size) / avg_time;

    std::cout << "Benchmark complete:" << std::endl;
    std::cout << "  - Avg run-time: " << avg_time << " seconds" << std::endl;
    std::cout << "  - Throughput:   " << cells_per_sec << " grid cells / sec" << std::endl;
    std::cout << "benchmark_thompson_cpp PASS!" << std::endl;
    return 0;
}
