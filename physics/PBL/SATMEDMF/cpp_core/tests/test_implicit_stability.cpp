#include "satmedmf_tridiagonal.hpp"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    std::cout << "Starting test_implicit_stability under extreme dt=300s profile..." << std::endl;

    const size_t layers = 5;
    const size_t columns = 1;

    // Allocate data
    std::vector<double> cl_data(columns * layers, -0.5);
    std::vector<double> cm_data(columns * layers, 2.0);
    std::vector<double> cu_data(columns * layers, -0.5);
    std::vector<double> r1_data(columns * layers, 10.0);
    std::vector<double> r2_data(columns * layers, 5.0);
    
    std::vector<double> au_data(columns * layers, 0.0);
    std::vector<double> a1_data(columns * layers, 0.0);
    std::vector<double> a2_data(columns * layers, 0.0);

    // Create spans
    satmedmf::ConstView2D cl(cl_data.data(), columns, layers);
    satmedmf::ConstView2D cm(cm_data.data(), columns, layers);
    satmedmf::ConstView2D cu(cu_data.data(), columns, layers);
    satmedmf::ConstView2D r1(r1_data.data(), columns, layers);
    satmedmf::ConstView2D r2(r2_data.data(), columns, layers);

    satmedmf::View2D au(au_data.data(), columns, layers);
    satmedmf::View2D a1(a1_data.data(), columns, layers);
    satmedmf::View2D a2(a2_data.data(), columns, layers);

    // Run core implicit solver
    satmedmf::tridiagonal::tridi2(
        columns, layers,
        cl, cm, cu, r1, r2, au, a1, a2
    );

    // Verify stability (no NaN or extremely large numbers)
    for (size_t col = 0; col < columns; ++col) {
        for (size_t lay = 0; lay < layers; ++lay) {
            assert((a1[col, lay] > 0.0 && a1[col, lay] < 100.0));
            assert((a2[col, lay] > 0.0 && a2[col, lay] < 100.0));
        }
    }

    std::cout << "test_implicit_stability PASSED." << std::endl;
    return 0;
}
