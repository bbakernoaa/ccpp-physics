#include "satmedmf_vdifq.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Starting SATMEDMF PBL Stride Verification Test..." << std::endl;

    const size_t layers = 4;
    const size_t columns = 3;

    // Allocate flat data array where each element contains its flat-memory index
    std::vector<double> flat_data(columns * layers);
    for (size_t i = 0; i < flat_data.size(); ++i) {
        flat_data[i] = static_cast<double>(i);
    }

    // Wrap in standard LayoutLeft mdspan view
    satmedmf::ConstView2D view(flat_data.data(), columns, layers);

    // Verify layout properties
    if (view.stride(0) != 1) {
        std::cerr << "Stride error: horizontal (columns) dimension stride is "
                  << view.stride(0) << " but must be 1 (contiguous Column-Major)." << std::endl;
        return 1;
    }

    if (view.stride(1) != columns) {
        std::cerr << "Stride error: vertical (layers) dimension stride is "
                  << view.stride(1) << " but must be " << columns << " (the number of columns)." << std::endl;
        return 1;
    }

    // Verify 1-to-1 index matching
    for (size_t lay = 0; lay < layers; ++lay) {
        for (size_t col = 0; col < columns; ++col) {
            const size_t expected_flat_index = col + columns * lay;
            const double expected_val = static_cast<double>(expected_flat_index);
            const double actual_val = view[col, lay];

            if (actual_val != expected_val) {
                std::cerr << "Verification failed: view[" << col << ", " << lay << "] returned "
                          << actual_val << " but expected " << expected_val << std::endl;
                return 1;
            }
        }
    }

    std::cout << "SATMEDMF PBL Stride Verification Test PASSED." << std::endl;
    return 0;
}
