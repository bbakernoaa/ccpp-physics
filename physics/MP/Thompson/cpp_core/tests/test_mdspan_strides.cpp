#include "thompson_microphysics.hpp"
#include <iostream>
#include <vector>
#include <cassert>

#if __has_include(<mdspan>)
#include <mdspan>
#define HAS_MDSPAN 1
#else
#define HAS_MDSPAN 0
#endif

int main() {
    std::cout << "Running test_mdspan_strides..." << std::endl;

    #if HAS_MDSPAN
        size_t layers = 128;
        size_t columns = 1000;
        std::vector<double> dummy(layers * columns, 0.0);

        // Populate sequential grid index values to verify layout
        for (size_t col = 0; col < columns; ++col) {
            for (size_t lay = 0; lay < layers; ++lay) {
                dummy[col + lay * columns] = static_cast<double>(col + lay * columns);
            }
        }

        // Wrap contiguous flat pointer in standard C++23 LayoutLeft (Column-Major) view (columns, layers)
        thompson::View2D view(dummy.data(), columns, layers);

        // Verify indexing maps exactly to Fortran column-major contiguous struts
        for (size_t col = 0; col < columns; ++col) {
            for (size_t lay = 0; lay < layers; ++lay) {
                double expected = static_cast<double>(col + lay * columns);
                double actual = view[col, lay];
                if (actual != expected) {
                    std::cerr << "test_mdspan_strides FAIL at indices (" << col << ", " << lay << ")" << std::endl;
                    std::cerr << "  Expected: " << expected << " | Actual: " << actual << std::endl;
                    return 1;
                }
            }
        }

        std::cout << "  - Standard C++23 LayoutLeft index stride verified with ZERO copy offsets." << std::endl;
        std::cout << "test_mdspan_strides PASS!" << std::endl;
        return 0;
    #else
        std::cerr << "test_mdspan_strides FAIL: Standard C++23 mdspan header missing!" << std::endl;
        return 1;
    #endif
}
