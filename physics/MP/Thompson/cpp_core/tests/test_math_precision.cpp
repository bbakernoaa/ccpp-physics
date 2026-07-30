#include "thompson_math_utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "Running test_math_precision..." << std::endl;

    std::vector<double> test_inputs = {0.0, -0.01, -0.05, -0.1, -0.2, -0.5, -0.8, -1.0, -2.5, -5.0, -10.0};

    for (double x : test_inputs) {
        double expected = std::exp(x);
        
        // 1. Verify Padé [2,2] rational approximation (always positive, very robust)
        double actual_pade = thompson::fast_exp_pade(x);
        double diff_pade = std::abs(expected - actual_pade);
        std::cout << "  [Padé]    x = " << x << " | expected: " << expected << " | actual: " << actual_pade << " | diff: " << diff_pade << std::endl;
        
        // Assert unphysical negative is never returned
        if (actual_pade < 0.0) {
            std::cerr << "test_math_precision FAIL: Padé returned unphysical negative value!" << std::endl;
            return 1;
        }

        // 2. Verify Chebyshev Minimax Polynomial
        double actual_mini = thompson::fast_exp_minimax(x);
        double diff_mini = std::abs(expected - actual_mini);
        std::cout << "  [Minimax] x = " << x << " | expected: " << expected << " | actual: " << actual_mini << " | diff: " << diff_mini << std::endl;

        // Tolerance gates
        if (x >= -1.0) {
            if (diff_mini > 2e-3) {
                std::cerr << "test_math_precision FAIL: Minimax exceeds 2e-3 limit at " << x << std::endl;
                return 1;
            }
            if (diff_pade > 2e-3) {
                std::cerr << "test_math_precision FAIL: Padé exceeds 2e-3 limit at " << x << std::endl;
                return 1;
            }
        }
    }

    std::cout << "test_math_precision PASS!" << std::endl;
    return 0;
}
