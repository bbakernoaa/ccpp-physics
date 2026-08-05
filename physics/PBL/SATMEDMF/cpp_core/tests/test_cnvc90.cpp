#include "cnvc90.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

int main() {
    std::cout << "Starting cnvc90 convective cloud cover diagnostic verification..." << std::endl;

    const size_t columns = 2;
    const size_t layers = 10;

    std::vector<double> rn = { 0.005, 0.0 }; // rain in col 0, no rain in col 1
    std::vector<int> kbot = { 2, 0 };
    std::vector<int> ktop = { 8, 0 };

    std::vector<double> prsi_data(columns * (layers + 1), 91000.0);
    for (size_t k = 0; k <= layers; ++k) {
        prsi_data[0 + columns * k] = 100000.0 - k * 5000.0;
        prsi_data[1 + columns * k] = 100000.0 - k * 5000.0;
    }

    satmedmf::View2D prsi(prsi_data.data(), columns, layers + 1);

    std::vector<double> acv(columns, 0.0);
    std::vector<double> acvb(columns, 100.0);
    std::vector<double> acvt(columns, 0.0);

    std::vector<double> cv(columns, 0.0);
    std::vector<double> cvb(columns, 0.0);
    std::vector<double> cvt(columns, 0.0);

    std::string errmsg;
    int errflg = 0;

    // 1. Reset call (clstp = 1000.0 -> LZ = 1)
    satmedmf::interstitials::cnvc90_run(
        1000.0, columns, layers,
        rn.data(), kbot.data(), ktop.data(),
        prsi,
        acv.data(), acvb.data(), acvt.data(),
        cv.data(), cvb.data(), cvt.data(),
        errmsg, errflg
    );
    assert(errflg == 0);
    assert(acvb[0] == 100.0);

    // 2. Accumulate call (clstp = 1100.0 -> LC = 1)
    satmedmf::interstitials::cnvc90_run(
        1100.0, columns, layers,
        rn.data(), kbot.data(), ktop.data(),
        prsi,
        acv.data(), acvb.data(), acvt.data(),
        cv.data(), cvb.data(), cvt.data(),
        errmsg, errflg
    );
    assert(errflg == 0);
    assert(acv[0] == 0.005);
    assert(acvb[0] == 2.0);
    assert(acvt[0] == 8.0);

    // 3. Normalization/Interpolation call (clstp = 1110.5 -> AH = 10.5, within range)
    satmedmf::interstitials::cnvc90_run(
        1110.5, columns, layers,
        rn.data(), kbot.data(), ktop.data(),
        prsi,
        acv.data(), acvb.data(), acvt.data(),
        cv.data(), cvb.data(), cvt.data(),
        errmsg, errflg
    );
    assert(errflg == 0);
    std::cout << "  - cvb[0] pressure (expected: 90000.0): " << cvb[0] << std::endl;
    std::cout << "  - cvt[0] pressure (expected: 55000.0): " << cvt[0] << std::endl;
    std::cout << "  - cv[0] cloud fraction (expected: ~0.540): " << cv[0] << std::endl;

    assert(std::abs(cvb[0] - 90000.0) < 1e-9);
    assert(std::abs(cvt[0] - 55000.0) < 1e-9);
    assert(cv[0] > 0.0);

    std::cout << "✓ GFS cnvc90 C++ translation PASSED!" << std::endl;
    return 0;
}
