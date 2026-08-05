#include "satmedmf_pbl_interstitials.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

int main() {
    std::cout << "Starting satmedmf_pbl_interstitials unit verification..." << std::endl;

    const size_t columns = 2;
    const size_t layers = 10;
    const size_t ntrac = 15;
    const size_t nvdiff = 12;

    std::vector<double> qgrs_data(columns * layers * ntrac, 0.0);
    std::vector<double> vdftra_data(columns * layers * nvdiff, 0.0);

    // Initialize qgrs with test values (increasing sequential concentrations)
    for (size_t n = 0; n < ntrac; ++n) {
        for (size_t k = 0; k < layers; ++k) {
            for (size_t i = 0; i < columns; ++i) {
                qgrs_data[i + columns * (k + layers * n)] = static_cast<double>(n + 1) * 0.001;
            }
        }
    }

    satmedmf::View3D qgrs(qgrs_data.data(), columns, layers, ntrac);
    satmedmf::View3D vdftra(vdftra_data.data(), columns, layers, nvdiff);

    // Setup active indices (using 1-based Fortran indices matching SCM)
    int ntqv = 1;
    int ntcw = 2;
    int ntiw = 3;
    int ntrw = 4;
    int ntsw = 5;
    int ntgl = 6;
    int ntlnc = 7;
    int ntinc = 8;
    int ntrnc = 9;
    int ntoz = 10;
    int ntwa = 11;
    int ntia = 12;
    int ntke = 13;
    int ntkev = 12;

    int imp_physics_thompson = 2;
    int imp_physics_wsm6 = 1;
    int imp_physics = imp_physics_thompson;

    std::string errmsg;
    int errflg = 0;
    int rtg_ozone_index = -1;

    // 1. Run GFS_PBL_generic_pre_run translation
    satmedmf::interstitials::pbl_generic_pre_run(
        columns, layers, nvdiff, ntrac,
        rtg_ozone_index,
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, 0, 0,
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, 0, false, 14, 1, // trans_aero = false
        0, 0, 0, 0, 0, 0, 0, 0,
        imp_physics, 0, imp_physics_thompson, imp_physics_wsm6,
        0, false, 0, false, 0, true, false, false, false, false,
        false, false, false,
        qgrs, vdftra,
        errmsg, errflg
    );

    assert(errflg == 0);
    assert(rtg_ozone_index == ntoz);

    // Verify preprocessing mapping correctness
    std::cout << "  - GFS_PBL_generic_pre_run output checks:" << std::endl;
    std::cout << "    vdftra[0, 0, 0] (qv expected: 0.001): " << vdftra[0, 0, 0] << std::endl;
    std::cout << "    vdftra[0, 0, 1] (qc expected: 0.002): " << vdftra[0, 0, 1] << std::endl;
    std::cout << "    vdftra[0, 0, 2] (qi expected: 0.003): " << vdftra[0, 0, 2] << std::endl;
    std::cout << "    vdftra[0, 0, 9] (oz expected: 0.010): " << vdftra[0, 0, 9] << std::endl;

    assert(std::abs(vdftra[0, 0, 0] - 0.001) < 1e-15);
    assert(std::abs(vdftra[0, 0, 1] - 0.002) < 1e-15);
    assert(std::abs(vdftra[0, 0, 2] - 0.003) < 1e-15);
    assert(std::abs(vdftra[0, 0, 9] - 0.010) < 1e-15);

    // 2. Setup GFS_PBL_generic_post_run arrays
    std::vector<double> ten_t_data(columns * layers, 0.0001);
    std::vector<double> ten_u_data(columns * layers, 0.0002);
    std::vector<double> ten_v_data(columns * layers, 0.0003);
    std::vector<double> ten_t_pbl_data(columns * layers, 0.0);
    std::vector<double> ten_q_pbl_data(columns * layers, 0.0);
    std::vector<double> ten_q_data(columns * layers * ntrac, 0.0);
    std::vector<double> gt0_data(columns * layers, 280.0);
    std::vector<double> gv0_data(columns * layers, 5.0);
    std::vector<double> gu0_data(columns * layers, 10.0);
    std::vector<double> gq0_data(columns * layers * ntrac, 0.001);
    std::vector<double> dudt_data(columns * layers, 0.0);
    std::vector<double> dvdt_data(columns * layers, 0.0);
    std::vector<double> dtdt_data(columns * layers, 0.0);
    std::vector<double> dqdt_data(columns * layers * ntrac, 0.0);

    satmedmf::View2D ten_t(ten_t_data.data(), columns, layers);
    satmedmf::View2D ten_u(ten_u_data.data(), columns, layers);
    satmedmf::View2D ten_v(ten_v_data.data(), columns, layers);
    satmedmf::View2D ten_t_pbl(ten_t_pbl_data.data(), columns, layers);
    satmedmf::View2D ten_q_pbl(ten_q_pbl_data.data(), columns, layers);
    satmedmf::View3D ten_q(ten_q_data.data(), columns, layers, ntrac);
    satmedmf::View2D gt0(gt0_data.data(), columns, layers);
    satmedmf::View2D gv0(gv0_data.data(), columns, layers);
    satmedmf::View2D gu0(gu0_data.data(), columns, layers);
    satmedmf::View3D gq0(gq0_data.data(), columns, layers, ntrac);
    satmedmf::View2D dudt(dudt_data.data(), columns, layers);
    satmedmf::View2D dvdt(dvdt_data.data(), columns, layers);
    satmedmf::View2D dtdt(dtdt_data.data(), columns, layers);
    satmedmf::View3D dqdt(dqdt_data.data(), columns, layers, ntrac);

    double dtf = 1.0;
    double dtp = 300.0;
    double rd = 287.05;
    double cp = 1004.6;
    double fvirt = 0.608;
    double hvap = 2.5e6;
    double huge = 9.9e36;

    std::vector<double> t1(columns, 280.0);
    std::vector<double> q1(columns, 0.005);
    std::vector<double> hflx(columns, 10.0);
    std::vector<double> oceanfrac(columns, 0.0);
    std::vector<double> prsl_data(columns * layers, 90000.0);
    std::vector<double> wind(columns, 5.0);
    std::vector<double> stress_wat(columns, 0.01);
    std::vector<double> hflx_wat(columns, 10.0);
    std::vector<double> evap_wat(columns, 0.0001);
    std::vector<double> ugrs1(columns, 5.0);
    std::vector<double> vgrs1(columns, 5.0);
    std::vector<double> hffac(columns, 1.0);

    satmedmf::View2D prsl(prsl_data.data(), columns, layers);

    bool wet[columns] = {false, false};
    bool dry[columns] = {true, true};
    bool icy[columns] = {false, false};

    // Run generic post-run translation under tend_opt_pbl = 1 (immediate apply)
    satmedmf::interstitials::pbl_generic_post_run(
        columns, layers, nvdiff, ntrac,
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, 0, 0,
        ntwa, ntia, ntgl, ntoz, ntke, ntkev, 0,
        1, false, 14, 1, 0, 0, 0, 0, 0, 0, 0, 0, // tend_opt_pbl = 1
        imp_physics, 0, imp_physics_thompson, imp_physics_wsm6,
        0, 0, false, false, 0, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false,
        false, false, false,
        vdftra, ten_t, ten_u, ten_v,
        ten_t_pbl, ten_q_pbl, ten_q,
        nullptr, nullptr, nullptr, nullptr, dtf, dtp,
        dudt, dvdt, dtdt, dqdt,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, satmedmf::ConstIntView2D(),
        0, 0, 0, 0,
        wet, dry, icy,
        nullptr, hffac.data(),
        rd, cp, fvirt, hvap,
        t1.data(), q1.data(), prsl, hflx.data(), oceanfrac.data(),
        wind.data(), stress_wat.data(), hflx_wat.data(), evap_wat.data(),
        ugrs1.data(), vgrs1.data(),
        nullptr, nullptr, nullptr, nullptr,
        false, nullptr, nullptr, nullptr, nullptr,
        gt0, gv0, gu0, gq0,
        1, huge,
        errmsg, errflg
    );

    assert(errflg == 0);

    // Verify postprocessing tendencies application
    std::cout << "  - GFS_PBL_generic_post_run output checks:" << std::endl;
    std::cout << "    ten_t_pbl[0, 0] (expected: 0.0001): " << ten_t_pbl[0, 0] << std::endl;
    std::cout << "    ten_q_pbl[0, 0] (expected: 0.001):  " << ten_q_pbl[0, 0] << std::endl;
    std::cout << "    gt0[0, 0] after apply (expected: 280.03): " << gt0[0, 0] << std::endl;

    assert(std::abs(ten_t_pbl[0, 0] - 0.0001) < 1e-15);
    assert(std::abs(ten_q_pbl[0, 0] - 0.001) < 1e-15);
    assert(std::abs(gt0[0, 0] - 280.03) < 1e-12);

    std::cout << "✓ ALL interstitial verification gates PASSED!" << std::endl;
    return 0;
}
