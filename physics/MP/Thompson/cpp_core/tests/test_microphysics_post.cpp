#include "thompson_microphysics_post.hpp"
#include "gfs_mp_generic_post.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

int main() {
    std::cout << "Starting Thompson & GFS Generic Microphysics postprocessing unit verification..." << std::endl;

    const size_t columns = 2;
    const size_t layers = 10;
    const size_t ntrac = 15;

    // 1. Verify Thompson-specific Post Limiters
    std::vector<double> dtgrs_data(columns * layers, 1.5); // temperature tendency = 1.5 K/s
    std::vector<double> tgrs_data(columns * layers, 280.0);
    std::vector<double> prslk_data(columns * layers, 0.9); // Exner function

    satmedmf::View2D dtgrs(dtgrs_data.data(), columns, layers);
    satmedmf::ConstView2D tgrs(tgrs_data.data(), columns, layers);
    satmedmf::ConstView2D prslk(prslk_data.data(), columns, layers);

    std::string errmsg;
    int errflg = 0;

    // Limit at ttendlim = 1.0 K/s (potential temperature tendency)
    double ttendlim = 1.0;
    thompson::post::thompson_post_init(ttendlim, errmsg, errflg);
    assert(errflg == 0);

    thompson::post::thompson_post_run(
        columns, layers,
        dtgrs, tgrs, prslk,
        300.0, ttendlim, 1,
        errmsg, errflg
    );
    assert(errflg == 0);

    // Initial potential tendency: dtgrs / prslk = 1.5 / 0.9 = 1.666 K/s.
    // Limited to 1.0 K/s, so limited dtgrs = 1.0 * prslk = 1.0 * 0.9 = 0.9 K/s.
    std::cout << "  - Thompson post_run output check:" << std::endl;
    std::cout << "    dtgrs[0, 0] after limit (expected 0.9): " << dtgrs[0, 0] << std::endl;
    assert(std::abs(dtgrs[0, 0] - 0.9) < 1e-15);

    thompson::post::thompson_post_finalize(errmsg, errflg);
    assert(errflg == 0);

    // 2. Verify GFS Generic Microphysics Post (GFS_MP_generic_post)
    std::vector<double> ten_t_data(columns * layers, 0.0001);
    std::vector<double> ten_u_data(columns * layers, 0.0);
    std::vector<double> ten_v_data(columns * layers, 0.0);
    std::vector<double> ten_q_data(columns * layers * ntrac, 0.0);

    // Setup input state and tendencies
    satmedmf::View2D ten_t(ten_t_data.data(), columns, layers);
    satmedmf::View2D ten_u(ten_u_data.data(), columns, layers);
    satmedmf::View2D ten_v(ten_v_data.data(), columns, layers);
    satmedmf::View3D ten_q(ten_q_data.data(), columns, layers, ntrac);

    std::vector<double> gt0_data(columns * layers, 280.0);
    std::vector<double> gu0_data(columns * layers, 10.0);
    std::vector<double> gv0_data(columns * layers, 5.0);
    std::vector<double> gq0_data(columns * layers * ntrac, 0.001);

    satmedmf::View2D gt0(gt0_data.data(), columns, layers);
    satmedmf::View2D gu0(gu0_data.data(), columns, layers);
    satmedmf::View2D gv0(gv0_data.data(), columns, layers);
    satmedmf::View3D gq0(gq0_data.data(), columns, layers, ntrac);

    std::vector<double> rainc(columns, 0.0005);
    std::vector<double> rain1(columns, 0.001);
    std::vector<double> rann_data(columns * 1, 0.0); // nominal random numbers
    satmedmf::View2D rann(rann_data.data(), columns, 1);

    std::vector<double> xlat(columns, 45.0);
    std::vector<double> xlon(columns, -90.0);
    std::vector<double> tsfc(columns, 280.0);

    std::vector<double> prsl_data(columns * layers, 90000.0);
    std::vector<double> prsi_data(columns * (layers + 1), 91000.0);
    std::vector<double> phii_data(columns * (layers + 1), 100.0);
    satmedmf::View2D prsl(prsl_data.data(), columns, layers);
    satmedmf::View2D prsi(prsi_data.data(), columns, layers + 1);
    satmedmf::View2D phii(phii_data.data(), columns, layers + 1);

    std::vector<double> ice(columns, 0.0);
    std::vector<double> snow(columns, 0.0);
    std::vector<double> graupel(columns, 0.0);
    std::vector<double> rain0(columns, 0.001);
    std::vector<double> ice0(columns, 0.0);
    std::vector<double> snow0(columns, 0.0005);
    std::vector<double> graupel0(columns, 0.0);

    std::vector<double> del_data(columns * layers, 1000.0);
    std::vector<double> phil_data(columns * layers, 500.0);
    satmedmf::View2D del_view(del_data.data(), columns, layers);
    satmedmf::View2D phil_view(phil_data.data(), columns, layers);

    std::vector<int> htop(columns, 0);

    std::vector<double> rain_out(columns, 0.0);
    std::vector<double> domr_diag(columns, 0.0);
    std::vector<double> domzr_diag(columns, 0.0);
    std::vector<double> domip_diag(columns, 0.0);
    std::vector<double> doms_diag(columns, 0.0);
    std::vector<double> tprcp(columns, 0.0);
    std::vector<double> srflag(columns, 0.0);
    std::vector<double> sr(columns, 0.0);

    std::vector<double> cnvprcp(columns, 0.0);
    std::vector<double> totprcp(columns, 0.0);
    std::vector<double> totice(columns, 0.0);
    std::vector<double> totsnw(columns, 0.0);
    std::vector<double> totgrp(columns, 0.0);
    std::vector<double> toticeb(columns, 0.0);
    std::vector<double> totsnwb(columns, 0.0);
    std::vector<double> totgrpb(columns, 0.0);
    std::vector<double> pwat(columns, 0.0);
    std::vector<double> cnvprcpb(columns * 1, 0.0);
    std::vector<double> totprcpb(columns * 1, 0.0);

    std::vector<double> frzr(columns, 0.0);
    std::vector<double> frzrb(columns, 0.0);
    std::vector<double> frozr(columns, 0.0);
    std::vector<double> frozrb(columns, 0.0);
    std::vector<double> tsnowp(columns, 0.0);
    std::vector<double> tsnowpb(columns, 0.0);
    std::vector<double> rhonewsn1(columns, 200.0);

    std::vector<double> drain_cpl(columns, 0.0);
    std::vector<double> dsnow_cpl(columns, 0.0);
    std::vector<double> rain_cpl(columns, 0.0);
    std::vector<double> snow_cpl(columns, 0.0);
    std::vector<double> rainc_cpl(columns, 0.0);

    double dtf = 1.0;
    double frain = 1.0;
    double con_g = 9.80665;
    double rhowater = 1000.0;
    double rainmin = 1e-10;
    double dtp = 300.0;
    int imp_physics_thompson = 2;

    satmedmf::interstitials::gfs_mp_generic_post_run(
        columns, layers, 1, 1, 1, 0, 0, ntrac,
        imp_physics_thompson, 1, imp_physics_thompson, 3, 4, 5, 6,
        false, true, false, false, false, true, true, // cal_pre = false, lssav = true
        con_g, rhowater, rainmin, dtf, frain,
        rainc.data(), rain1.data(), rann, xlat.data(), xlon.data(),
        ten_t, ten_u, ten_v, ten_q,
        satmedmf::View2D(), satmedmf::View2D(), satmedmf::View2D(), satmedmf::View3D(),
        gt0, gu0, gv0, gq0,
        prsl, prsi, phii, tsfc.data(),
        ice.data(), snow.data(), graupel.data(),
        rain0.data(), ice0.data(), snow0.data(), graupel0.data(),
        del_view, phil_view, htop.data(), satmedmf::View2D(),
        0, 0, 0, 0, 0,
        273.15, rain_out.data(), domr_diag.data(), domzr_diag.data(), domip_diag.data(), doms_diag.data(),
        tprcp.data(), srflag.data(), sr.data(),
        cnvprcp.data(), totprcp.data(), totice.data(), totsnw.data(), totgrp.data(),
        toticeb.data(), totsnwb.data(), totgrpb.data(), pwat.data(),
        cnvprcpb.data(), totprcpb.data(), 1,
        frzr.data(), frzrb.data(), frozr.data(), frozrb.data(),
        tsnowp.data(), tsnowpb.data(), rhonewsn1.data(),
        drain_cpl.data(), dsnow_cpl.data(), rain_cpl.data(), snow_cpl.data(), rainc_cpl.data(),
        0, 1, 2,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, dtp,
        errmsg, errflg
    );

    assert(errflg == 0);

    // Verify diagnostic updates
    std::cout << "  - GFS_MP_generic_post_run output checks:" << std::endl;
    std::cout << "    gt0[0, 0] after apply (expected 280.03): " << gt0[0, 0] << std::endl;
    std::cout << "    tprcp[0] (expected 0.0015): " << tprcp[0] << std::endl;
    std::cout << "    totice[0] (expected 0.0): " << totice[0] << std::endl;

    assert(std::abs(gt0[0, 0] - 280.03) < 1e-12);
    assert(std::abs(tprcp[0] - 0.0015) < 1e-15);
    assert(std::abs(totice[0] - 0.0) < 1e-15);

    std::cout << "✓ ALL microphysics postprocessing verification gates PASSED!" << std::endl;
    return 0;
}
