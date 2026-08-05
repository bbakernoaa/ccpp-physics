#include "sfc_diff.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * @file sfc_diff.cpp
 * @brief High-performance C++23 std::mdspan implementation of GFS Surface Diffusion similarity solver (sfc_diff.f).
 */

namespace sfc {
namespace diff {

// 1. Aerodynamic roughness length formulations over water
// Option 6: Wang (2017) momentum roughness length
inline double znot_m_v6(double uref) {
    const double p13 = -1.296521881682694e-02;
    const double p12 =  2.855780863283819e-01;
    const double p11 = -1.597898515251717e+00;
    const double p10 = -8.396975715683501e+00;

    const double p25 =  3.790846746036765e-10;
    const double p24 =  3.281964357650687e-09;
    const double p23 =  1.962282433562894e-07;
    const double p22 = -1.240239171056262e-06;
    const double p21 =  1.739759082358234e-07;
    const double p20 =  2.147264020369413e-05;

    const double p35 =  1.840430200185075e-07;
    const double p34 = -2.793849676757154e-05;
    const double p33 =  1.735308193700643e-03;
    const double p32 = -6.139315534216305e-02;
    const double p31 =  1.255457892775006e+00;
    const double p30 = -1.663993561652530e+01;

    const double p40 =  4.579369142033410e-04;

    if (uref >= 0.0 && uref <= 6.5) {
        return std::exp(p10 + uref * (p11 + uref * (p12 + uref * p13)));
    } else if (uref > 6.5 && uref <= 15.7) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 15.7 && uref <= 53.0) {
        return std::exp(p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35)))));
    } else {
        return p40;
    }
}

// Option 6: Wang (2017) scalar roughness length for heat
inline double znot_t_v6(double uref) {
    const double p00 =  1.100000000000000e-04;
    const double p15 = -9.144581627678278e-10;
    const double p14 =  7.020346616456421e-08;
    const double p13 = -2.155602086883837e-06;
    const double p12 =  3.333848806567684e-05;
    const double p11 = -2.628501274963990e-04;
    const double p10 =  8.634221567969181e-04;

    const double p25 = -8.654513012535990e-12;
    const double p24 =  1.232380050058077e-09;
    const double p23 = -6.837922749505057e-08;
    const double p22 =  1.871407733439947e-06;
    const double p21 = -2.552246987137160e-05;
    const double p20 =  1.428968311457630e-04;

    const double p35 =  3.207515102100162e-12;
    const double p34 = -2.945761895342535e-10;
    const double p33 =  8.788972147364181e-09;
    const double p32 = -3.814457439412957e-08;
    const double p31 = -2.448983648874671e-06;
    const double p30 =  3.436721779020359e-05;

    const double p45 = -3.530687797132211e-11;
    const double p44 =  3.939867958963747e-09;
    const double p43 = -1.227668406985956e-08;
    const double p42 = -1.367469811838390e-05;
    const double p41 =  5.988240863928883e-04;
    const double p40 = -7.746288511324971e-03;

    const double p56 = -1.187982453329086e-13;
    const double p55 =  4.801984186231693e-11;
    const double p54 = -8.049200462388188e-09;
    const double p53 =  7.169872601310186e-07;
    const double p52 = -3.581694433758150e-05;
    const double p51 =  9.503919224192534e-04;
    const double p50 = -1.036679430885215e-02;

    const double p60 =  4.751256171799112e-05;

    if (uref >= 0.0 && uref < 5.9) {
        return p00;
    } else if (uref >= 5.9 && uref <= 15.4) {
        return p10 + uref * (p11 + uref * (p12 + uref * (p13 + uref * (p14 + uref * p15))));
    } else if (uref > 15.4 && uref <= 21.6) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 21.6 && uref <= 42.2) {
        return p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35))));
    } else if (uref > 42.2 && uref <= 53.3) {
        return p40 + uref * (p41 + uref * (p42 + uref * (p43 + uref * (p44 + uref * p45))));
    } else if (uref > 53.3 && uref <= 80.0) {
        return p50 + uref * (p51 + uref * (p52 + uref * (p53 + uref * (p54 + uref * (p55 + uref * p56)))));
    } else {
        return p60;
    }
}

// Option 7: Wang (2018) momentum roughness length
inline double znot_m_v7(double uref) {
    const double p13 = -1.296521881682694e-02;
    const double p12 =  2.855780863283819e-01;
    const double p11 = -1.597898515251717e+00;
    const double p10 = -8.396975715683501e+00;

    const double p25 =  3.790846746036765e-10;
    const double p24 =  3.281964357650687e-09;
    const double p23 =  1.962282433562894e-07;
    const double p22 = -1.240239171056262e-06;
    const double p21 =  1.739759082358234e-07;
    const double p20 =  2.147264020369413e-05;

    const double p35 =  1.897534489606422e-07;
    const double p34 = -3.019495980684978e-05;
    const double p33 =  1.931392924987349e-03;
    const double p32 = -6.797293095862357e-02;
    const double p31 =  1.346757797103756e+00;
    const double p30 = -1.707846930193362e+01;

    const double p40 =  3.371427455376717e-04;

    if (uref >= 0.0 && uref <= 6.5) {
        return std::exp(p10 + uref * (p11 + uref * (p12 + uref * p13)));
    } else if (uref > 6.5 && uref <= 15.7) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 15.7 && uref <= 53.0) {
        return std::exp(p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35)))));
    } else {
        return p40;
    }
}

// Option 7: Wang (2018) scalar roughness length for heat
inline double znot_t_v7(double uref) {
    const double p00 =  1.100000000000000e-04;
    const double p15 = -9.193764479895316e-10;
    const double p14 =  7.052217518653943e-08;
    const double p13 = -2.163419217747114e-06;
    const double p12 =  3.342963077911962e-05;
    const double p11 = -2.633566691328004e-04;
    const double p10 =  8.644979973037803e-04;

    const double p25 = -9.402722450219142e-12;
    const double p24 =  1.325396583616614e-09;
    const double p23 = -7.299148051141852e-08;
    const double p22 =  1.982901461144764e-06;
    const double p21 = -2.680293455916390e-05;
    const double p20 =  1.484341646128200e-04;

    const double p35 =  7.921446674311864e-12;
    const double p34 = -1.019028029546602e-09;
    const double p33 =  5.251986927351103e-08;
    const double p32 = -1.337841892062716e-06;
    const double p31 =  1.659454106237737e-05;
    const double p30 = -7.558911792344770e-05;

    const double p45 = -2.694370426850801e-10;
    const double p44 =  5.817362913967911e-08;
    const double p43 = -5.000813324746342e-06;
    const double p42 =  2.143803523428029e-04;
    const double p41 = -4.588070983722060e-03;
    const double p40 =  3.924356617245624e-02;

    const double p56 = -1.663918773476178e-13;
    const double p55 =  6.724854483077447e-11;
    const double p54 = -1.127030176632823e-08;
    const double p53 =  1.003683177025925e-06;
    const double p52 = -5.012618091180904e-05;
    const double p51 =  1.329762020689302e-03;
    const double p50 = -1.450062148367566e-02;
    const double p60 =  6.840803042788488e-05;

    if (uref >= 0.0 && uref < 5.9) {
        return p00;
    } else if (uref >= 5.9 && uref <= 15.4) {
        return p10 + uref * (p11 + uref * (p12 + uref * (p13 + uref * (p14 + uref * p15))));
    } else if (uref > 15.4 && uref <= 21.6) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 21.6 && uref <= 42.6) {
        return p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35))));
    } else if (uref > 42.6 && uref <= 53.0) {
        return p40 + uref * (p41 + uref * (p42 + uref * (p43 + uref * (p44 + uref * p45))));
    } else if (uref > 53.0 && uref <= 80.0) {
        return p50 + uref * (p51 + uref * (p52 + uref * (p53 + uref * (p54 + uref * (p55 + uref * p56)))));
    } else {
        return p60;
    }
}

// 2. High-fidelity GFS stability similarity solver (stability subroutine)
inline void stability_exact(
    double z1, double zvfun, double gdx, double tv1, double thv1, double wind, double z0max, double ztmax, double tvs, double grav,
    bool thsfc_loc,
    double& rb, double& fm, double& fh, double& fm10, double& fh2, double& cm, double& ch, double& stress, double& ustar
) {
    const double alpha = 5.0;
    const double a0 = -3.975;
    const double a1 = 12.32;
    const double alpha4 = 20.0;
    const double b1 = -7.755;
    const double b2 = 6.041;
    const double xkrefsqr = 0.3;
    const double xkmin = 0.05;
    const double xkgdx = 3000.0;
    const double a0p = -7.941;
    const double a1p = 24.75;
    const double b1p = -8.705;
    const double b2p = 7.899;
    const double zolmin = -10.0;
    const double zero = 0.0;
    const double one = 1.0;

    double z1i = one / z1;

    double xkzo;
    if (gdx >= xkgdx) {
        xkzo = one;
    } else {
        xkzo = gdx / xkgdx;
    }

    double tem1 = tv1 - tvs;
    if (tem1 > zero) {
        double tem2 = xkzo * zvfun;
        xkzo = std::min(std::max(tem2, xkmin), xkzo);
    }

    double zolmax = xkrefsqr / std::sqrt(xkzo);

    double dtv = thv1 - tvs;
    double adtv = std::max(std::abs(dtv), 0.001);
    dtv = (dtv >= 0.0 ? 1.0 : -1.0) * adtv;

    if (thsfc_loc) {
        rb = std::max(-5000.0, (grav + grav) * dtv * z1 / ((thv1 + tvs) * wind * wind));
    } else {
        rb = std::max(-5000.0, grav * dtv * z1 / (tv1 * wind * wind));
    }

    double temp1 = one / z0max;
    double temp2 = one / ztmax;
    fm = std::log((z0max + z1) * temp1);
    fh = std::log((ztmax + z1) * temp2);
    fm10 = std::log((z0max + 10.0) * temp1);
    fh2 = std::log((ztmax + 2.0) * temp2);
    double hlinf = rb * fm * fm / fh;
    hlinf = std::min(std::max(hlinf, zolmin), zolmax);

    double pm = 0.0, ph = 0.0, pm10 = 0.0, ph2 = 0.0;

    if (dtv >= zero) {
        double hl1 = hlinf;
        if (hlinf > 0.25) {
            double tem1_s = hlinf * z1i;
            double hl0inf = z0max * tem1_s;
            double hltinf = ztmax * tem1_s;
            double aa = std::sqrt(one + alpha4 * hlinf);
            double aa0 = std::sqrt(one + alpha4 * hl0inf);
            double bb = aa;
            double bb0 = std::sqrt(one + alpha4 * hltinf);
            pm = aa0 - aa + std::log((aa + one) / (aa0 + one));
            ph = bb0 - bb + std::log((bb + one) / (bb0 + one));
            double fms = fm - pm;
            double fhs = fh - ph;
            hl1 = fms * fms * rb / fhs;
            hl1 = std::min(hl1, zolmax);
        }

        double tem1_s = hl1 * z1i;
        double hl0 = z0max * tem1_s;
        double hlt = ztmax * tem1_s;
        double aa = std::sqrt(one + alpha4 * hl1);
        double aa0 = std::sqrt(one + alpha4 * hl0);
        double bb = aa;
        double bb0 = std::sqrt(one + alpha4 * hlt);
        pm = aa0 - aa + std::log((one + aa) / (one + aa0));
        ph = bb0 - bb + std::log((one + bb) / (one + bb0));

        double hl110 = hl1 * 10.0 * z1i;
        aa = std::sqrt(one + alpha4 * hl110);
        pm10 = aa0 - aa + std::log((one + aa) / (one + aa0));

        double hl12 = (hl1 + hl1) * z1i;
        bb = std::sqrt(one + alpha4 * hl12);
        ph2 = bb0 - bb + std::log((one + bb) / (one + bb0));
    } else {
        double olinf = z1 / hlinf;
        double tem1_s = 50.0 * z0max;
        if (std::abs(olinf) <= tem1_s) {
            hlinf = -z1 / tem1_s;
            hlinf = std::max(hlinf, zolmin);
        }

        if (hlinf >= -0.5) {
            double hl1 = hlinf;
            pm = (a0 + a1 * hl1) * hl1 / (one + (b1 + b2 * hl1) * hl1);
            ph = (a0p + a1p * hl1) * hl1 / (one + (b1p + b2p * hl1) * hl1);
            double hl110 = hl1 * 10.0 * z1i;
            pm10 = (a0 + a1 * hl110) * hl110 / (one + (b1 + b2 * hl110) * hl110);
            double hl12 = (hl1 + hl1) * z1i;
            ph2 = (a0p + a1p * hl12) * hl12 / (one + (b1p + b2p * hl12) * hl12);
        } else {
            double hl1 = -hlinf;
            double tem1_s = one / std::sqrt(hl1);
            pm = std::log(hl1) + 2.0 * std::sqrt(tem1_s) - 0.8776;
            ph = std::log(hl1) + 0.5 * tem1_s + 1.386;

            double hl110 = hl1 * 10.0 * z1i;
            pm10 = std::log(hl110) + 2.0 / std::sqrt(std::sqrt(hl110)) - 0.8776;

            double hl12 = (hl1 + hl1) * z1i;
            ph2 = std::log(hl12) + 0.5 / std::sqrt(hl12) + 1.386;
        }
    }

    fm = fm - pm;
    fh = fh - ph;
    fm10 = fm10 - pm10;
    fh2 = fh2 - ph2;

    const double ca = constants::karman; // 0.4
    cm = ca * ca / (fm * fm);
    ch = ca * ca / (fm * fh);

    double temp_min = 0.00001 / z1;
    cm = std::max(cm, temp_min);
    ch = std::max(ch, temp_min);

    stress = cm * wind * wind;
    ustar = std::sqrt(stress);
}

/**
 * @brief Standalone physical C++23 solver for GFS Surface Diffusion exchange coefficients.
 */
void sfc_diff_run(
    size_t columns,
    SurfaceSounding sounding,
    const Real* z0,
    int sfc_z0_type,
    View1D cm,
    View1D ch,
    View1D ustar,
    View1D stress
) {
    // Thread-safe parallel execution across horizontal columns
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < columns; ++i) {
        // Initialize outputs
        cm[i] = 0.0;
        ch[i] = 0.0;
        ustar[i] = 0.0;
        stress[i] = 0.0;

        // 1. Compute wind speed and clamp to prevent division by zero
        double u_val = sounding.u1[i];
        double v_val = sounding.v1[i];
        double wind = std::sqrt(u_val * u_val + v_val * v_val);
        double wind_clamped = std::max(1e-4, wind);

        // 2. Compute dynamically updated roughness lengths (z0_val, zt_val) over ocean
        double z0_val = z0[i];
        double zt_val = z0[i] * 0.1; // Default heat roughness scale

        if (sfc_z0_type == 6) {
            z0_val = znot_m_v6(wind_clamped);
            zt_val = znot_t_v6(wind_clamped);
        } else if (sfc_z0_type == 7) {
            z0_val = znot_m_v7(wind_clamped);
            zt_val = znot_t_v7(wind_clamped);
        }

        // 3. Compute bulk stability coefficients using the exact GFS stability solver
        double rb_dummy = 0.0;
        double fm_dummy = 0.0;
        double fh_dummy = 0.0;
        double fm10_dummy = 0.0;
        double fh2_dummy = 0.0;
        double cm_val = 0.0;
        double ch_val = 0.0;
        double stress_val = 0.0;
        double ustar_val = 0.0;

        // Use potential temperature reference (thsfc_loc = false)
        stability_exact(
            sounding.z1[i], 1.0, 3000.0, sounding.t1[i], sounding.t1[i], wind_clamped, z0_val, zt_val, sounding.tskin[i], constants::grav,
            false,
            rb_dummy, fm_dummy, fh_dummy, fm10_dummy, fh2_dummy, cm_val, ch_val, stress_val, ustar_val
        );

        cm[i] = cm_val;
        ch[i] = ch_val;
        ustar[i] = ustar_val;

        // Compute surface wind stress using air density
        double air_density = sounding.ps[i] / (constants::rd * sounding.t1[i]);
        stress[i] = air_density * ustar[i] * ustar[i];
    }
}

} // namespace diff
} // namespace sfc

extern "C" {

/**
 * @brief Flat C ABI entry point for the C++23 GFS Surface Diffusion similarity solver.
 */
void c_sfc_diff_run(
    size_t columns,
    const double* u1, const double* v1, const double* t1, const double* q1,
    const double* z1, const double* ps, const double* tskin,
    const double* z0,
    int sfc_z0_type,
    double* cm, double* ch, double* ustar, double* stress
) {
    sfc::ConstView1D u1_view(u1, columns);
    sfc::ConstView1D v1_view(v1, columns);
    sfc::ConstView1D t1_view(t1, columns);
    sfc::ConstView1D q1_view(q1, columns);
    sfc::ConstView1D z1_view(z1, columns);
    sfc::ConstView1D ps_view(ps, columns);
    sfc::ConstView1D tskin_view(tskin, columns);

    sfc::View1D cm_view(cm, columns);
    sfc::View1D ch_view(ch, columns);
    sfc::View1D ustar_view(ustar, columns);
    sfc::View1D stress_view(stress, columns);

    sfc::SurfaceSounding sounding{
        u1_view,
        v1_view,
        t1_view,
        q1_view,
        z1_view,
        ps_view,
        tskin_view
    };

    sfc::diff::sfc_diff_run(
        columns,
        sounding,
        z0,
        sfc_z0_type,
        cm_view,
        ch_view,
        ustar_view,
        stress_view
    );
}

}
