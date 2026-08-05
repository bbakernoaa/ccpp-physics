#ifndef CNVC90_HPP
#define CNVC90_HPP

#include "satmedmf_types.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

/**
 * @file cnvc90.hpp
 * @brief High-performance C++23 std::mdspan implementation of the GFS Convective Cloud Diagnostics scheme (cnvc90.f).
 *
 * This module computes the convective cloud fraction as well as the pressure 
 * at the base and top of the convective cloud between radiation calls.
 */

namespace satmedmf {
namespace interstitials {

/**
 * @brief Computes convective cloud fraction and cloud base/top pressures.
 *
 * Translates `cnvc90_run` into C++23. It processes:
 * -# **Reset / Zero-out Step**: Resets diagnostics if `clstp >= 1000.0` (start of accum interval).
 * -# **Accumulation Step**: If convective rain is falling, accumulates rainfall amounts and updates maximum cloud top and minimum cloud base heights.
 * -# **Normalization / Interpolation Step**: If within the radiation time limit, converts layer height indices to pressures from `prsi`,
 *    calculates a scaled convective rain proxy, and linearly interpolates the convective cloud fraction from empirical tables.
 *
 * @param clstp Convective accumulation time step interval counter.
 * @param im Number of horizontal grid columns.
 * @param km Number of vertical model layers.
 * @param rn Convective rainfall amount profile [columns].
 * @param kbot Convective cloud base layer indices [columns].
 * @param ktop Convective cloud top layer indices [columns].
 * @param prsi Pressure values at vertical layer interfaces [columns, layers + 1].
 * @param acv Accumulated convective rainfall amount [columns] (in-out).
 * @param acvb Accumulated minimum convective cloud base level indices [columns] (in-out).
 * @param acvt Accumulated maximum convective cloud top level indices [columns] (in-out).
 * @param cv Output convective cloud area fraction [columns].
 * @param cvb Output pressure at convective cloud base [columns].
 * @param cvt Output pressure at convective cloud top [columns].
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void cnvc90_run(
    double clstp, size_t im, size_t km,
    const double* rn, const int* kbot, const int* ktop,
    ConstView2D prsi,
    double* acv, double* acvb, double* acvt,
    double* cv, double* cvb, double* cvt,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;

    const double cons_100 = 100.0;
    const double cvb0 = 100.0;
    const size_t ncc = 9;

    double cc[ncc] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    double p[ncc]  = { 0.14, 0.31, 0.70, 1.6, 3.4, 7.7, 17.0, 38.0, 85.0 };

    int lz = 0;
    int lc = 0;

    if (clstp >= 1000.0) {
        lz = 1;
    }
    if (clstp >= 1100.0 || (clstp < 1000.0 && clstp >= 100.0)) {
        lc = 1;
    }

    double ah = std::fmod(clstp, cons_100);

    if (lz != 0) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            acv[i]  = 0.0;
            acvb[i] = cvb0;
            acvt[i] = 0.0;
        }
    }

    if (lc != 0) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (rn[i] > 0.0) {
                acv[i]  += rn[i];
                acvb[i]  = std::min(acvb[i], static_cast<double>(kbot[i]));
                acvt[i]  = std::max(acvt[i], static_cast<double>(ktop[i]));
            }
        }
    }

    if (ah > 0.01 && ah < 99.99) {
        std::vector<double> pmd(im, 0.0);
        std::vector<int> nmd(im, 0);

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (acv[i] > 0.0) {
                // Convert floating-point accumulated level indices to standard rounded layer bounds
                int itop = std::round(acvt[i]);
                int ibot = std::round(acvb[i]);
                
                // Map to interface pressure levels
                cvt[i] = prsi[i, itop + 1];
                cvb[i] = prsi[i, ibot];
            } else {
                cvb[i] = 0.0;
                cvt[i] = 0.0;
            }
            pmd[i] = acv[i] * (24000.0 / ah);
        }

        // Empirical range matching sweeps
        for (size_t n = 0; n < ncc; ++n) {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                if (pmd[i] > p[n]) {
                    nmd[i] = static_cast<int>(n + 1); // 1-based index emulation
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (nmd[i] == 0) {
                cv[i]  = 0.0;
                cvb[i] = 0.0;
                cvt[i] = 0.0;
            } else if (nmd[i] == static_cast<int>(ncc)) {
                cv[i]  = cc[ncc - 1];
            } else {
                size_t idx = static_cast<size_t>(nmd[i] - 1); // convert to 0-based
                double cc1 = cc[idx];
                double cc2 = cc[idx + 1];
                double p1  = p[idx];
                double p2  = p[idx + 1];
                cv[i] = cc1 + (cc2 - cc1) * (pmd[i] - p1) / (p2 - p1);
            }
        }
    }
}

} // namespace interstitials
} // namespace satmedmf

#endif // CNVC90_HPP
