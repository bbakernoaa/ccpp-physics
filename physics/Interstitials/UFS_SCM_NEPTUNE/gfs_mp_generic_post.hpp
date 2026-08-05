#ifndef GFS_MP_GENERIC_POST_HPP
#define GFS_MP_GENERIC_POST_HPP

#include "satmedmf_types.hpp"
#include "satmedmf_math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

/**
 * @file gfs_mp_generic_post.hpp
 * @brief High-performance C++23 std::mdspan GFS Microphysics generic postprocessing interstitial scheme.
 *
 * This file contains the complete postprocessing translation of GFS_MP_generic_post.F90, solving
 * microphysics tendencies application, rain/snow density conversions, and diagnostic sums.
 */

// Fortran mangled calpreciptype subroutine signature
extern "C" {
    void calpreciptype_(
        const int* kdt, const int* nrcm, const int* im, const int* ix, const int* lm, const int* lp1,
        const double* randomno, const double* xlat, const double* xlon,
        const double* gt0, const double* gq0, const double* prsl, const double* prsi,
        const double* prec, const double* phii, const double* tskin,
        double* domr, double* domzr, double* domip, double* doms
    );
}

namespace satmedmf {
namespace interstitials {

/**
 * @brief Post-microphysics generic diagnostic and state application routine.
 *
 * This routine translates `GFS_MP_generic_post_run` into C++23. It processes:
 * -# **Tendencies Application**: Immediately applies or sums temperature, wind, and tracer tendencies (`tend_opt_mp`).
 * -# **Precipitation Partitioning**: Automatically partitions precipitation into rain and snow classes based on surface temperature and microphysics species.
 * -# **Reflectivity Merging**: Merges subgrid convective cloud reflectivity with resolved Thompson/NSSL explicitly simulated reflectivity.
 * -# **Snow/Ice Density**: Evaluates falling frozen precipitation densities (`rhonewsn1`) based on temperature and snow/graupel ratios.
 * -# **Coupled Fluxes**: Accumulates surface coupled fluxes to land/ocean boundary arrays if `cplflx` is enabled.
 * -# **Diagnostic Buckets**: Aggregates output precipitation diagnostics to global model output arrays.
 *
 * @param im Number of horizontal grid columns.
 * @param levs Number of vertical model layers.
 * @param kdt Timestep integration step counter.
 * @param tend_opt_mp Control flag identifying how microphysics tendencies are applied to state fields (1-4).
 * @param nrcm Total convective random number columns.
 * @param nncl Number of cloud species.
 * @param ntcw Cloud water start index.
 * @param ntrac Total number of model tracers.
 * @param imp_physics Active microphysics scheme parameter.
 * @param cal_pre Active dominant precipitation type calculation flag.
 * @param cplflx-cpllnd Coupling flags (Ocean, Chemistry, Land).
 * @param progsigma Use prognostic updraft area fraction.
 * @param con_g Gravitational constant.
 * @param rhowater Density of liquid water.
 * @param rainmin Minimum threshold for precipitation detection.
 * @param dtf Time step fractional scaling factor (s).
 * @param frain Time step explicit precip scaling factor.
 * @param rainc Convective precipitation amount.
 * @param rain1 Explicit precipitation amount.
 * @param rann Random number array for precip type calculations.
 * @param xlat-xlon Grid column latitudes and longitudes.
 * @param ten_t-ten_v Computed raw microphysical tendencies.
 * @param dudt-dqdt Accumulated state tendency update arrays.
 * @param gt0-gq0 Applied real-time state arrays (temperature, wind, tracers).
 * @param prsl-phil Coordinate layer and interface pressure and geopotential arrays.
 * @param tsfc Skin/surface temperature profile (K).
 * @param ice-graupel Surface precipitation species diagnostic arrays.
 * @param rain Output explicit surface rainfall diagnostic profile.
 * @param domr_diag-doms_diag Output precipitation type diagnosed sums.
 * @param tprcp Output total precipitation amount.
 * @param srflag Output snow/rain fraction (0 for all rain, 1 for all snow, or ratio).
 * @param totprcp-totgrp Surface accumulators for model output.
 * @param dtp Primary physics timestep (s).
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void gfs_mp_generic_post_run(
    size_t im, size_t levs, int kdt, int tend_opt_mp, int nrcm, int nncl, int ntcw, int ntrac,
    int imp_physics, int imp_physics_gfdl, int imp_physics_thompson, int imp_physics_tempo,
    int imp_physics_nssl, int imp_physics_mg, int imp_physics_fer_hires,
    bool cal_pre, bool cplflx, bool cplchm, bool cpllnd, bool progsigma, bool exticeden, bool lssav,
    double con_g, double rhowater, double rainmin, double dtf, double frain,
    double* rainc, const double* rain1, ConstView2D rann, const double* xlat, const double* xlon,
    ConstView2D ten_t, ConstView2D ten_u, ConstView2D ten_v, ConstView3D ten_q,
    View2D dudt, View2D dvdt, View2D dtdt, View3D dqdt,
    View2D gt0, View2D gu0, View2D gv0, View3D gq0,
    ConstView2D prsl, ConstView2D prsi, ConstView2D phii, const double* tsfc,
    double* ice, double* snow, double* graupel,
    const double* rain0, const double* ice0, const double* snow0, const double* graupel0,
    ConstView2D del, ConstView2D phil, const int* htop, View2D refl_10cm,
    int imfshalcnv, int imfshalcnv_gf, int imfdeepcnv, int imfdeepcnv_gf, int imfdeepcnv_samf,
    double con_t0c, double* rain, double* domr_diag, double* domzr_diag, double* domip_diag, double* doms_diag,
    double* tprcp, double* srflag, const double* sr,
    double* cnvprcp, double* totprcp, double* totice, double* totsnw, double* totgrp,
    double* toticeb, double* totsnwb, double* totgrpb, double* pwat,
    double* cnvprcpb, double* totprcpb, int num_diag_buckets,
    double* frzr, double* frzrb, double* frozr, double* frozrb,
    double* tsnowp, double* tsnowpb, double* rhonewsn1,
    double* drain_cpl, double* dsnow_cpl, double* rain_cpl, double* snow_cpl, double* rainc_cpl,
    int lsm, int lsm_ruc, int lsm_noahmp,
    double* raincprv, double* rainncprv, double* iceprv, double* snowprv, double* graupelprv,
    double* draincprv, double* drainncprv, double* diceprv, double* dsnowprv, double* dgraupelprv,
    double* dqdt_qmicro, double* prevsq, double dtp,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;

    const double con_p001 = 0.001;
    const double con_day = 86400.0;
    const double p850 = 85000.0;
    const double zero = 0.0;
    const double one = 1.0;
    const double qmin = 1.0e-8;
    const double dbzmin = -20.0;

    auto f2c = [](int idx) { return idx - 1; };

    double onebg = one / con_g;

    // Save initial temperature state (needed if overwritten by radar DA)
    std::vector<double> save_t_data(im * levs);
    View2D save_t(save_t_data.data(), im, levs);
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < levs; ++k) {
        for (size_t i = 0; i < im; ++i) {
            save_t[i, k] = gt0[i, k];
        }
    }

    // 1. Tendency Application Controls
    if (tend_opt_mp == 1) { // immediately apply
        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                gt0[i, k] += dtp * ten_t[i, k];
                gu0[i, k] += dtp * ten_u[i, k];
                gv0[i, k] += dtp * ten_v[i, k];
                for (size_t n = 0; n < ntrac; ++n) {
                    gq0[i, k, n] += dtp * ten_q[i, k, n];
                }
            }
        }
    } else if (tend_opt_mp == 2) { // add to sum
        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                dtdt[i, k] += ten_t[i, k];
                dudt[i, k] += ten_u[i, k];
                dvdt[i, k] += ten_v[i, k];
                for (size_t n = 0; n < ntrac; ++n) {
                    dqdt[i, k, n] += ten_q[i, k, n];
                }
            }
        }
    } else if (tend_opt_mp == 3) { // sum and apply
        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                gt0[i, k] += dtp * (dtdt[i, k] + ten_t[i, k]);
                dtdt[i, k] = 0.0;
                gu0[i, k] += dtp * (dudt[i, k] + ten_u[i, k]);
                dudt[i, k] = 0.0;
                gv0[i, k] += dtp * (dvdt[i, k] + ten_v[i, k]);
                dvdt[i, k] = 0.0;
                for (size_t n = 0; n < ntrac; ++n) {
                    gq0[i, k, n] += dtp * (dqdt[i, k, n] + ten_q[i, k, n]);
                    dqdt[i, k, n] = 0.0;
                }
            }
        }
    } else if (tend_opt_mp == 4) {
        // Current state unchanged
    } else {
        errflg = 1;
        errmsg = "A tendency application control was outside of the acceptable range (1-4)";
        return;
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        rain[i] = rainc[i] + frain * rain1[i];
    }

    // 2. Combine Convective & Microphysics Radar Reflectivity (logarithmic sum)
    if ((imp_physics == imp_physics_thompson || imp_physics == imp_physics_tempo || imp_physics == imp_physics_nssl) &&
        (imfdeepcnv == imfdeepcnv_samf || imfdeepcnv == imfdeepcnv_gf || imfshalcnv == imfshalcnv_gf)) {
        
        std::vector<double> zfrz(im);
        std::vector<double> factor(im, 0.0);
        std::vector<double> zo_data(im * levs);
        View2D zo(zo_data.data(), im, levs);

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            bool lfrz = true;
            zfrz[i] = phil[i, 0] * onebg;
            for (int k = static_cast<int>(levs) - 1; k >= 0; --k) {
                zo[i, k] = phil[i, k] * onebg;
                if (gt0[i, k] >= con_t0c && lfrz) {
                    zfrz[i] = zo[i, k];
                    lfrz = false;
                }
            }
            if (rainc[i] > 0.0 && htop[i] > 0) {
                factor[i] = -2.0 / std::max(1000.0, zo[i, f2c(htop[i])] - zfrz[i]);
            }
        }

        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                if (rainc[i] > 0.0 && static_cast<int>(k) <= f2c(htop[i])) {
                    double fctz = 0.0;
                    double delz = zo[i, k] - zfrz[i];
                    if (delz < 0.0) {
                        fctz = 1.0;
                    } else {
                        fctz = std::pow(10.0, factor[i] * delz);
                    }
                    double cuprate = rainc[i] * 3.6e6 / dtp; // rate mm/h
                    double ze_conv = 300.0 * std::pow(cuprate, 1.4);
                    ze_conv = fctz * ze_conv;
                    double ze_mp = std::pow(10.0, 0.1 * refl_10cm[i, k]);
                    double dbz_sum = std::max(dbzmin, 10.0 * std::log10(ze_mp + ze_conv));
                    refl_10cm[i, k] = dbz_sum;
                }
            }
        }
    }

    // 3. Falling frozen precipitation and density conversions
    if (imp_physics == imp_physics_gfdl || imp_physics == imp_physics_thompson ||
        imp_physics == imp_physics_tempo || imp_physics == imp_physics_nssl) {
        
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (gt0[i, 0] <= 273.15) {
                frzr[i] += rain0[i];
                frzrb[i] += rain0[i];
            }
            tsnowp[i] += snow0[i];
            tsnowpb[i] += snow0[i];
            frozr[i] += graupel0[i];
            frozrb[i] += graupel0[i];
        }

        if (exticeden) {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                rhonewsn1[i] = 200.0;
                double prcpncfr = rain1[i] * sr[i];
                double prcpcufr = 0.0;
                if (sr[i] > 0.0 && gt0[i, 0] < 273.15) {
                    prcpcufr = std::max(0.0, rainc[i] * sr[i]);
                } else if (gt0[i, 0] < 273.15) {
                    prcpcufr = std::max(0.0, rainc[i]);
                }

                double tot_frz = prcpncfr + prcpcufr;
                if (tot_frz > 0.0) {
                    double snowrat = std::min(1.0, std::max(0.0, snow0[i] / tot_frz));
                    double grauprat = std::min(1.0, std::max(0.0, graupel0[i] / tot_frz));
                    double icerat = std::min(1.0, std::max(0.0, (prcpncfr - snow0[i] - graupel0[i]) / tot_frz));
                    double curat = std::min(1.0, std::max(0.0, prcpcufr / tot_frz));

                    double rhonewsnow = std::min(125.0, 1000.0 / std::max(8.0, 17.0 * std::tanh((276.65 - gt0[i, 0]) * 0.15)));
                    double rhonewgr = std::min(500.0, rhowater / std::max(2.0, 3.5 * std::tanh((274.15 - gt0[i, 0]) * 0.3333)));
                    double rhonewice = rhonewsnow;

                    double rhoprcpice = std::min(500.0, std::max(58.8, rhonewsnow * snowrat + rhonewgr * grauprat +
                                                                       rhonewice * icerat + rhonewgr * curat));
                    rhonewsn1[i] = rhoprcpice;
                }
            }
        }
    }

    // 4. Surface Partitioning on physics timestep
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        if (imp_physics == imp_physics_gfdl) {
            tprcp[i] = std::max(zero, rain[i]);
            graupel[i] = graupel0[i];
            ice[i] = ice0[i];
            snow[i] = snow0[i];
        } else if (imp_physics == imp_physics_thompson || imp_physics == imp_physics_tempo || imp_physics == imp_physics_nssl) {
            tprcp[i] = std::max(zero, rainc[i] + frain * rain1[i]);
            graupel[i] = frain * graupel0[i];
            ice[i] = frain * ice0[i];
            snow[i] = frain * snow0[i];
        } else if (imp_physics == imp_physics_fer_hires) {
            tprcp[i] = std::max(zero, rain[i]);
            ice[i] = frain * rain1[i] * sr[i];
        }
    }

    // 5. LSM Previous steps mapping
    if (lsm == lsm_ruc || lsm == lsm_noahmp) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            raincprv[i] = rainc[i];
            rainncprv[i] = frain * rain1[i];
            iceprv[i] = ice[i];
            snowprv[i] = snow[i];
            graupelprv[i] = graupel[i];

            if (lsm == lsm_noahmp) {
                double tem_val = one / (dtp * con_p001);
                draincprv[i] = tem_val * raincprv[i];
                drainncprv[i] = tem_val * rainncprv[i];
                dsnowprv[i] = tem_val * snowprv[i];
                dgraupelprv[i] = tem_val * graupelprv[i];
                diceprv[i] = tem_val * iceprv[i];
            }
        }
    }

    // 6. Precipitation type calculation call
    if (cal_pre) {
        int im_i = static_cast<int>(im);
        int lm_i = static_cast<int>(levs);
        int lp1_i = lm_i + 1;
        std::vector<double> domr(im), domzr(im), domip(im), doms(im);

        // Fortran linkage using column-major, contiguous parameters
        calpreciptype_(&kdt, &nrcm, &im_i, &im_i, &lm_i, &lp1_i,
                       rann.data_handle(), xlat, xlon, gt0.data_handle(),
                       gq0.data_handle(), prsl.data_handle(), prsi.data_handle(),
                       rain, phii.data_handle(), tsfc,
                       domr.data(), domzr.data(), domip.data(), doms.data());

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (imp_physics != imp_physics_gfdl && imp_physics != imp_physics_thompson &&
                imp_physics != imp_physics_tempo && imp_physics != imp_physics_nssl) {
                tprcp[i] = std::max(zero, rain[i]);
                if (doms[i] > zero || domip[i] > zero) {
                    srflag[i] = one;
                } else {
                    srflag[i] = zero;
                }
            }
            if (lssav) {
                domr_diag[i] += domr[i] * dtf;
                domzr_diag[i] += domzr[i] * dtf;
                domip_diag[i] += domip[i] * dtf;
                doms_diag[i] += doms[i] * dtf;
            }
        }
    }

    // 7. Calculate 850 hPa temperature
    std::vector<double> t850(im);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        t850[i] = gt0[i, 0];
    }
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < levs - 1; ++k) {
        for (size_t i = 0; i < im; ++i) {
            if (prsl[i, k] > p850 && prsl[i, k + 1] <= p850) {
                t850[i] = gt0[i, k] - (prsl[i, k] - p850) / (prsl[i, k] - prsl[i, k + 1]) * (gt0[i, k] - gt0[i, k + 1]);
            }
        }
    }

    // 8. Determine explicit/convective rain snow fractions / srflag
    if (imp_physics == imp_physics_gfdl || imp_physics == imp_physics_thompson ||
        imp_physics == imp_physics_tempo || imp_physics == imp_physics_nssl) {
        if (lsm != lsm_ruc) {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                srflag[i] = zero;
                double crain = rainc[i];
                double csnow = zero;
                if (tsfc[i] < 273.15) {
                    crain = zero;
                    csnow = rainc[i];
                }
                double total_precip = snow0[i] + ice0[i] + graupel0[i] + rain0[i] + rainc[i];
                if (total_precip > rainmin) {
                    srflag[i] = (snow0[i] + ice0[i] + graupel0[i] + csnow) / total_precip;
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                srflag[i] = sr[i];
            }
        }
    } else if (!cal_pre) {
        if (imp_physics == imp_physics_mg) {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                if (rain[i] > rainmin) {
                    double tem1 = std::max(zero, rain[i] - rainc[i]) * sr[i];
                    double tem2 = one / rain[i];
                    if (t850[i] > 273.16) {
                        srflag[i] = std::max(zero, std::min(one, tem1 * tem2));
                    } else {
                        srflag[i] = std::max(zero, std::min(one, (tem1 + rainc[i]) * tem2));
                    }
                } else {
                    srflag[i] = zero;
                    rain[i] = zero;
                    rainc[i] = zero;
                }
                tprcp[i] = std::max(zero, rain[i]);
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                tprcp[i] = std::max(zero, rain[i]);
                srflag[i] = sr[i];
            }
        }
    }

    // 9. Saving diagnostic aggregations
    if (lssav) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            cnvprcp[i] += rainc[i];
            totprcp[i] += rain[i];
            totice[i]  += ice[i];
            totsnw[i]  += snow[i];
            totgrp[i]  += graupel[i];

            toticeb[i] += ice[i];
            totsnwb[i] += snow[i];
            totgrpb[i] += graupel[i];
        }
        for (int ib = 0; ib < num_diag_buckets; ++ib) {
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < im; ++i) {
                cnvprcpb[i * num_diag_buckets + ib] += rainc[i];
                totprcpb[i * num_diag_buckets + ib] += rain[i];
            }
        }
    }

    // 10. Prog updraft area fraction tendencies
    if (progsigma) {
        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                dqdt_qmicro[i * levs + k] = ten_q[i, k, 0];
            }
        }
    }

    // 11. Coupled Boundary fluxes mapping
    if (cplflx || cplchm || cpllnd) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            dsnow_cpl[i] = std::max(zero, rain[i] * srflag[i]);
            drain_cpl[i] = std::max(zero, rain[i] - dsnow_cpl[i]);
            rain_cpl[i] += drain_cpl[i];
            snow_cpl[i] += dsnow_cpl[i];
        }
    }

    if (cplchm || cpllnd) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            rainc_cpl[i] += rainc[i];
        }
    }

    // 12. Calculate precipitable water profile (pwat)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        pwat[i] = zero;
    }
    for (size_t k = 0; k < levs; ++k) {
        std::vector<double> work1(im, 0.0);
        if (nncl > 0) {
            for (int ic = ntcw; ic < ntcw + nncl; ++ic) {
                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < im; ++i) {
                    work1[i] += gq0[i, k, f2c(ic)];
                }
            }
        }
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            pwat[i] += del[i, k] * (gq0[i, k, 0] + work1[i]);
        }
    }
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < im; ++i) {
        pwat[i] *= onebg;
    }

    if (progsigma) {
        #pragma omp parallel for schedule(static)
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                prevsq[i * levs + k] = gq0[i, k, 0];
            }
        }
    }
}

} // namespace interstitials
} // namespace satmedmf

#endif // GFS_MP_GENERIC_POST_HPP
