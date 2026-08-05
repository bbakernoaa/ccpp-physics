#ifndef SATMEDMF_PBL_INTERSTITIALS_HPP
#define SATMEDMF_PBL_INTERSTITIALS_HPP

#include "satmedmf_types.hpp"
#include "satmedmf_math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

/**
 * @file satmedmf_pbl_interstitials.hpp
 * @brief High-performance C++23 std::mdspan implementation of the GFS PBL generic pre- and post-processing interstitial schemes.
 *
 * This file contains the complete physical mappings, microphysics index translations, and mass-conservative
 * boundary layer interstitial routines designed to bridge CCPP state vectors with the SATMEDMF solver.
 */

namespace satmedmf {
namespace interstitials {

/**
 * @brief Helper routine to retrieve the index offset in the diffused tracer array where aerosol/chemistry tracers begin.
 *
 * Matches the logic of set_aerosol_tracer_index from GFS_PBL_generic_common.F90.
 *
 * @param imp_physics Selected microphysics scheme option index.
 * @param imp_physics_wsm6 WSM6 microphysics index.
 * @param imp_physics_thompson Thompson microphysics index.
 * @param ltaerosol Active aerosol-aware microphysics flag.
 * @param mraerosol MERRA2-aerosol-aware microphysics flag.
 * @param imp_physics_tempo Tempo microphysics index.
 * @param lthailaware Active hail-aware microphysics flag.
 * @param imp_physics_mg Morrison-Gettelman (MG) microphysics index.
 * @param ntgl Number of graupel tracers.
 * @param imp_physics_gfdl GFDL microphysics index.
 * @param imp_physics_nssl NSSL microphysics index.
 * @param nssl_hail_on Active NSSL hail-aware flag.
 * @param nssl_ccn_on Active NSSL CCN-aware flag.
 * @param nssl_3moment Active NSSL 3-moment flag.
 * @param kk Output calculated tracer offset.
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void set_aerosol_tracer_index(
    int imp_physics, int imp_physics_wsm6,
    int imp_physics_thompson, bool ltaerosol, bool mraerosol,
    int imp_physics_tempo, bool lthailaware,
    int imp_physics_mg, int ntgl, int imp_physics_gfdl,
    int imp_physics_nssl,
    bool nssl_hail_on, bool nssl_ccn_on, bool nssl_3moment,
    int& kk,
    std::string& errmsg, int& errflg
) {
    errflg = 0;
    if (imp_physics == imp_physics_wsm6) {
        kk = 4;
    } else if (imp_physics == imp_physics_thompson) {
        if (ltaerosol) {
            kk = 12;
        } else if (mraerosol) {
            kk = 10;
        } else {
            kk = 9;
        }
    } else if (imp_physics == imp_physics_tempo) {
        kk = 9;
        if (ltaerosol) {
            kk += 3;
        }
        if (lthailaware) {
            kk += 2;
        }
    } else if (imp_physics == imp_physics_mg) {
        if (ntgl > 0) {
            kk = 12;
        } else {
            kk = 10;
        }
    } else if (imp_physics == imp_physics_gfdl) {
        kk = 7;
    } else if (imp_physics == imp_physics_nssl) {
        if (nssl_hail_on) {
            kk = 16;
            if (nssl_3moment) kk += 3;
        } else {
            kk = 13;
            if (nssl_3moment) kk += 2;
        }
        if (nssl_ccn_on) kk += 1;
    } else {
        errmsg = "Logic error: unknown microphysics option in set_aerosol_tracer_index";
        kk = -999;
        errflg = 1;
    }
}

/**
 * @brief Preprocessing interstitial scheme (CCPP-compliant GFS_PBL_generic_pre_run).
 *
 * Maps state tracer arrays to the vertically diffused tracer array `vdftra` before executing local PBL physics.
 * Mapping configurations are fully adapted based on the active microphysics scheme options.
 *
 * @param im Number of horizontal grid columns.
 * @param levs Number of vertical model layers.
 * @param nvdiff Total number of diffused tracers.
 * @param ntrac Total number of atmospheric state tracers.
 * @param rtg_ozone_index Output tracked index inside diffused tracer array where Ozone resides.
 * @param ntqv-nthz Subgrid indices identifying water vapor, hydrometeors, TKE, and chemical species.
 * @param imp_physics-mraerosol Microphysics flags and option parameters.
 * @param qgrs Input atmospheric state tracers multidimensional View [columns, layers, ntrac].
 * @param vdftra Output vertically diffused tracer View [columns, layers, nvdiff].
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void pbl_generic_pre_run(
    size_t im, size_t levs, size_t nvdiff, size_t ntrac,
    int& rtg_ozone_index,
    int ntqv, int ntcw, int ntiw, int ntrw, int ntsw, int ntlnc, int ntinc, int ntrnc, int ntsnc, int ntgnc,
    int ntwa, int ntia, int ntgl, int ntoz, int ntke, int ntkev, int nqrimef, bool trans_aero, int ntchs, int ntchm,
    int ntccn, int nthl, int nthnc, int ntgv, int nthv, int ntrz, int ntgz, int nthz,
    int imp_physics, int imp_physics_gfdl, int imp_physics_thompson, int imp_physics_wsm6,
    int imp_physics_tempo, bool lthailaware, int imp_physics_mg, bool imp_physics_fer_hires,
    int imp_physics_nssl, bool ltaerosol, bool mraerosol, bool nssl_ccn_on, bool nssl_hail_on, bool nssl_3moment,
    bool hybedmf, bool do_shoc, bool satmedmf,
    ConstView3D qgrs,
    View3D vdftra,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;
    rtg_ozone_index = -1;

    // Helper lambda to map 1-based Fortran indices safely to 0-based C++ coordinates
    auto f2c = [](int idx) { return idx - 1; };

    if (nvdiff == ntrac && (hybedmf || do_shoc || satmedmf)) {
        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < ntrac; ++n) {
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, n] = qgrs[i, k, n];
                }
            }
        }
        rtg_ozone_index = ntoz;
    } else {
        if (imp_physics == imp_physics_wsm6) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, 0] = qgrs[i, k, f2c(ntqv)];
                    vdftra[i, k, 1] = qgrs[i, k, f2c(ntcw)];
                    vdftra[i, k, 2] = qgrs[i, k, f2c(ntiw)];
                    vdftra[i, k, 3] = qgrs[i, k, f2c(ntoz)];
                }
            }
            rtg_ozone_index = 4;
        } else if (imp_physics == imp_physics_fer_hires) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, 0] = qgrs[i, k, f2c(ntqv)];
                    vdftra[i, k, 1] = qgrs[i, k, f2c(ntcw)];
                    vdftra[i, k, 2] = qgrs[i, k, f2c(ntiw)];
                    vdftra[i, k, 3] = qgrs[i, k, f2c(ntrw)];
                    vdftra[i, k, 4] = qgrs[i, k, f2c(nqrimef)];
                    vdftra[i, k, 5] = qgrs[i, k, f2c(ntoz)];
                }
            }
            rtg_ozone_index = 6;
        } else if (imp_physics == imp_physics_thompson) {
            if (ltaerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntoz)];
                        vdftra[i, k, 10] = qgrs[i, k, f2c(ntwa)];
                        vdftra[i, k, 11] = qgrs[i, k, f2c(ntia)];
                    }
                }
                rtg_ozone_index = 10;
            } else if (mraerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntoz)];
                    }
                }
                rtg_ozone_index = 10;
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0] = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1] = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2] = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3] = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4] = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5] = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6] = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 7] = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 8] = qgrs[i, k, f2c(ntoz)];
                    }
                }
                rtg_ozone_index = 9;
            }
        } else if (imp_physics == imp_physics_tempo) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, 0] = qgrs[i, k, f2c(ntqv)];
                    vdftra[i, k, 1] = qgrs[i, k, f2c(ntcw)];
                    vdftra[i, k, 2] = qgrs[i, k, f2c(ntiw)];
                    vdftra[i, k, 3] = qgrs[i, k, f2c(ntrw)];
                    vdftra[i, k, 4] = qgrs[i, k, f2c(ntsw)];
                    vdftra[i, k, 5] = qgrs[i, k, f2c(ntgl)];
                    vdftra[i, k, 6] = qgrs[i, k, f2c(ntinc)];
                    vdftra[i, k, 7] = qgrs[i, k, f2c(ntrnc)];
                    vdftra[i, k, 8] = qgrs[i, k, f2c(ntoz)];
                }
            }
            rtg_ozone_index = 9;
            int n = 9;
            if (ltaerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, n] = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntwa)];
                        vdftra[i, k, n + 2] = qgrs[i, k, f2c(ntia)];
                    }
                }
                n = 12;
            }
            if (lthailaware) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, n] = qgrs[i, k, f2c(ntgnc)];
                        vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntgv)];
                    }
                }
            }
        } else if (imp_physics == imp_physics_mg) {
            if (ntgl > 0) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntsnc)];
                        vdftra[i, k, 10] = qgrs[i, k, f2c(ntgnc)];
                        vdftra[i, k, 11] = qgrs[i, k, f2c(ntoz)];
                    }
                }
                rtg_ozone_index = 12;
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntsnc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntoz)];
                    }
                }
                rtg_ozone_index = 10;
            }
        } else if (imp_physics == imp_physics_gfdl) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, 0] = qgrs[i, k, f2c(ntqv)];
                    vdftra[i, k, 1] = qgrs[i, k, f2c(ntcw)];
                    vdftra[i, k, 2] = qgrs[i, k, f2c(ntiw)];
                    vdftra[i, k, 3] = qgrs[i, k, f2c(ntrw)];
                    vdftra[i, k, 4] = qgrs[i, k, f2c(ntsw)];
                    vdftra[i, k, 5] = qgrs[i, k, f2c(ntgl)];
                    vdftra[i, k, 6] = qgrs[i, k, f2c(ntoz)];
                }
            }
            rtg_ozone_index = 7;
        } else if (imp_physics == imp_physics_nssl) {
            if (nssl_hail_on) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(nthl)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 10] = qgrs[i, k, f2c(ntsnc)];
                        vdftra[i, k, 11] = qgrs[i, k, f2c(ntgnc)];
                        vdftra[i, k, 12] = qgrs[i, k, f2c(nthnc)];
                        vdftra[i, k, 13] = qgrs[i, k, f2c(ntgv)];
                        vdftra[i, k, 14] = qgrs[i, k, f2c(nthv)];
                        vdftra[i, k, 15] = qgrs[i, k, f2c(ntoz)];
                    }
                }
                int n = 15;
                if (nssl_ccn_on) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntccn)];
                        }
                    }
                    n += 1;
                }
                if (nssl_3moment) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntrz)];
                            vdftra[i, k, n + 2] = qgrs[i, k, f2c(ntgz)];
                            vdftra[i, k, n + 3] = qgrs[i, k, f2c(nthz)];
                        }
                    }
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, 0]  = qgrs[i, k, f2c(ntqv)];
                        vdftra[i, k, 1]  = qgrs[i, k, f2c(ntcw)];
                        vdftra[i, k, 2]  = qgrs[i, k, f2c(ntiw)];
                        vdftra[i, k, 3]  = qgrs[i, k, f2c(ntrw)];
                        vdftra[i, k, 4]  = qgrs[i, k, f2c(ntsw)];
                        vdftra[i, k, 5]  = qgrs[i, k, f2c(ntgl)];
                        vdftra[i, k, 6]  = qgrs[i, k, f2c(ntlnc)];
                        vdftra[i, k, 7]  = qgrs[i, k, f2c(ntinc)];
                        vdftra[i, k, 8]  = qgrs[i, k, f2c(ntrnc)];
                        vdftra[i, k, 9]  = qgrs[i, k, f2c(ntsnc)];
                        vdftra[i, k, 10] = qgrs[i, k, f2c(ntgnc)];
                        vdftra[i, k, 11] = qgrs[i, k, f2c(ntgv)];
                        vdftra[i, k, 12] = qgrs[i, k, f2c(ntoz)];
                    }
                }
                int n = 12;
                if (nssl_ccn_on) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntccn)];
                        }
                    }
                    n += 1;
                }
                if (nssl_3moment) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            vdftra[i, k, n + 1] = qgrs[i, k, f2c(ntrz)];
                            vdftra[i, k, n + 2] = qgrs[i, k, f2c(ntgz)];
                        }
                    }
                }
            }
        }

        int kk = 0;
        if (trans_aero) {
            set_aerosol_tracer_index(
                imp_physics, imp_physics_wsm6, imp_physics_thompson, ltaerosol, mraerosol,
                imp_physics_tempo, lthailaware, imp_physics_mg, ntgl, imp_physics_gfdl, imp_physics_nssl,
                nssl_hail_on, nssl_ccn_on, nssl_3moment, kk, errmsg, errflg
            );
            if (errflg != 0) return;

            int k1 = kk;
            for (int n = ntchs; n < ntchm + ntchs; ++n) {
                k1++;
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        vdftra[i, k, f2c(k1)] = qgrs[i, k, f2c(n)];
                    }
                }
            }
        }

        if (ntke > 0) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    vdftra[i, k, f2c(ntkev)] = qgrs[i, k, f2c(ntke)];
                }
            }
        }
    }
}

/**
 * @brief Postprocessing interstitial scheme (CCPP-compliant GFS_PBL_generic_post_run).
 *
 * Recovers diffused tracer increments, maps them back into standard atmospheric chemistry state,
 * and updates temperature, moisture, and momentum tendencies. Computes surface coupled ocean/tile fluxes.
 *
 * @param im Number of horizontal grid columns.
 * @param levs Number of vertical model layers.
 * @param nvdiff Total number of diffused tracers.
 * @param ntrac Total number of atmospheric state tracers.
 * @param ntqv-nthz Subgrid tracer indices.
 * @param tend_opt_pbl Control flag identifying how PBL tendencies are accumulated and applied to the state (options 1-4).
 * @param trans_aero Active aerosol-aware chemistry flag.
 * @param cplflx-cplaqm Direct coupled model coupling flags (Ocean, Chemistry, Air Quality).
 * @param dtf Time step fractional scaling factor (s).
 * @param dtp Primary physics integration step (s).
 * @param rd Gas constant for dry air.
 * @param cp Specific heat at constant pressure.
 * @param fvirt Virtual potential temperature expansion coefficient.
 * @param hvap Latent heat of vaporization.
 * @param t1-q1 Surface temperature and humidity values.
 * @param prsl Center-layer pressure View.
 * @param hflx Surface sensible heat flux.
 * @param oceanfrac Subgrid ocean fraction profile.
 * @param kdt Timestep integration step counter.
 * @param wind-evap_wat Boundary layer variables at ocean interfaces.
 * @param qgrs Input atmospheric state tracers View.
 * @param ugrs-tgrs Input atmospheric state wind components and temperature View.
 * @param dvdftra Input vertically diffused tracer concentrations View [columns, layers, nvdiff].
 * @param ten_t-ten_v Input computed raw PBL tendencies.
 * @param ten_t_pbl-ten_q_pbl Output final temperature and moisture tendencies.
 * @param ten_q Output updated tracer tendencies View [columns, layers, ntrac].
 * @param dusfc1-dqsfc1 Surface friction flux increments.
 * @param dudt-dtdt Accumulated and applied atmospheric tendencies.
 * @param dqdt Output accumulated and applied tracer tendencies View.
 * @param gt0-gv0 Real-time applied state variables (temperature, wind).
 * @param gq0 Output real-time applied tracer concentrations View.
 * @param errmsg Output diagnostic error message.
 * @param errflg Output diagnostic error flag (0 for success, 1 for failure).
 */
inline void pbl_generic_post_run(
    size_t im, size_t levs, size_t nvdiff, size_t ntrac,
    int ntqv, int ntcw, int ntiw, int ntrw, int ntsw, int ntlnc, int ntinc, int ntrnc, int ntsnc, int ntgnc,
    int ntwa, int ntia, int ntgl, int ntoz, int ntke, int ntkev, int nqrimef,
    int tend_opt_pbl, bool trans_aero, int ntchs, int ntchm, int ntccn, int nthl, int nthnc, int ntgv, int nthv,
    int ntrz, int ntgz, int nthz, int imp_physics, int imp_physics_gfdl, int imp_physics_thompson, int imp_physics_wsm6,
    int imp_physics_mg, int imp_physics_tempo, bool lthailaware, bool imp_physics_fer_hires, int imp_physics_nssl,
    bool nssl_ccn_on, bool ltaerosol, bool mraerosol, bool nssl_hail_on, bool nssl_3moment,
    bool cplflx, bool cplaqm, bool cplchm, bool lssav, bool flag_for_pbl_generic_tend, bool ldiag3d, bool lsidea,
    bool hybedmf, bool do_shoc, bool satmedmf, bool shinhong, bool do_ysu,
    ConstView3D dvdftra,
    ConstView2D ten_t, ConstView2D ten_u, ConstView2D ten_v,
    View2D ten_t_pbl, View2D ten_q_pbl,
    View3D ten_q,
    const double* dusfc1, const double* dvsfc1, const double* dtsfc1, const double* dqsfc1,
    double dtf, double dtp,
    View2D dudt, View2D dvdt, View2D dtdt,
    View3D dqdt,
    double* dusfc_cpl, double* dvsfc_cpl, double* dtsfc_cpl, double* dqsfc_cpl,
    double* dusfci_cpl, double* dvsfci_cpl, double* dtsfci_cpl, double* dqsfci_cpl,
    double* dusfc_diag, double* dvsfc_diag, double* dtsfc_diag, double* dqsfc_diag,
    double* dusfci_diag, double* dvsfci_diag, double* dtsfci_diag, double* dqsfci_diag,
    double* dtend_data, ConstIntView2D dtidx,
    int index_of_temperature, int index_of_x_wind, int index_of_y_wind, int index_of_process_pbl,
    const bool* wet, const bool* dry, const bool* icy,
    double* ushfsfci,
    const double* hffac,
    double rd, double cp, double fvirt, double hvap,
    const double* t1, const double* q1, ConstView2D prsl, const double* hflx, const double* oceanfrac,
    const double* wind, const double* stress_wat, const double* hflx_wat, const double* evap_wat,
    const double* ugrs1, const double* vgrs1,
    const double* dusfc_cice, const double* dvsfc_cice, const double* dtsfc_cice, const double* dqsfc_cice,
    bool use_med_flux,
    const double* dtsfc_med, const double* dqsfc_med, const double* dusfc_med, const double* dvsfc_med,
    View2D gt0, View2D gv0, View2D gu0, View3D gq0,
    int kdt, double huge,
    std::string& errmsg, int& errflg
) {
    errmsg = "";
    errflg = 0;
    const double qmin = 1.0e-8;

    // Helper lambda to map 1-based Fortran indices safely to 0-based C++ coordinates
    auto f2c = [](int idx) { return idx - 1; };

    #pragma omp parallel for schedule(static)
    for (size_t n = 0; n < ntrac; ++n) {
        for (size_t k = 0; k < levs; ++k) {
            for (size_t i = 0; i < im; ++i) {
                ten_q[i, k, n] = 0.0;
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < levs; ++k) {
        for (size_t i = 0; i < im; ++i) {
            ten_t_pbl[i, k] = 0.0;
            ten_q_pbl[i, k] = 0.0;
        }
    }

    if (nvdiff == ntrac && (hybedmf || do_shoc || satmedmf)) {
        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < ntrac; ++n) {
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    ten_q[i, k, n] = dvdftra[i, k, n];
                }
            }
        }
    } else if (nvdiff != ntrac && !shinhong && !do_ysu) {
        if (ntke > 0) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    ten_q[i, k, f2c(ntke)] = dvdftra[i, k, f2c(ntkev)];
                }
            }
        }

        int kk = 0;
        if (trans_aero) {
            set_aerosol_tracer_index(
                imp_physics, imp_physics_wsm6, imp_physics_thompson, ltaerosol, mraerosol,
                imp_physics_tempo, lthailaware, imp_physics_mg, ntgl, imp_physics_gfdl, imp_physics_nssl,
                nssl_hail_on, nssl_ccn_on, nssl_3moment, kk, errmsg, errflg
            );
            if (errflg != 0) return;

            int k1 = kk;
            for (int n = ntchs; n < ntchm + ntchs; ++n) {
                k1++;
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(n)] = dvdftra[i, k, f2c(k1)];
                    }
                }
            }
        }

        if (imp_physics == imp_physics_wsm6) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    ten_q[i, k, f2c(ntqv)] = dvdftra[i, k, 0];
                    ten_q[i, k, f2c(ntcw)] = dvdftra[i, k, 1];
                    ten_q[i, k, f2c(ntiw)] = dvdftra[i, k, 2];
                    ten_q[i, k, f2c(ntoz)] = dvdftra[i, k, 3];
                }
            }
        } else if (imp_physics == imp_physics_fer_hires) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    ten_q[i, k, f2c(ntqv)]     = dvdftra[i, k, 0];
                    ten_q[i, k, f2c(ntcw)]     = dvdftra[i, k, 1];
                    ten_q[i, k, f2c(ntiw)]     = dvdftra[i, k, 2];
                    ten_q[i, k, f2c(ntrw)]     = dvdftra[i, k, 3];
                    ten_q[i, k, f2c(nqrimef)]  = dvdftra[i, k, 4];
                    ten_q[i, k, f2c(ntoz)]     = dvdftra[i, k, 5];
                }
            }
        } else if (imp_physics == imp_physics_thompson) {
            if (ltaerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 9];
                        ten_q[i, k, f2c(ntwa)]  = dvdftra[i, k, 10];
                        ten_q[i, k, f2c(ntia)]  = dvdftra[i, k, 11];
                    }
                }
            } else if (mraerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 9];
                    }
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 8];
                    }
                }
            }
        } else if (imp_physics == imp_physics_tempo) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    dqdt[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                    dqdt[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                    dqdt[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                    dqdt[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                    dqdt[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                    dqdt[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                    dqdt[i, k, f2c(ntinc)] = dvdftra[i, k, 6];
                    dqdt[i, k, f2c(ntrnc)] = dvdftra[i, k, 7];
                    dqdt[i, k, f2c(ntoz)]  = dvdftra[i, k, 8];
                }
            }
            int n = 9;
            if (ltaerosol) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dqdt[i, k, f2c(ntlnc)] = dvdftra[i, k, n];
                        dqdt[i, k, f2c(ntwa)]  = dvdftra[i, k, n + 1];
                        dqdt[i, k, f2c(ntia)]  = dvdftra[i, k, n + 2];
                    }
                }
                n = 12;
            }
            if (lthailaware) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dqdt[i, k, f2c(ntgnc)] = dvdftra[i, k, n];
                        dqdt[i, k, f2c(ntgv)]  = dvdftra[i, k, n + 1];
                    }
                }
            }
        } else if (imp_physics == imp_physics_mg) {
            if (ntgl > 0) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, 0]          = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntsnc)] = dvdftra[i, k, 9];
                        ten_q[i, k, f2c(ntgnc)] = dvdftra[i, k, 10];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 11];
                    }
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, 0]          = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntsnc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 9];
                    }
                }
            }
        } else if (imp_physics == imp_physics_gfdl) {
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < levs; ++k) {
                for (size_t i = 0; i < im; ++i) {
                    ten_q[i, k, f2c(ntqv)] = dvdftra[i, k, 0];
                    ten_q[i, k, f2c(ntcw)] = dvdftra[i, k, 1];
                    ten_q[i, k, f2c(ntiw)] = dvdftra[i, k, 2];
                    ten_q[i, k, f2c(ntrw)] = dvdftra[i, k, 3];
                    ten_q[i, k, f2c(ntsw)] = dvdftra[i, k, 4];
                    ten_q[i, k, f2c(ntgl)] = dvdftra[i, k, 5];
                    ten_q[i, k, f2c(ntoz)] = dvdftra[i, k, 6];
                }
            }
        } else if (imp_physics == imp_physics_nssl) {
            if (nssl_hail_on) {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(nthl)]  = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 9];
                        ten_q[i, k, f2c(ntsnc)] = dvdftra[i, k, 10];
                        ten_q[i, k, f2c(ntgnc)] = dvdftra[i, k, 11];
                        ten_q[i, k, f2c(nthnc)] = dvdftra[i, k, 12];
                        ten_q[i, k, f2c(ntgv)]  = dvdftra[i, k, 13];
                        ten_q[i, k, f2c(nthv)]  = dvdftra[i, k, 14];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 15];
                    }
                }
                int n = 15;
                if (nssl_ccn_on) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            ten_q[i, k, f2c(ntccn)] = dvdftra[i, k, n + 1];
                        }
                    }
                    n += 1;
                }
                if (nssl_3moment) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            ten_q[i, k, f2c(ntrz)] = dvdftra[i, k, n + 1];
                            ten_q[i, k, f2c(ntgz)] = dvdftra[i, k, n + 2];
                            ten_q[i, k, f2c(nthz)] = dvdftra[i, k, n + 3];
                        }
                    }
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        ten_q[i, k, f2c(ntqv)]  = dvdftra[i, k, 0];
                        ten_q[i, k, f2c(ntcw)]  = dvdftra[i, k, 1];
                        ten_q[i, k, f2c(ntiw)]  = dvdftra[i, k, 2];
                        ten_q[i, k, f2c(ntrw)]  = dvdftra[i, k, 3];
                        ten_q[i, k, f2c(ntsw)]  = dvdftra[i, k, 4];
                        ten_q[i, k, f2c(ntgl)]  = dvdftra[i, k, 5];
                        ten_q[i, k, f2c(ntlnc)] = dvdftra[i, k, 6];
                        ten_q[i, k, f2c(ntinc)] = dvdftra[i, k, 7];
                        ten_q[i, k, f2c(ntrnc)] = dvdftra[i, k, 8];
                        ten_q[i, k, f2c(ntsnc)] = dvdftra[i, k, 9];
                        ten_q[i, k, f2c(ntgnc)] = dvdftra[i, k, 10];
                        ten_q[i, k, f2c(ntgv)]  = dvdftra[i, k, 11];
                        ten_q[i, k, f2c(ntoz)]  = dvdftra[i, k, 12];
                    }
                }
                int n = 12;
                if (nssl_ccn_on) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            ten_q[i, k, f2c(ntccn)] = dvdftra[i, k, n + 1];
                        }
                    }
                    n += 1;
                }
                if (nssl_3moment) {
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            ten_q[i, k, f2c(ntrz)] = dvdftra[i, k, n + 1];
                            ten_q[i, k, f2c(ntgz)] = dvdftra[i, k, n + 2];
                        }
                    }
                }
            }
        }
    }

    if (tend_opt_pbl == 1) { // immediately apply tendencies
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
    } else if (tend_opt_pbl == 2) { // add tendencies to sum
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
    } else if (tend_opt_pbl == 3) { // add tendencies to sum and apply
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
    } else if (tend_opt_pbl == 4) { // Current state unchanged
        // do nothing, exit
    } else {
        errflg = 1;
        errmsg = "A tendency application control was outside of the acceptable range (1-4)";
        return;
    }

    // Surface Coupled Fluxes Calculation
    if (cplflx) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (oceanfrac[i] > 0.0) { // Ocean only, NO LAKES
                if (!wet[i]) { // no open water
                    if (kdt > 1) { // use results from CICE
                        dusfci_cpl[i] = dusfc_cice[i];
                        dvsfci_cpl[i] = dvsfc_cice[i];
                        dtsfci_cpl[i] = dtsfc_cice[i];
                        dqsfci_cpl[i] = dqsfc_cice[i];
                    } else { // use PBL fluxes when CICE is unavailable
                        dusfci_cpl[i] = dusfc1[i];
                        dvsfci_cpl[i] = dvsfc1[i];
                        dtsfci_cpl[i] = dtsfc1[i] * hffac[i];
                        dqsfci_cpl[i] = dqsfc1[i];
                    }
                } else if (icy[i] || dry[i]) { // use stress_ocean from sfc_diff for open water component
                    double rho_val = prsl[i, 0] / (rd * t1[i] * (1.0 + fvirt * std::max(q1[i], qmin)));
                    if (wind[i] > 0.0) {
                        double tem_val = -rho_val * stress_wat[i] / wind[i];
                        dusfci_cpl[i] = tem_val * ugrs1[i];
                        dvsfci_cpl[i] = tem_val * vgrs1[i];
                    } else {
                        dusfci_cpl[i] = 0.0;
                        dvsfci_cpl[i] = 0.0;
                    }
                    dtsfci_cpl[i] = cp * rho_val * hflx_wat[i];
                    dqsfci_cpl[i] = hvap * rho_val * evap_wat[i];
                } else { // 100% open ocean
                    if (use_med_flux && kdt > 1) { // use results from CMEPS mediator
                        dusfci_cpl[i] = dusfc_med[i];
                        dvsfci_cpl[i] = dvsfc_med[i];
                        dtsfci_cpl[i] = dtsfc_med[i];
                        dqsfci_cpl[i] = dqsfc_med[i];
                    } else { // use results from PBL scheme
                        dusfci_cpl[i] = dusfc1[i];
                        dvsfci_cpl[i] = dvsfc1[i];
                        dtsfci_cpl[i] = dtsfc1[i] * hffac[i];
                        dqsfci_cpl[i] = dqsfc1[i];
                    }
                }
                dusfc_cpl[i] += dusfci_cpl[i] * dtf;
                dvsfc_cpl[i] += dvsfci_cpl[i] * dtf;
                dtsfc_cpl[i] += dtsfci_cpl[i] * dtf;
                dqsfc_cpl[i] += dqsfci_cpl[i] * dtf;
            } else {
                dusfc_cpl[i] = huge;
                dvsfc_cpl[i] = huge;
                dtsfc_cpl[i] = huge;
                dqsfc_cpl[i] = huge;
            }
        }
    }

    if (cplchm) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (oceanfrac[i] > 0.0) {
                ushfsfci[i] = dtsfci_cpl[i];
            } else {
                double rho_val = prsl[i, 0] / (rd * t1[i] * (1.0 + fvirt * std::max(q1[i], qmin)));
                ushfsfci[i] = cp * rho_val * hflx[i];
            }
        }
    }

    if (cplaqm) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            if (oceanfrac[i] > 0.0) {
                if (!cplflx) {
                    dtsfci_cpl[i] = dtsfc1[i] * hffac[i];
                    dqsfci_cpl[i] = dqsfc1[i];
                }
            } else { // land
                dtsfci_cpl[i] = dtsfc1[i] * hffac[i];
                dqsfci_cpl[i] = dqsfc1[i];
            }
        }
    }

    if (lssav) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < im; ++i) {
            dusfc_diag[i]  += dusfc1[i] * dtf;
            dvsfc_diag[i]  += dvsfc1[i] * dtf;
            dusfci_diag[i] = dusfc1[i];
            dvsfci_diag[i] = dvsfc1[i];
            dtsfci_diag[i] = dtsfc1[i] * hffac[i];
            dqsfci_diag[i] = dqsfc1[i];
            dtsfc_diag[i]  += dtsfci_diag[i] * dtf;
            dqsfc_diag[i]  += dqsfci_diag[i] * dtf;
        }

        if (ldiag3d && flag_for_pbl_generic_tend && dtend_data != nullptr) {
            View3D dtend(dtend_data, im, levs, static_cast<size_t>(dtidx.extent(1)));

            if (lsidea) {
                int idtend_t = dtidx[f2c(index_of_temperature), f2c(index_of_process_pbl)];
                if (idtend_t >= 1) {
                    size_t id_t = f2c(idtend_t);
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            dtend[i, k, id_t] += ten_t[i, k] * dtf;
                        }
                    }
                }
            } else {
                int idtend_t = dtidx[f2c(index_of_temperature), f2c(index_of_process_pbl)];
                if (idtend_t >= 1) {
                    size_t id_t = f2c(idtend_t);
                    #pragma omp parallel for schedule(static)
                    for (size_t k = 0; k < levs; ++k) {
                        for (size_t i = 0; i < im; ++i) {
                            dtend[i, k, id_t] += ten_t[i, k] * dtp;
                        }
                    }
                }
            }
            int idtend_u = dtidx[f2c(index_of_x_wind), f2c(index_of_process_pbl)];
            if (idtend_u >= 1) {
                size_t id_u = f2c(idtend_u);
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dtend[i, k, id_u] += ten_u[i, k] * dtp;
                    }
                }
            }
            int idtend_v = dtidx[f2c(index_of_y_wind), f2c(index_of_process_pbl)];
            if (idtend_v >= 1) {
                size_t id_v = f2c(idtend_v);
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dtend[i, k, id_v] += ten_v[i, k] * dtp;
                    }
                }
            }
            int idtend_qv = dtidx[f2c(100 + ntqv), f2c(index_of_process_pbl)];
            if (idtend_qv >= 1) {
                size_t id_qv = f2c(idtend_qv);
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dtend[i, k, id_qv] += ten_q[i, k, f2c(ntqv)] * dtp;
                    }
                }
            }
            int idtend_oz = dtidx[f2c(100 + ntoz), f2c(index_of_process_pbl)];
            if (idtend_oz >= 1) {
                size_t id_oz = f2c(idtend_oz);
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dtend[i, k, id_oz] += ten_q[i, k, f2c(ntoz)] * dtp;
                    }
                }
            }
            int idtend_ke = dtidx[f2c(100 + ntke), f2c(index_of_process_pbl)];
            if (idtend_ke >= 1) {
                size_t id_ke = f2c(idtend_ke);
                #pragma omp parallel for schedule(static)
                for (size_t k = 0; k < levs; ++k) {
                    for (size_t i = 0; i < im; ++i) {
                        dtend[i, k, id_ke] += ten_q[i, k, f2c(ntke)] * dtp;
                    }
                }
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < levs; ++k) {
        for (size_t i = 0; i < im; ++i) {
            ten_t_pbl[i, k] = ten_t[i, k];
            ten_q_pbl[i, k] = ten_q[i, k, f2c(ntqv)];
        }
    }
}

} // namespace interstitials
} // namespace satmedmf

#endif // SATMEDMF_PBL_INTERSTITIALS_HPP
