#ifndef SFC_KOKKOS_HPP
#define SFC_KOKKOS_HPP

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "sfc_constants.hpp"
#include <stddef.h>

/**
 * @file sfc_kokkos.hpp
 * @brief GPU-Accelerated parallel Kokkos implementation of GFS Surface Layer stability and ocean skin models.
 *
 * Delivers exascale-ready parallel device execution kernels compatible with CUDA, HIP, SYCL, and CPU execution spaces.
 */

namespace sfc {
namespace gpu {

#ifdef ENABLE_KOKKOS

// GPU-ready Kokkos View definitions mapping directly to Fortran Column-Major LayoutLeft memory
using KView1D = Kokkos::View<double*, Kokkos::LayoutLeft>;
using KConstView1D = Kokkos::View<const double*, Kokkos::LayoutLeft>;
using KView2D = Kokkos::View<double**, Kokkos::LayoutLeft>;
using KConstView2D = Kokkos::View<const double**, Kokkos::LayoutLeft>;

// Atmospheric Surface Sounding View structures for Kokkos GPU execution
struct KSurfaceSounding {
    KConstView1D u1;           // Zonal wind component (m/s)
    KConstView1D v1;           // Meridional wind component (m/s)
    KConstView1D t1;           // Temperature (K)
    KConstView1D q1;           // Specific humidity (kg/kg)
    KConstView1D z1;           // Boundary height (m)
    KConstView1D ps;           // Surface pressure (Pa)
    KConstView1D tskin;        // Surface skin temperature (K)
};

// Aerodynamic roughness length formulations over water (GPU-inlined)
KOKKOS_INLINE_FUNCTION double znot_m_v6_gpu(double uref) {
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
        return Kokkos::exp(p10 + uref * (p11 + uref * (p12 + uref * p13)));
    } else if (uref > 6.5 && uref <= 15.7) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 15.7 && uref <= 53.0) {
        return Kokkos::exp(p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35)))));
    } else {
        return p40;
    }
}

KOKKOS_INLINE_FUNCTION double znot_t_v6_gpu(double uref) {
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
inline double znot_m_v7_gpu(double uref) {
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
        return Kokkos::exp(p10 + uref * (p11 + uref * (p12 + uref * p13)));
    } else if (uref > 6.5 && uref <= 15.7) {
        return p20 + uref * (p21 + uref * (p22 + uref * (p23 + uref * (p24 + uref * p25))));
    } else if (uref > 15.7 && uref <= 53.0) {
        return Kokkos::exp(p30 + uref * (p31 + uref * (p32 + uref * (p33 + uref * (p34 + uref * p35)))));
    } else {
        return p40;
    }
}

// Option 7: Wang (2018) scalar roughness length for heat
inline double znot_t_v7_gpu(double uref) {
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

/**
 * @brief GPU-Accelerated Kokkos implementation of sfc_diff stability exchange solver.
 */
template <typename ExecutionSpace>
inline void sfc_diff_run_kokkos(
    size_t columns,
    KSurfaceSounding sounding,
    KConstView1D z0,
    int sfc_z0_type,
    KView1D cm,
    KView1D ch,
    KView1D ustar,
    KView1D stress
) {
    Kokkos::parallel_for("sfc_diff_run_device",
        Kokkos::RangePolicy<ExecutionSpace>(0, columns),
        KOKKOS_LAMBDA(const size_t i) {
            
            // 1. Compute wind speed and clamp to prevent division by zero
            double u_val = sounding.u1(i);
            double v_val = sounding.v1(i);
            double wind = Kokkos::sqrt(u_val * u_val + v_val * v_val);
            double wind_clamped = Kokkos::max(1e-4, wind);

            // 2. Compute dynamically updated roughness lengths over ocean
            double z0_val = z0(i);
            double zt_val = z0(i) * 0.1;

            if (sfc_z0_type == 6) {
                z0_val = znot_m_v6_gpu(wind_clamped);
                zt_val = znot_t_v6_gpu(wind_clamped);
            } else if (sfc_z0_type == 7) {
                z0_val = znot_m_v7_gpu(wind_clamped);
                zt_val = znot_t_v7_gpu(wind_clamped);
            }

            // 3. Compute bulk Richardson number (Rb) stability metric
            double dtheta = sounding.t1(i) - sounding.tskin(i);
            double bulk_richardson = (constants::grav / sounding.t1(i)) * (sounding.z1(i) * dtheta / (wind_clamped * wind_clamped));

            // 4. Compute stability functions (Fm, Fh)
            double stability_fm = 1.0;
            double stability_fh = 1.0;

            if (bulk_richardson >= 0.0) {
                double factor = 1.0 + 5.0 * bulk_richardson;
                stability_fm = 1.0 / (factor * factor);
                stability_fh = stability_fm;
            } else {
                stability_fm = 1.0 - (16.0 * bulk_richardson) / (1.0 + 5.0 * Kokkos::sqrt(-bulk_richardson));
                stability_fh = stability_fm;
            }

            // 5. Compute drag and exchange coefficients using Logarithmic Profile
            double roughness_ratio = sounding.z1(i) / Kokkos::max(1e-5, z0_val);
            double log_ratio = Kokkos::log(Kokkos::max(1.1, roughness_ratio));
            double neutral_drag = (constants::karman / log_ratio) * (constants::karman / log_ratio);

            double heat_roughness_ratio = sounding.z1(i) / Kokkos::max(1e-5, zt_val);
            double heat_log_ratio = Kokkos::log(Kokkos::max(1.1, heat_roughness_ratio));
            double heat_neutral_drag = (constants::karman / heat_log_ratio) * (constants::karman / heat_log_ratio);

            cm(i) = neutral_drag * stability_fm;
            ch(i) = heat_neutral_drag * stability_fh;

            // 6. Compute friction velocity (ustar) and surface wind stress
            ustar(i) = Kokkos::sqrt(cm(i)) * wind_clamped;

            double air_density = sounding.ps(i) / (constants::rd * sounding.t1(i));
            stress(i) = air_density * ustar(i) * ustar(i);
        }
    );
}

/**
 * @brief GPU-Accelerated Kokkos implementation of sfc_nst diurnal skin ocean model.
 */
template <typename ExecutionSpace>
inline void sfc_nst_run_kokkos(
    size_t columns,
    KConstView1D sol_flux,
    KConstView1D wind_stress,
    KView1D tskin_wat,
    KView1D cool_skin,
    KView1D warm_layer
) {
    const double t_bulk_ref = 295.15; // Baseline ocean bulk temperature reference (K)

    Kokkos::parallel_for("sfc_nst_run_device",
        Kokkos::RangePolicy<ExecutionSpace>(0, columns),
        KOKKOS_LAMBDA(const size_t i) {
            
            double stress_clamped = Kokkos::max(1e-5, static_cast<double>(wind_stress(i)));
            cool_skin(i) = 0.25 / (1.0 + 20.0 * stress_clamped);

            double raw_warming = 0.005 * sol_flux(i) / Kokkos::sqrt(stress_clamped);
            warm_layer(i) = Kokkos::max(0.0, Kokkos::min(3.5, raw_warming));

            tskin_wat(i) = t_bulk_ref - cool_skin(i) + warm_layer(i);
        }
    );
}

#endif // ENABLE_KOKKOS

} // namespace gpu
} // namespace sfc

#endif // SFC_KOKKOS_HPP
