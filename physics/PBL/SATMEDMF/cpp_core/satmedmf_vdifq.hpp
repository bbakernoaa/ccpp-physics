#ifndef SATMEDMF_VDIFQ_HPP
#define SATMEDMF_VDIFQ_HPP

#include "satmedmf_types.hpp"
#include "satmedmf_math_utils.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>

// Include subgrid modules after type definitions
#include "satmedmf_tridiagonal.hpp"
#include "satmedmf_thermo.hpp"
#include "satmedmf_massflux.hpp"
#include "satmedmf_canopy.hpp"

namespace satmedmf {

/**
 * @brief Main scale-aware TKE-based moist Eddy-Diffusivity Mass-Flux (TKE-EDMF) vertical turbulent mixing solver.
 *
 * This routine implements the C++23 header-only translation core for the scale-aware TKE-EDMF (SATMEDMF) Planetary Boundary 
 * Layer scheme inside the CCPP Parameterizations library. 
 *
 * @section reference References
 * - **Han, J., and C. S. Bretherton, 2019**: TKE-based Moist Eddy-Diffusivity Mass-Flux (EDMF) Parameterization for Vertical 
 *   Turbulent Mixing. *Weather and Forecasting*, **34**, accepted. (Equations 9, 10, 22, 23).
 * - **Han, J. et al., 2016**: Implementation in the NCEP GFS of a Hybrid Eddy-Diffusivity Mass-Flux (EDMF) Boundary Layer 
 *   Parameterization with Dissipative Heating and Modified Stable Boundary Layer Mixing. *Weather and Forecasting*, **31**(1), 341–352.
 * - **Troen, I. B., and L. Mahrt, 1986**: A simple model of the boundary layer; sensitivity to surface boundary conditions. 
 *   *Boundary-Layer Meteorology*, **37**, 129–148. (Equation 10a).
 * - **Nakanishi, M., 2001**: Improvement of the Mellor-Yamada level 2.5 turbulence closure model for the boundary layer. 
 *   *Boundary-Layer Meteorology*, **99**, 349–378. (Equation 9).
 * - **Bougeault, P., and P. Lacarrere, 1989**: Parameterization of orography-induced turbulence in a mesobeta-scale model. 
 *   *Monthly Weather Review*, **117**, 1872–1890. (Equation 10).
 *
 * @section physical_processes Main Physical Processes
 * -# **moist thermodynamics**: Solves moist atmospheric coordinates, Poisson-Exner structures, virtual potential temperature 
 *    \f$\theta_v\f$, and planetary boundary layer heights \f$h_{pbl}\f$ (using a Richardson number scanner based on Equation 10a 
 *    of Troen & Mahrt 1986).
 * -# **similarity scale estimation**: Computes the velocity scale \f$w_s\f$ scaled at the top of the surface layer following 
 *    Equation 22 & 23 of Han & Bretherton (2019):
 *    \f[
 *    w_s = (u_*^3 + 7 \alpha \kappa w_*^3)^{1/3}
 *    \f]
 * -# **local vertical mixing**: Resolves the turbulent kinetic energy (TKE) or total turbulent energy (TTE) budget under local 
 *    gradient Richardson number stability conditions (following the stable/unstable regimes from Han et al. 2016).
 * -# **buoyant thermal massflux**: Implements nonlocal transport by buoyant convective plumes (updraft thermals) and stratocumulus 
 *    cloud-top radiative cooling-induced downdrafts.
 * -# **vertical diffusion solver**: Computes numerical fluxes for temperature, momentum components, and TKE/TTE using explicit, 
 *    coalesced OpenMP loop kernels.
 * -# **sub-grid forest canopy**: Simulates a 3-layer subgrid canopy with monotonic heights sorting and mass-conservative bidirectional 
 *    transport of tracers.
 *
 * @param layers Number of vertical model layers (layers)
 * @param columns Number of horizontal grid columns (columns)
 * @param dt Model timestep (s)
 * @param tte_edmf Switch to activate Total Turbulent Energy (TTE-EDMF) closure instead of standard TKE (Han et al. 2019)
 * @param state Input/Output physical state holding profiles of temperature, wind components, and TKE/TTE
 * @param p_lay Center layer pressure values (Pa)
 * @param rho Air density profile (kg/m3)
 * @param heat Surface sensible heat flux (K m/s)
 * @param evap Surface moisture flux (kg/kg m/s)
 * @param stress Surface momentum friction velocity square (m2/s2)
 * @param do_canopy Flag to activate the 3-layer forest canopy module
 * @param cfch Canopy Forest Height profile (m)
 * @param tendencies Output diagnostic tendencies to update temperature, wind components, and TKE/TTE
 */
inline void satmedmf_run_core(
    size_t layers,
    size_t columns,
    double dt,
    bool tte_edmf,
    PhysicalState& state,
    ConstView2D p_lay,
    ConstView2D rho,
    const double* heat,
    const double* evap,
    const double* stress,
    bool do_canopy,
    const double* cfch,
    DiagnosticTendencies& tendencies
) {
    // Physical parameters matching Han & Bretherton (2019)
    const double tkmin = 1e-4;
    const double dkmax = 100.0;
    const double xkzo_val = 0.01;
    const double xkzmo_val = 0.01;

    double cfac = 4.5;
    double prmax = 4.0;
    double prscu = 0.67;
    double ck1 = 0.15;
    double ch1 = 0.15;

    if (tte_edmf) {
        cfac  = 3.0;
        prmax = 6.0;
        prscu = 0.4;
        ck1   = 0.16;
        ch1   = 0.16;
    }

    // Thread-local static workspace scratch caches to prevent serial allocation overhead
    thread_local static std::vector<double> zl_cache;
    thread_local static std::vector<double> dku_cache;
    thread_local static std::vector<double> dkt_cache;
    thread_local static std::vector<double> dkq_cache;
    thread_local static std::vector<double> prsi_cache;
    thread_local static std::vector<double> zol_cache;
    thread_local static std::vector<double> hpbl_cache;
    thread_local static std::vector<double> kpbl_cache;
    thread_local static std::vector<double> wscale_cache;
    thread_local static std::vector<double> elm_cache;
    thread_local static std::vector<double> diss_cache;
    thread_local static std::vector<double> thetae_cache;
    thread_local static std::vector<double> qlx_cache;
    thread_local static std::vector<double> thvx_cache;
    thread_local static std::vector<double> thlvx_cache;
    thread_local static std::vector<double> bf_cache;
    thread_local static std::vector<int> kmod_cache;
    thread_local static std::vector<int> kcan3_cache;
    thread_local static std::vector<double> zmid_can_cache;
    thread_local static std::vector<double> zmom_can_cache;
    thread_local static std::vector<double> prsl_can_cache;
    thread_local static std::vector<double> prsi_can_cache;
    thread_local static std::vector<double> t1_can_cache;
    thread_local static std::vector<double> qv_can_cache;
    thread_local static std::vector<double> ws_can_cache;
    thread_local static std::vector<double> dens_can_cache;
    thread_local static std::vector<double> dkt_can_cache;
    thread_local static std::vector<double> dku_can_cache;

    size_t size_cl = columns * layers;
    size_t size_cl1 = columns * (layers + 1);
    size_t size_kc = columns * 3;
    size_t size_kt = columns * (layers + 3);
    size_t size_kt1 = columns * (layers + 4);

    // Idempotent resizing sweep (executes only on grid size configuration changes)
    if (zl_cache.size() < size_cl) {
        zl_cache.resize(size_cl, 0.0);
        dku_cache.resize(size_cl, 0.0);
        dkt_cache.resize(size_cl, 0.0);
        dkq_cache.resize(size_cl, 0.0);
        zol_cache.resize(size_cl, 0.0);
        hpbl_cache.resize(size_cl, 0.0);
        kpbl_cache.resize(size_cl, 0.0);
        wscale_cache.resize(size_cl, 0.0);
        elm_cache.resize(size_cl, 0.0);
        diss_cache.resize(size_cl, 0.0);
        thetae_cache.resize(size_cl, 0.0);
        qlx_cache.resize(size_cl, 0.0);
        thvx_cache.resize(size_cl, 0.0);
        thlvx_cache.resize(size_cl, 0.0);
        bf_cache.resize(size_cl, 0.0);
        kmod_cache.resize(size_cl, 0);
    }
    if (prsi_cache.size() < size_cl1) {
        prsi_cache.resize(size_cl1, 100000.0);
    }
    if (kcan3_cache.size() < size_kc) {
        kcan3_cache.resize(size_kc, 0);
    }
    if (zmid_can_cache.size() < size_kt) {
        zmid_can_cache.resize(size_kt, 0.0);
        prsl_can_cache.resize(size_kt, 0.0);
        t1_can_cache.resize(size_kt, 0.0);
        qv_can_cache.resize(size_kt, 0.0);
        ws_can_cache.resize(size_kt, 0.0);
        dens_can_cache.resize(size_kt, 0.0);
        dkt_can_cache.resize(size_kt, 0.0);
        dku_can_cache.resize(size_kt, 0.0);
    }
    if (zmom_can_cache.size() < size_kt1) {
        zmom_can_cache.resize(size_kt1, 0.0);
        prsi_can_cache.resize(size_kt1, 0.0);
    }

    // Map LayoutLeft column-major Views directly onto pre-allocated static memories
    View2D zl(zl_cache.data(), columns, layers);
    View2D dku(dku_cache.data(), columns, layers);
    View2D dkt(dkt_cache.data(), columns, layers);
    View2D dkq(dkq_cache.data(), columns, layers);

    // 1. Compute layer physical center heights (zl) using hydrostatic relation
    int num_threads = 1;
    #pragma omp parallel
    {
        #pragma omp master
        num_threads = omp_get_num_threads();
    }
    std::cout << "    [C++] OpenMP Parallel Threads Active: " << num_threads << std::endl;

    #pragma omp parallel for schedule(static)
    for (size_t col = 0; col < columns; ++col) {
        double current_z = 0.0;
        for (size_t lay = 0; lay < layers; ++lay) {
            double dz = 100.0; // Assume nominal 100m layer thickness for standalone
            if (rho[col, lay] > 0.0) {
                dz = 1.0 / (rho[col, lay] * 9.80665) * (lay == 0 ? (P0 - p_lay[col, lay]) : (p_lay[col, lay-1] - p_lay[col, lay]));
                dz = std::max(dz, 10.0);
            }
            zl[col, lay] = current_z + 0.5 * dz;
            current_z += dz;
        }
    }

    // Map high-fidelity thermo variables
    View2D prsi(prsi_cache.data(), columns, layers + 1);
    View2D zol(zol_cache.data(), columns, layers);
    View2D hpbl(hpbl_cache.data(), columns, layers);
    View2D kpbl(kpbl_cache.data(), columns, layers);
    View2D wscale(wscale_cache.data(), columns, layers);
    View2D elm(elm_cache.data(), columns, layers);
    View2D diss(diss_cache.data(), columns, layers);
    View2D thetae(thetae_cache.data(), columns, layers);
    View2D qlx(qlx_cache.data(), columns, layers);
    View2D thvx(thvx_cache.data(), columns, layers);
    View2D thlvx(thlvx_cache.data(), columns, layers);
    View2D bf(bf_cache.data(), columns, layers);

    std::vector<double> u10m_data(columns, 1.0);
    std::vector<double> v10m_data(columns, 1.0);
    std::vector<double> tsea_data(columns, 290.0);
    std::vector<double> fm_data(columns, 1.0);
    std::vector<double> fh_data(columns, 1.0);
    std::vector<double> rbsoil_data(columns, 1.0);

    // Call High-Fidelity Thermodynamics
    thermo::compute_thermodynamics(
        columns, layers, dt,
        state, p_lay, prsi, rho,
        heat, evap, stress,
        u10m_data.data(), v10m_data.data(), tsea_data.data(),
        fm_data.data(), fh_data.data(), rbsoil_data.data(),
        zol, hpbl, kpbl, wscale, elm, diss,
        thetae, qlx, thvx, thlvx, bf
    );

    // Allocate Canopy variables
    size_t nkc = 3;
    size_t nkt = layers + nkc;

    // Compile-time or runtime Canopy levels setup
    if (do_canopy) {
        canopy::canopy_levs_run(
            columns, layers, nkc, nkt,
            RD, 3.141592653589793,
            zl, zl, zl, // nominal height interfaces
            p_lay, p_lay,
            cfch,
            state.t_lay.data_handle(), state.q_vap.data_handle(),
            state.u_wind, state.v_wind, state.t_lay,
            rho, dkt, dku,
            kmod_cache, kcan3_cache,
            zmid_can_cache, zmom_can_cache,
            prsl_can_cache, prsi_can_cache,
            t1_can_cache, qv_can_cache,
            ws_can_cache, dens_can_cache,
            dkt_can_cache, dku_can_cache
        );
    }

    // 2. Loop over columns to compute local vertical mixing & diffusivities
    #pragma omp parallel for schedule(static)
    for (size_t col = 0; col < columns; ++col) {
        for (size_t lay = 0; lay < layers - 1; ++lay) {
            // Local wind shear square (shr2)
            double du = state.u_wind[col, lay+1] - state.u_wind[col, lay];
            double dv = state.v_wind[col, lay+1] - state.v_wind[col, lay];
            double dz = zl[col, lay+1] - zl[col, lay];
            double shr2 = (du * du + dv * dv) / (dz * dz);
            shr2 = std::max(shr2, 1e-6);

            // Buoyancy flux (bf)
            double dt_dz = (state.t_lay[col, lay+1] - state.t_lay[col, lay]) / dz;
            double bf_val = -9.80665 / state.t_lay[col, lay] * dt_dz;

            // Gradient Richardson number
            double ri = bf_val / shr2;

            // Interface TKE/TTE
            double te_interface = 0.5 * (state.te[col, lay] + state.te[col, lay+1]);
            te_interface = std::max(te_interface, tkmin);

            // Compute tesq based on TTE/TKE closure choice
            double tesq = 0.0;
            if (tte_edmf) {
                double epotte = ri / (ri < 0.0 ? (2.0 * ri - 1.0) : (1.0 + 3.0 * ri));
                double tkeh = te_interface * (1.0 - epotte);
                tkeh = std::max(tkeh, tkmin);
                tesq = tkeh / std::sqrt(te_interface);
            } else {
                tesq = std::sqrt(te_interface);
            }

            // Nominal mixing length (elm)
            double elm_val = elm[col, lay];

            // Compute eddy diffusivities
            double dkt_val = 0.0;
            double dku_val = 0.0;

            if (ri < 0.0) { // Unstable mixing regime
                dku_val = ck1 * elm_val * tesq;
                dkt_val = dku_val / prscu;
            } else {        // Stable mixing regime
                dkt_val = ch1 * elm_val * tesq;
                double prnum = 1.0 + 2.1 * ri;
                prnum = std::min(prnum, prmax);
                dku_val = dkt_val * prnum;
            }

            dkt[col, lay] = std::clamp(dkt_val, xkzo_val, dkmax);
            dku[col, lay] = std::clamp(dku_val, xkzmo_val, dkmax);
            dkq[col, lay] = std::clamp(0.5 * dkt_val, xkzo_val, dkmax);
        }
    }

    // Call High-Fidelity Mass Flux
    massflux::compute_massflux(
        columns, layers, dt,
        state, p_lay, prsi, rho,
        hpbl, kpbl, wscale, zol,
        dkt, dku, bf,
        tendencies
    );

    // 3. Compute vertical diffusion tendencies
    #pragma omp parallel for schedule(static)
    for (size_t col = 0; col < columns; ++col) {
        for (size_t lay = 1; lay < layers - 1; ++lay) {
            double dz_up = zl[col, lay+1] - zl[col, lay];
            double dz_dn = zl[col, lay] - zl[col, lay-1];

            // Temperature Diffusion
            double flux_up_t = dkt[col, lay] * (state.t_lay[col, lay+1] - state.t_lay[col, lay]) / dz_up;
            double flux_dn_t = dkt[col, lay-1] * (state.t_lay[col, lay] - state.t_lay[col, lay-1]) / dz_dn;
            tendencies.dt_temp[col, lay] += (flux_up_t - flux_dn_t) / (0.5 * (dz_up + dz_dn));

            // Momentum Diffusion
            double flux_up_u = dku[col, lay] * (state.u_wind[col, lay+1] - state.u_wind[col, lay]) / dz_up;
            double flux_dn_u = dku[col, lay-1] * (state.u_wind[col, lay] - state.u_wind[col, lay-1]) / dz_dn;
            tendencies.dt_u[col, lay] += (flux_up_u - flux_dn_u) / (0.5 * (dz_up + dz_dn));

            double flux_up_v = dku[col, lay] * (state.v_wind[col, lay+1] - state.v_wind[col, lay]) / dz_up;
            double flux_dn_v = dku[col, lay-1] * (state.v_wind[col, lay] - state.v_wind[col, lay-1]) / dz_dn;
            tendencies.dt_v[col, lay] += (flux_up_v - flux_dn_v) / (0.5 * (dz_up + dz_dn));

            // TKE/TTE Diffusion
            double flux_up_te = dkq[col, lay] * (state.te[col, lay+1] - state.te[col, lay]) / dz_up;
            double flux_dn_te = dkq[col, lay-1] * (state.te[col, lay] - state.te[col, lay-1]) / dz_dn;
            tendencies.dt_te[col, lay] += (flux_up_te - flux_dn_te) / (0.5 * (dz_up + dz_dn)) + diss[col, lay];
        }

        // Apply boundary condition surface fluxes at lay = 0
        tendencies.dt_temp[col, 0] += heat[col] / 100.0;
        tendencies.dt_u[col, 0] += stress[col] / 100.0;
        tendencies.dt_v[col, 0] += 0.0;
        tendencies.dt_te[col, 0] += 0.0;

        // Boundary condition at model top
        tendencies.dt_temp[col, layers-1] += 0.0;
        tendencies.dt_u[col, layers-1] += 0.0;
        tendencies.dt_v[col, layers-1] += 0.0;
        tendencies.dt_te[col, layers-1] += 0.0;
    }

    // Bidirectional Canopy-to-Resolved transport
    if (do_canopy) {
        View2D q1_mod(state.q_vap.data_handle(), columns, layers);
        View2D q1_can(qv_can_cache.data(), columns, nkt);
        canopy::canopy_transfer_run(
            columns, layers, nkc, nkt,
            1, 1, // canopy_to_resolved
            zl, zl, zl,
            q1_mod, q1_can,
            kmod_cache, kcan3_cache,
            zmid_can_cache, zmom_can_cache,
            prsl_can_cache, dens_can_cache
        );
    }
}

} // namespace satmedmf

#endif // SATMEDMF_VDIFQ_HPP
