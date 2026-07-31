#include "satmedmf_interface.hpp"
#include "satmedmf_vdifq.hpp"

#ifdef ENABLE_KOKKOS
#include "satmedmf_kokkos.hpp"
#endif

extern "C" {

void c_satmedmf_run(
    size_t layers,
    size_t columns,
    double dt,
    int tte_edmf,
    double* t_lay,
    const double* p_lay,
    const double* rho,
    double* u_wind,
    double* v_wind,
    double* q_vap,
    double* te,
    const double* heat,
    const double* evap,
    const double* stress,
    int do_canopy,
    const double* cfch,
    double* dt_temp,
    double* dt_u,
    double* dt_v,
    double* dt_te
) {
#ifdef ENABLE_KOKKOS
    if (!Kokkos::is_initialized()) {
        Kokkos::initialize();
    }
    {
        // Wrap raw pointers directly into unmanaged Kokkos Host Views
        using KHostView2D = Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using KConstHostView2D = Kokkos::View<const double**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using KHostView1D = Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using KConstHostView1D = Kokkos::View<const double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        KHostView2D t_lay_h(t_lay, columns, layers);
        KHostView2D u_wind_h(u_wind, columns, layers);
        KHostView2D v_wind_h(v_wind, columns, layers);
        KHostView2D q_vap_h(q_vap, columns, layers);
        KHostView2D te_h(te, columns, layers);

        KConstHostView2D p_lay_h(p_lay, columns, layers);
        KConstHostView2D rho_h(rho, columns, layers);

        KHostView2D dt_temp_h(dt_temp, columns, layers);
        KHostView2D dt_u_h(dt_u, columns, layers);
        KHostView2D dt_v_h(dt_v, columns, layers);
        KHostView2D dt_te_h(dt_te, columns, layers);

        KConstHostView1D heat_h(heat, columns);
        KConstHostView1D evap_h(evap, columns);
        KConstHostView1D stress_h(stress, columns);
        KConstHostView1D cfch_h(cfch, columns);

        // If execution is on a GPU, we must copy to Device Views and copy back.
        // If execution is on CPU (e.g., OpenMP, Serial), we can run directly on Host Views.
        if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace::execution_space>) {
            // CPU execution space, run directly on Host Views with zero copy
            satmedmf::gpu::KPhysicalState state{ t_lay_h, u_wind_h, v_wind_h, q_vap_h, te_h };
            satmedmf::gpu::KDiagnosticTendencies tendencies{ dt_temp_h, dt_u_h, dt_v_h, dt_te_h };

            satmedmf::gpu::satmedmf_run_kokkos<Kokkos::DefaultExecutionSpace>(
                layers, columns, dt, tte_edmf != 0,
                state, p_lay_h, rho_h,
                heat_h, evap_h, stress_h,
                do_canopy != 0, cfch_h,
                tendencies
            );
        } else {
            // GPU execution space, allocate managed Device Views and copy
            using KDeviceView2D = Kokkos::View<double**, Kokkos::LayoutLeft, typename Kokkos::DefaultExecutionSpace::memory_space>;
            using KConstDeviceView2D = Kokkos::View<const double**, Kokkos::LayoutLeft, typename Kokkos::DefaultExecutionSpace::memory_space>;
            using KDeviceView1D = Kokkos::View<double*, Kokkos::LayoutLeft, typename Kokkos::DefaultExecutionSpace::memory_space>;
            using KConstDeviceView1D = Kokkos::View<const double*, Kokkos::LayoutLeft, typename Kokkos::DefaultExecutionSpace::memory_space>;

            KDeviceView2D t_lay_d("t_lay_d", columns, layers);
            KDeviceView2D u_wind_d("u_wind_d", columns, layers);
            KDeviceView2D v_wind_d("v_wind_d", columns, layers);
            KDeviceView2D q_vap_d("q_vap_d", columns, layers);
            KDeviceView2D te_d("te_d", columns, layers);

            KDeviceView2D p_lay_d("p_lay_d", columns, layers);
            KDeviceView2D rho_d("rho_d", columns, layers);

            KDeviceView2D dt_temp_d("dt_temp_d", columns, layers);
            KDeviceView2D dt_u_d("dt_u_d", columns, layers);
            KDeviceView2D dt_v_d("dt_v_d", columns, layers);
            KDeviceView2D dt_te_d("dt_te_d", columns, layers);

            KDeviceView1D heat_d("heat_d", columns);
            KDeviceView1D evap_d("evap_d", columns);
            KDeviceView1D stress_d("stress_d", columns);
            KConstDeviceView1D cfch_d("cfch_d", columns);

            // Copy to Device
            Kokkos::deep_copy(t_lay_d, t_lay_h);
            Kokkos::deep_copy(u_wind_d, u_wind_h);
            Kokkos::deep_copy(v_wind_d, v_wind_h);
            Kokkos::deep_copy(q_vap_d, q_vap_h);
            Kokkos::deep_copy(te_d, te_h);
            Kokkos::deep_copy(p_lay_d, p_lay_h);
            Kokkos::deep_copy(rho_d, rho_h);
            Kokkos::deep_copy(dt_temp_d, dt_temp_h);
            Kokkos::deep_copy(dt_u_d, dt_u_h);
            Kokkos::deep_copy(dt_v_d, dt_v_h);
            Kokkos::deep_copy(dt_te_d, dt_te_h);
            Kokkos::deep_copy(heat_d, heat_h);
            Kokkos::deep_copy(evap_d, evap_h);
            Kokkos::deep_copy(stress_d, stress_h);
            Kokkos::deep_copy(cfch_d, cfch_h);

            satmedmf::gpu::KPhysicalState state{ t_lay_d, u_wind_d, v_wind_d, q_vap_d, te_d };
            satmedmf::gpu::KDiagnosticTendencies tendencies{ dt_temp_d, dt_u_d, dt_v_d, dt_te_d };

            satmedmf::gpu::satmedmf_run_kokkos<Kokkos::DefaultExecutionSpace>(
                layers, columns, dt, tte_edmf != 0,
                state, p_lay_d, rho_d,
                heat_d, evap_d, stress_d,
                do_canopy != 0, cfch_d,
                tendencies
            );

            // Copy back results
            Kokkos::deep_copy(t_lay_h, t_lay_d);
            Kokkos::deep_copy(u_wind_h, u_wind_d);
            Kokkos::deep_copy(v_wind_h, v_wind_d);
            Kokkos::deep_copy(q_vap_h, q_vap_d);
            Kokkos::deep_copy(te_h, te_d);
            Kokkos::deep_copy(dt_temp_h, dt_temp_d);
            Kokkos::deep_copy(dt_u_h, dt_u_d);
            Kokkos::deep_copy(dt_v_h, dt_v_d);
            Kokkos::deep_copy(dt_te_h, dt_te_d);
        }
    }
#else
    // Construct LayoutLeft column-major mdspan views
    satmedmf::View2D t_lay_view(t_lay, columns, layers);
    satmedmf::View2D u_wind_view(u_wind, columns, layers);
    satmedmf::View2D v_wind_view(v_wind, columns, layers);
    satmedmf::View2D q_vap_view(q_vap, columns, layers);
    satmedmf::View2D te_view(te, columns, layers);

    satmedmf::ConstView2D p_lay_view(p_lay, columns, layers);
    satmedmf::ConstView2D rho_view(rho, columns, layers);

    satmedmf::View2D dt_temp_view(dt_temp, columns, layers);
    satmedmf::View2D dt_u_view(dt_u, columns, layers);
    satmedmf::View2D dt_v_view(dt_v, columns, layers);
    satmedmf::View2D dt_te_view(dt_te, columns, layers);

    satmedmf::PhysicalState state{
        t_lay_view,
        u_wind_view,
        v_wind_view,
        q_vap_view,
        te_view
    };

    satmedmf::DiagnosticTendencies tendencies{
        dt_temp_view,
        dt_u_view,
        dt_v_view,
        dt_te_view
    };

    satmedmf::satmedmf_run_core(
        layers,
        columns,
        dt,
        tte_edmf != 0,
        state,
        p_lay_view,
        rho_view,
        heat,
        evap,
        stress,
        do_canopy != 0,
        cfch,
        tendencies
    );
#endif
}

}
