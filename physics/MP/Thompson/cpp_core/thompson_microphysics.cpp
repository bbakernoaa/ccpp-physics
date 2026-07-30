#include "thompson_microphysics_interface.hpp"
#include "thompson_microphysics.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void c_thompson_microphysics_run(
    size_t layers,
    size_t columns,
    double dt,
    double* t_lay,
    const double* p_lay,
    const double* rho,
    double* qv,
    double* qc,
    double* qr,
    double* qi,
    double* qs,
    double* qg,
    double* ni,
    double* nr,
    double* ns,
    double* ng,
    double* precip
) {
    // Wrap standard contiguous flat pointers in C++23 std::mdspan LayoutLeft views (columns, layers)
    thompson::View2D t_lay_view(t_lay, columns, layers);
    thompson::ConstView2D p_lay_view(p_lay, columns, layers);
    thompson::ConstView2D rho_view(rho, columns, layers);
    thompson::View2D qv_view(qv, columns, layers);
    thompson::View2D qc_view(qc, columns, layers);
    thompson::View2D qr_view(qr, columns, layers);
    thompson::View2D qi_view(qi, columns, layers);
    thompson::View2D qs_view(qs, columns, layers);
    thompson::View2D qg_view(qg, columns, layers);
    thompson::View2D ni_view(ni, columns, layers);
    thompson::View2D nr_view(nr, columns, layers);
    thompson::View2D ns_view(ns, columns, layers);
    thompson::View2D ng_view(ng, columns, layers);

    // Execute Thompson Microphysics C++23 calculations
    thompson::thompson_microphysics_run_core(
        layers, columns, dt,
        t_lay_view, p_lay_view, rho_view,
        qv_view, qc_view, qr_view, qi_view, qs_view, qg_view,
        ni_view, nr_view, ns_view, ng_view,
        precip
    );
}

#ifdef __cplusplus
}
#endif
