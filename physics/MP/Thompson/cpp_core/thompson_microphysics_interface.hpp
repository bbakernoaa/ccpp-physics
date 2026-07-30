#ifndef THOMPSON_MICROPHYSICS_INTERFACE_HPP
#define THOMPSON_MICROPHYSICS_INTERFACE_HPP

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Public flat C ABI entry point for the C++23 Thompson Microphysics solver.
 * Wraps contiguous pointers natively in mdspan views inside the implementation.
 *
 * @param layers    Number of vertical grid layers (passed by value)
 * @param columns   Number of horizontal columns (passed by value)
 * @param t_lay     Layer temperatures array pointer (Kelvin)
 * @param p_lay     Layer pressures array pointer (Pascals)
 * @param rho       Layer air density array pointer (kg/m3)
 * @param qv        Water vapor mixing ratio pointer (kg/kg)
 * @param qc        Cloud water mixing ratio pointer (kg/kg)
 * @param qr        Rain water mixing ratio pointer (kg/kg)
 * @param qi        Cloud ice mixing ratio pointer (kg/kg)
 * @param qs        Snow mixing ratio pointer (kg/kg)
 * @param qg        Graupel mixing ratio pointer (kg/kg)
 * @param ni        Cloud ice number concentration pointer (1/kg)
 * @param nr        Rain number concentration pointer (1/kg)
 * @param ns        Snow number concentration pointer (1/kg)
 * @param ng        Graupel number concentration pointer (1/kg)
 * @param precip    Surface precipitation rates pointer [rain, snow, graupel] (m/s)
 */
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
);

#ifdef __cplusplus
}
#endif

#endif // THOMPSON_MICROPHYSICS_INTERFACE_HPP
