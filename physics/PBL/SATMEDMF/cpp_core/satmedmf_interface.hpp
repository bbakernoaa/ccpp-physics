#ifndef SATMEDMF_INTERFACE_HPP
#define SATMEDMF_INTERFACE_HPP

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Public C ABI entrypoint to invoke the C++23 SATMEDMF planetary boundary layer solver.
 * Wraps flat contiguous pointers into std::mdspan layout-left views on the C++ boundary.
 *
 * @param layers    Number of resolved vertical layers (km)
 * @param columns   Number of horizontal grid columns (im)
 * @param dt        Physics timestep (s)
 * @param tte_edmf  Total Turbulent Energy flag (0 = TKE, 1 = TTE)
 * @param t_lay     Temperature array [columns * layers] (K)
 * @param p_lay     Mean layer pressure [columns * layers] (Pa)
 * @param rho       Air density [columns * layers] (kg/m3)
 * @param u_wind    Zonal wind [columns * layers] (m/s)
 * @param v_wind    Meridional wind [columns * layers] (m/s)
 * @param q_vap     Water vapor specific humidity [columns * layers] (kg/kg)
 * @param te        TKE/TTE concentration [columns * layers] (m2/s2)
 * @param heat      Kinematic upward surface sensible heat flux [columns] (K m/s)
 * @param evap      Kinematic upward surface latent heat flux [columns] (kg/kg m/s)
 * @param stress    Surface wind stress [columns] (m2/s2)
 * @param do_canopy Forest canopy active flag (0 = OFF, 1 = ON)
 * @param cfch      Canopy forest height [columns] (m)
 * @param dt_temp   OUT: Temperature tendency [columns * layers] (K/s)
 * @param dt_u      OUT: U-wind tendency [columns * layers] (m/s2)
 * @param dt_v      OUT: V-wind tendency [columns * layers] (m/s2)
 * @param dt_te     OUT: TKE/TTE tendency [columns * layers] (m2/s3)
 */
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
);

#ifdef __cplusplus
}
#endif

#endif // SATMEDMF_INTERFACE_HPP
