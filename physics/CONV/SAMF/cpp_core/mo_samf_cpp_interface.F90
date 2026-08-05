module mo_samf_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_samf_deep_convection_run(columns, layers, dt, &
        t_lay, q_vap, u_wind, v_wind, p_lay, p_int, z_lay, z_int, &
        dt_t, dt_q, dt_u, dt_v, &
        ud_mf, dd_mf, dt_mf, cnvw, cnvc, &
        kbot, ktop, kcnv, rain) &
        bind(C, name="c_samf_deep_convection_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      real(c_double), value    :: dt
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: q_vap(*)
      real(c_double), intent(in)    :: u_wind(*)
      real(c_double), intent(in)    :: v_wind(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: p_int(*)
      real(c_double), intent(in)    :: z_lay(*)
      real(c_double), intent(in)    :: z_int(*)
      real(c_double), intent(out)   :: dt_t(*)
      real(c_double), intent(out)   :: dt_q(*)
      real(c_double), intent(out)   :: dt_u(*)
      real(c_double), intent(out)   :: dt_v(*)
      real(c_double), intent(out)   :: ud_mf(*)
      real(c_double), intent(out)   :: dd_mf(*)
      real(c_double), intent(out)   :: dt_mf(*)
      real(c_double), intent(out)   :: cnvw(*)
      real(c_double), intent(out)   :: cnvc(*)
      integer(c_int), intent(out)   :: kbot(*)
      integer(c_int), intent(out)   :: ktop(*)
      integer(c_int), intent(out)   :: kcnv(*)
      real(c_double), intent(out)   :: rain(*)
    end subroutine c_samf_deep_convection_run

    subroutine c_samf_shallow_convection_run(columns, layers, dt, &
        t_lay, q_vap, u_wind, v_wind, p_lay, p_int, z_lay, z_int, &
        dt_t, dt_q, dt_u, dt_v, &
        ud_mf, dt_mf, kbot, ktop) &
        bind(C, name="c_samf_shallow_convection_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      real(c_double), value    :: dt
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: q_vap(*)
      real(c_double), intent(in)    :: u_wind(*)
      real(c_double), intent(in)    :: v_wind(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: p_int(*)
      real(c_double), intent(in)    :: z_lay(*)
      real(c_double), intent(in)    :: z_int(*)
      real(c_double), intent(out)   :: dt_t(*)
      real(c_double), intent(out)   :: dt_q(*)
      real(c_double), intent(out)   :: dt_u(*)
      real(c_double), intent(out)   :: dt_v(*)
      real(c_double), intent(out)   :: ud_mf(*)
      real(c_double), intent(out)   :: dt_mf(*)
      integer(c_int), intent(out)   :: kbot(*)
      integer(c_int), intent(out)   :: ktop(*)
    end subroutine c_samf_shallow_convection_run
  end interface

end module mo_samf_cpp_interface
