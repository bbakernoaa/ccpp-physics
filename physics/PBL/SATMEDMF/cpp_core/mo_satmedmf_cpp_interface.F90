module mo_satmedmf_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_satmedmf_run(layers, columns, dt, tte_edmf, t_lay, p_lay, rho, &
        u_wind, v_wind, q_vap, te, heat, evap, stress, do_canopy, cfch, &
        dt_temp, dt_u, dt_v, dt_te) &
        bind(C, name="c_satmedmf_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: layers
      integer(c_size_t), value :: columns
      real(c_double), value    :: dt
      integer(c_int), value    :: tte_edmf
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: rho(*)
      real(c_double), intent(in)    :: u_wind(*)
      real(c_double), intent(in)    :: v_wind(*)
      real(c_double), intent(in)    :: q_vap(*)
      real(c_double), intent(in)    :: te(*)
      real(c_double), intent(in)    :: heat(*)
      real(c_double), intent(in)    :: evap(*)
      real(c_double), intent(in)    :: stress(*)
      integer(c_int), value    :: do_canopy
      real(c_double), intent(in)    :: cfch(*)
      real(c_double), intent(out)   :: dt_temp(*)
      real(c_double), intent(out)   :: dt_u(*)
      real(c_double), intent(out)   :: dt_v(*)
      real(c_double), intent(out)   :: dt_te(*)
    end subroutine c_satmedmf_run
  end interface

end module mo_satmedmf_cpp_interface
