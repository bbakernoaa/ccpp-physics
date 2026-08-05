module mo_photochem_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_run_o3prog_2015(columns, layers, con_1ovg, dt, &
        t_lay, p_lay, dp, oz, ozpl, do3_dt_prd, do3_dt_temp) &
        bind(C, name="c_run_o3prog_2015")
      import :: c_size_t, c_double
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      real(c_double), value         :: con_1ovg
      real(c_double), value         :: dt
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: dp(*)
      real(c_double), intent(in)    :: oz(*)
      real(c_double), intent(in)    :: ozpl(*)
      real(c_double), intent(out)   :: do3_dt_prd(*)
      real(c_double), intent(out)   :: do3_dt_temp(*)
    end subroutine c_run_o3prog_2015

    subroutine c_run_h2ophys(columns, layers, dt, &
        t_lay, p_lay, dp, h2o, h2opltc, dqv_dt_prd, dqv_dt_qv) &
        bind(C, name="c_run_h2ophys")
      import :: c_size_t, c_double
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      real(c_double), value         :: dt
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: dp(*)
      real(c_double), intent(in)    :: h2o(*)
      real(c_double), intent(in)    :: h2opltc(*)
      real(c_double), intent(out)   :: dqv_dt_prd(*)
      real(c_double), intent(out)   :: dqv_dt_qv(*)
    end subroutine c_run_h2ophys
  end interface

end module mo_photochem_cpp_interface
