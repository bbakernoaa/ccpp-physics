module mo_radiation_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_rrtmg_sw_radiation_run(columns, layers, bands, &
        t_lay, q_vap, o3_vap, cld_frac, p_lay, p_int, albedo, &
        cos_solar_zenith, sw_heating_rate, sw_flux_down, sw_flux_up) &
        bind(C, name="c_rrtmg_sw_radiation_run")
      import :: c_size_t, c_double
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      integer(c_size_t), value :: bands
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: q_vap(*)
      real(c_double), intent(in)    :: o3_vap(*)
      real(c_double), intent(in)    :: cld_frac(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: p_int(*)
      real(c_double), intent(in)    :: albedo(*)
      real(c_double), intent(in)    :: cos_solar_zenith(*)
      real(c_double), intent(out)   :: sw_heating_rate(*)
      real(c_double), intent(out)   :: sw_flux_down(*)
      real(c_double), intent(out)   :: sw_flux_up(*)
    end subroutine c_rrtmg_sw_radiation_run

    subroutine c_rrtmg_lw_radiation_run(columns, layers, bands, &
        t_lay, q_vap, o3_vap, cld_frac, p_lay, p_int, emissivity, &
        lw_heating_rate, lw_flux_down, lw_flux_up) &
        bind(C, name="c_rrtmg_lw_radiation_run")
      import :: c_size_t, c_double
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      integer(c_size_t), value :: bands
      real(c_double), intent(in)    :: t_lay(*)
      real(c_double), intent(in)    :: q_vap(*)
      real(c_double), intent(in)    :: o3_vap(*)
      real(c_double), intent(in)    :: cld_frac(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: p_int(*)
      real(c_double), intent(in)    :: emissivity(*)
      real(c_double), intent(out)   :: lw_heating_rate(*)
      real(c_double), intent(out)   :: lw_flux_down(*)
      real(c_double), intent(out)   :: lw_flux_up(*)
    end subroutine c_rrtmg_lw_radiation_run
  end interface

end module mo_radiation_cpp_interface
