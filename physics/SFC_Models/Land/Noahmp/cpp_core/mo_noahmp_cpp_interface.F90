module mo_noahmp_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function c_noahmp_sflx_run(columns, soil_layers, dt, &
        stc, smc, sh2o, sldpst, tg, tv, &
        sfctmp, sfcprs, q2, soldn, lwdn, &
        wind, ch, is_glacier, &
        idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc, &
        iopt_frz, iopt_inf, iopt_rad, iopt_alb, iopt_snf, &
        iopt_tbot, iopt_stc, iopt_trs, iopt_diag, iopt_rsf, &
        iopt_soil, iopt_pedo, iopt_crop, iopt_gla, iopt_z0m, &
        sheat, eta, gflux, runoff) result(status) &
        bind(C, name="c_noahmp_sflx_run")
      import :: c_size_t, c_double, c_int
      integer(c_int)                :: status
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: soil_layers
      real(c_double), value         :: dt
      real(c_double), intent(in)    :: stc(*)
      real(c_double), intent(in)    :: smc(*)
      real(c_double), intent(in)    :: sh2o(*)
      real(c_double), intent(in)    :: sldpst(*)
      real(c_double), intent(in)    :: tg(*)
      real(c_double), intent(in)    :: tv(*)
      real(c_double), intent(in)    :: sfctmp(*)
      real(c_double), intent(in)    :: sfcprs(*)
      real(c_double), intent(in)    :: q2(*)
      real(c_double), intent(in)    :: soldn(*)
      real(c_double), intent(in)    :: lwdn(*)
      real(c_double), intent(in)    :: wind(*)
      real(c_double), intent(in)    :: ch(*)
      integer(c_int), intent(in)    :: is_glacier(*)
      
      ! 19 physics options passed individually to keep interface flat
      integer(c_int), value         :: idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc
      integer(c_int), value         :: iopt_frz, iopt_inf, iopt_rad, iopt_alb, iopt_snf
      integer(c_int), value         :: iopt_tbot, iopt_stc, iopt_trs, iopt_diag, iopt_rsf
      integer(c_int), value         :: iopt_soil, iopt_pedo, iopt_crop, iopt_gla, iopt_z0m
      
      real(c_double), intent(out)   :: sheat(*)
      real(c_double), intent(out)   :: eta(*)
      real(c_double), intent(out)   :: gflux(*)
      real(c_double), intent(out)   :: runoff(*)
    end function c_noahmp_sflx_run
  end interface

end module mo_noahmp_cpp_interface
