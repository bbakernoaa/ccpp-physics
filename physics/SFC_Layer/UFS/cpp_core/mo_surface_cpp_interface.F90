module mo_surface_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_sfc_diff_run(columns, u1, v1, t1, q1, z1, ps, tskin, z0, &
        sfc_z0_type, cm, ch, ustar, stress) &
        bind(C, name="c_sfc_diff_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: columns
      real(c_double), intent(in)    :: u1(*)
      real(c_double), intent(in)    :: v1(*)
      real(c_double), intent(in)    :: t1(*)
      real(c_double), intent(in)    :: q1(*)
      real(c_double), intent(in)    :: z1(*)
      real(c_double), intent(in)    :: ps(*)
      real(c_double), intent(in)    :: tskin(*)
      real(c_double), intent(in)    :: z0(*)
      integer(c_int), value         :: sfc_z0_type
      real(c_double), intent(out)   :: cm(*)
      real(c_double), intent(out)   :: ch(*)
      real(c_double), intent(out)   :: ustar(*)
      real(c_double), intent(out)   :: stress(*)
    end subroutine c_sfc_diff_run

    subroutine c_sfc_nst_run(columns, sol_flux, wind_stress, &
        tskin_wat, cool_skin, warm_layer) &
        bind(C, name="c_sfc_nst_run")
      import :: c_size_t, c_double
      integer(c_size_t), value :: columns
      real(c_double), intent(in)    :: sol_flux(*)
      real(c_double), intent(in)    :: wind_stress(*)
      real(c_double), intent(out)   :: tskin_wat(*)
      real(c_double), intent(out)   :: cool_skin(*)
      real(c_double), intent(out)   :: warm_layer(*)
    end subroutine c_sfc_nst_run
  end interface

end module mo_surface_cpp_interface
