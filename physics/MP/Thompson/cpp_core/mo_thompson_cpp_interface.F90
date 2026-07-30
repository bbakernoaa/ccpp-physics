module mo_thompson_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_thompson_microphysics_run(layers, columns, dt, t_lay, p_lay, rho, &
        qv, qc, qr, qi, qs, qg, ni, nr, ns, ng, precip) &
        bind(C, name="c_thompson_microphysics_run")
      import :: c_size_t, c_double
      integer(c_size_t), value :: layers
      integer(c_size_t), value :: columns
      real(c_double), value    :: dt
      real(c_double), intent(inout) :: t_lay(*)
      real(c_double), intent(in)    :: p_lay(*)
      real(c_double), intent(in)    :: rho(*)
      real(c_double), intent(inout) :: qv(*)
      real(c_double), intent(inout) :: qc(*)
      real(c_double), intent(inout) :: qr(*)
      real(c_double), intent(inout) :: qi(*)
      real(c_double), intent(inout) :: qs(*)
      real(c_double), intent(inout) :: qg(*)
      real(c_double), intent(inout) :: ni(*)
      real(c_double), intent(inout) :: nr(*)
      real(c_double), intent(inout) :: ns(*)
      real(c_double), intent(inout) :: ng(*)
      real(c_double), intent(out)   :: precip(*)
    end subroutine c_thompson_microphysics_run
  end interface

end module mo_thompson_cpp_interface
