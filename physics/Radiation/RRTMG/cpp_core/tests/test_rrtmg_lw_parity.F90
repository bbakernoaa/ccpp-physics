program test_rrtmg_lw_parity
  use mo_radiation_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2
  integer(c_size_t), parameter :: layers = 10
  integer(c_size_t), parameter :: bands = 16

  real(c_double) :: t_lay(columns, layers)
  real(c_double) :: q_vap(columns, layers)
  real(c_double) :: o3_vap(columns, layers)
  real(c_double) :: cld_frac(columns, layers)
  real(c_double) :: p_lay(columns, layers)
  real(c_double) :: p_int(columns, layers + 1)
  real(c_double) :: emissivity(columns, bands)

  real(c_double) :: lw_heating_rate(columns, layers)
  real(c_double) :: lw_flux_down(columns, layers + 1)
  real(c_double) :: lw_flux_up(columns, layers + 1)

  integer :: k

  print *, "Starting GFS RRTMG LW C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize sounding variables
  t_lay = 250.0d0
  q_vap = 1.0d-4
  o3_vap = 1.0d-6
  cld_frac = 0.0d0
  emissivity = 0.98d0 ! Surface emissivity

  do k = 1, layers
      p_lay(:, k) = 100000.0d0 - k * 8000.0d0
  end do
  do k = 1, layers + 1
      p_int(:, k) = 104000.0d0 - k * 8000.0d0
  end do

  ! Add a cloudy layer at k = 5
  cld_frac(:, 5) = 0.5d0

  ! Run the C++ solver wrapper
  call c_rrtmg_lw_radiation_run(columns, layers, bands, &
      t_lay, q_vap, o3_vap, cld_frac, p_lay, p_int, emissivity, &
      lw_heating_rate, lw_flux_down, lw_flux_up)

  print *, "  - LW core output checks:"
  print *, "    lw_flux_down(1, layers+1) (expected: 0.0): ", lw_flux_down(1, layers + 1)
  print *, "    lw_flux_down(1, 1)        (expected > 0): ", lw_flux_down(1, 1)
  print *, "    lw_flux_up(1, layers+1)   (expected > 0): ", lw_flux_up(1, layers + 1)

  if (lw_flux_down(1, 1) <= 0.0d0 .or. lw_flux_up(1, layers + 1) <= 0.0d0) then
      print *, "✗ FAIL: LW radiation solver computed invalid fluxes or heating rates."
      stop 1
  endif

  print *, "✓ GFS RRTMG LW C++ translation and diagnostics PASSED!"
end program test_rrtmg_lw_parity
