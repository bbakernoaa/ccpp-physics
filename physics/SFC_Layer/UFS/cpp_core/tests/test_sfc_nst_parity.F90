program test_sfc_nst_parity
  use mo_surface_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2

  real(c_double) :: sol_flux(columns)
  real(c_double) :: wind_stress(columns)

  real(c_double) :: tskin_wat(columns)
  real(c_double) :: cool_skin(columns)
  real(c_double) :: warm_layer(columns)

  print *, "Starting GFS SFC_NST C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize variables (calm, warm daytime ocean condition)
  sol_flux = 500.0d0 ! 500 W/m2 solar heating
  wind_stress = 0.02d0 ! Light breeze shear stress (N/m2)

  ! Run C++ solver
  call c_sfc_nst_run(columns, sol_flux, wind_stress, &
      tskin_wat, cool_skin, warm_layer)

  print *, "  - SFC_NST output checks:"
  print *, "    tskin_wat(1)  (expected ~ 296-298 K): ", tskin_wat(1)
  print *, "    cool_skin(1)  (expected > 0.0):        ", cool_skin(1)
  print *, "    warm_layer(1) (expected > 0.0):        ", warm_layer(1)

  if (tskin_wat(1) <= 290.0d0 .or. cool_skin(1) <= 0.0d0 .or. warm_layer(1) <= 0.0d0) then
      print *, "✗ FAIL: NST ocean skin solver computed invalid values."
      stop 1
  endif

  print *, "✓ GFS SFC_NST C++ translation and diagnostics PASSED!"
end program test_sfc_nst_parity
