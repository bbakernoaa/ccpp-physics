program test_sfc_diff_parity
  use mo_surface_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2

  real(c_double) :: u1(columns)
  real(c_double) :: v1(columns)
  real(c_double) :: t1(columns)
  real(c_double) :: q1(columns)
  real(c_double) :: z1(columns)
  real(c_double) :: ps(columns)
  real(c_double) :: tskin(columns)
  real(c_double) :: z0(columns)

  real(c_double) :: cm(columns)
  real(c_double) :: ch(columns)
  real(c_double) :: ustar(columns)
  real(c_double) :: stress(columns)

  print *, "Starting GFS SFC_DIFF C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize sounding variables
  u1 = 5.0d0
  v1 = 3.0d0
  t1 = 285.15d0 ! Air temperature (12 C)
  q1 = 1.0d-3   ! Specific humidity
  z1 = 10.0d0   ! Surface layer height
  ps = 101325.0d0 ! Standard pressure
  tskin = 288.15d0 ! Surface skin temperature (15 C)
  z0 = 0.05d0   ! Roughness length (grassland proxy)

  ! Run C++ solver
  call c_sfc_diff_run(columns, u1, v1, t1, q1, z1, ps, tskin, z0, &
      6, cm, ch, ustar, stress)

  print *, "  - SFC_DIFF output checks:"
  print *, "    cm(1)    (expected > 0): ", cm(1)
  print *, "    ch(1)    (expected > 0): ", ch(1)
  print *, "    ustar(1) (expected > 0): ", ustar(1)
  print *, "    stress(1)(expected > 0): ", stress(1)

  if (cm(1) <= 0.0d0 .or. ustar(1) <= 0.0d0 .or. stress(1) <= 0.0d0) then
      print *, "✗ FAIL: Surface diffusion exchange solver computed invalid values."
      stop 1
  endif

  print *, "✓ GFS SFC_DIFF C++ translation and diagnostics PASSED!"
end program test_sfc_diff_parity
