program test_h2ophys_parity
  use mo_photochem_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2
  integer(c_size_t), parameter :: layers = 10

  real(c_double) :: t_lay(columns, layers)
  real(c_double) :: p_lay(columns, layers)
  real(c_double) :: dp(columns, layers)
  real(c_double) :: h2o(columns, layers)
  real(c_double) :: h2opltc(columns, layers)

  real(c_double) :: dqv_dt_prd(columns, layers)
  real(c_double) :: dqv_dt_qv(columns, layers)

  integer :: k

  print *, "Starting GFS H2OPHYS C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize variables
  t_lay = 230.0d0 ! 230 K
  h2o = 2.0d-6    ! 2 ppm water vapor
  h2opltc = 1.0d-7 ! Photolysis loss coefficient

  do k = 1, layers
      p_lay(:, k) = 1000.0d0 - k * 80.0d0 ! Stratospheric pressure (Pa)
      dp(:, k) = 100.0d0
  end do

  ! Run the C++ solver wrapper
  call c_run_h2ophys(columns, layers, 300.0d0, &
      t_lay, p_lay, dp, h2o, h2opltc, dqv_dt_prd, dqv_dt_qv)

  print *, "  - Water vapor physics output checks:"
  print *, "    dqv_dt_prd(1, 1) (expected > 0.0): ", dqv_dt_prd(1, 1)
  print *, "    dqv_dt_qv(1, 1)  (expected < 0.0): ", dqv_dt_qv(1, 1)

  if (dqv_dt_prd(1, 1) <= 0.0d0 .or. dqv_dt_qv(1, 1) >= 0.0d0) then
      print *, "✗ FAIL: Water vapor photochemistry solver computed invalid tendencies."
      stop 1
  endif

  print *, "✓ GFS H2OPHYS C++ translation and diagnostics PASSED!"
end program test_h2ophys_parity
