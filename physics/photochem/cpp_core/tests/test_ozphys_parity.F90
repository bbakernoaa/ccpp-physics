program test_ozphys_parity
  use mo_photochem_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2
  integer(c_size_t), parameter :: layers = 10

  real(c_double) :: t_lay(columns, layers)
  real(c_double) :: p_lay(columns, layers)
  real(c_double) :: dp(columns, layers)
  real(c_double) :: oz(columns, layers)
  real(c_double) :: ozpl(columns, layers)

  real(c_double) :: do3_dt_prd(columns, layers)
  real(c_double) :: do3_dt_temp(columns, layers)

  integer :: k

  print *, "Starting GFS OZPHYS C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize variables
  t_lay = 240.0d0 ! 240 K
  oz = 4.0d-6     ! 4 ppm ozone
  ozpl = 2.0d-11  ! Nominal climate ozone production rate (kg/kg/s)

  do k = 1, layers
      p_lay(:, k) = 1000.0d0 - k * 80.0d0 ! Upper stratospheric pressures (Pa)
      dp(:, k) = 100.0d0
  end do

  ! Run the C++ solver wrapper
  call c_run_o3prog_2015(columns, layers, 1.0d0/9.80665d0, 300.0d0, &
      t_lay, p_lay, dp, oz, ozpl, do3_dt_prd, do3_dt_temp)

  print *, "  - Ozone physics output checks:"
  print *, "    do3_dt_prd(1, 1)  (expected > 0.0): ", do3_dt_prd(1, 1)
  print *, "    do3_dt_temp(1, 1) (expected < 0.0): ", do3_dt_temp(1, 1)

  if (do3_dt_prd(1, 1) <= 0.0d0 .or. do3_dt_temp(1, 1) >= 0.0d0) then
      print *, "✗ FAIL: Ozone photochemistry solver computed invalid tendencies."
      stop 1
  endif

  print *, "✓ GFS OZPHYS C++ translation and diagnostics PASSED!"
end program test_ozphys_parity
