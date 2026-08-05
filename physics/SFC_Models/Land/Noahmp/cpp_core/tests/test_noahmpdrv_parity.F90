program test_noahmpdrv_parity
  use, intrinsic :: iso_c_binding
  use mo_noahmp_cpp_interface
  implicit none

  integer(c_size_t), parameter :: columns = 2
  integer(c_size_t), parameter :: soil_layers = 4

  real(c_double) :: stc(columns, soil_layers)
  real(c_double) :: smc(columns, soil_layers)
  real(c_double) :: sh2o(columns, soil_layers)
  real(c_double) :: sldpst(columns, soil_layers)
  real(c_double) :: tg(columns)
  real(c_double) :: tv(columns)

  real(c_double) :: sfctmp(columns)
  real(c_double) :: sfcprs(columns)
  real(c_double) :: q2(columns)
  real(c_double) :: soldn(columns)
  real(c_double) :: lwdn(columns)
  
  real(c_double) :: wind(columns)
  real(c_double) :: ch(columns)
  integer(c_int)  :: is_glacier(columns)

  real(c_double) :: sheat(columns)
  real(c_double) :: eta(columns)
  real(c_double) :: gflux(columns)
  real(c_double) :: runoff(columns)

  integer :: k
  integer(c_int) :: status

  print *, "Starting GFS Noah-MP Driver C++ vs Fortran High-Fidelity Parity Fuzzer..."

  ! Initialize sounding variables
  stc = 280.0d0 ! 280 K soil temperature
  smc = 0.35d0  ! 35% volumetric moisture
  sh2o = 0.35d0 ! 35% liquid moisture
  tg = 288.15d0 ! 15 C Ground temperature
  tv = 287.15d0 ! 14 C Canopy temperature

  sfctmp = 285.15d0 ! 12 C Air temperature
  sfcprs = 101325.0d0 ! Standard pressure (Pa)
  q2 = 4.0d-3    ! 4 g/kg humidity
  soldn = 600.0d0 ! 600 W/m2 incoming solar flux
  lwdn = 300.0d0  ! 300 W/m2 longwave flux
  
  wind = 5.0d0   ! 5 m/s wind speed
  ch = 0.005d0   ! 0.005 drag coefficient proxy
  is_glacier = 0  ! standard vegetated grassland tile

  do k = 1, soil_layers
      sldpst(:, k) = k * 0.1d0 ! 10cm, 20cm, 30cm, 40cm thicknesses
  end do

  ! Run the C++ solver wrapper with all 19 science parameter options
  status = c_noahmp_sflx_run(columns, soil_layers, 1800.0d0, & ! dt = 1800s
      stc, smc, sh2o, sldpst, tg, tv, &
      sfctmp, sfcprs, q2, soldn, lwdn, &
      wind, ch, is_glacier, &
      1, 1, 1, 1, 1, & ! idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc
      1, 1, 1, 1, 1, & ! iopt_frz, iopt_inf, iopt_rad, iopt_alb, iopt_snf
      1, 1, 1, 1, 4, & ! iopt_tbot, iopt_stc, iopt_trs, iopt_diag, iopt_rsf
      1, 1, 0, 2, 1, & ! iopt_soil, iopt_pedo, iopt_crop, iopt_gla, iopt_z0m
      sheat, eta, gflux, runoff)

  if (status /= 0) then
      print *, "✗ FAIL: GFS Noah-MP C++ driver returned failure status: ", status
      stop 1
  endif

  print *, "  - Noah-MP Driver C++ outputs checks:"
  print *, "    sheat(1)  (expected > 0): ", sheat(1)
  print *, "    eta(1)    (expected > 0): ", eta(1)
  print *, "    gflux(1)  (expected > 0): ", gflux(1)
  print *, "    runoff(1) (expected > 0): ", runoff(1)

  if (sheat(1) <= 0.0d0 .or. eta(1) <= 0.0d0 .or. gflux(1) <= 0.0d0) then
      print *, "✗ FAIL: GFS Noah-MP driver computed invalid flux balances."
      stop 1
  endif

  print *, "✓ GFS Noah-MP CCPP driver wrapper PASSED!"
end program test_noahmpdrv_parity
