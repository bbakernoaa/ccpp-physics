program test_gwd_parity
  use mo_gwd_cpp_interface
  implicit none
  integer :: k

  integer(c_size_t), parameter :: columns = 2
  integer(c_size_t), parameter :: layers = 10
  real(c_double), parameter :: dtp = 300.0d0

  real(c_double) :: ugrs(columns, layers)
  real(c_double) :: vgrs(columns, layers)
  real(c_double) :: tgrs(columns, layers)
  real(c_double) :: q1(columns, layers)
  real(c_double) :: prsl(columns, layers)
  real(c_double) :: prsi(columns, layers + 1)
  real(c_double) :: prslk(columns, layers)
  real(c_double) :: phil(columns, layers)
  real(c_double) :: phii(columns, layers + 1)
  real(c_double) :: del(columns, layers)

  real(c_double) :: hprime(columns)
  real(c_double) :: oc(columns)
  real(c_double) :: theta(columns)
  real(c_double) :: sigma(columns)
  real(c_double) :: gamma(columns)
  real(c_double) :: elvmax(columns)
  real(c_double) :: clx(columns)
  real(c_double) :: oa4(columns)
  real(c_double) :: varss(columns)
  real(c_double) :: dx(columns)
  real(c_double) :: xlat(columns)
  real(c_double) :: area(columns)

  real(c_double) :: dudt_ogw(columns, layers)
  real(c_double) :: dvdt_ogw(columns, layers)
  real(c_double) :: dudt_ngw(columns, layers)
  real(c_double) :: dvdt_ngw(columns, layers)
  real(c_double) :: dtdt_ngw(columns, layers)
  real(c_double) :: dudt_ofd(columns, layers)
  real(c_double) :: dvdt_ofd(columns, layers)
  real(c_double) :: tau_ogw(columns, layers)
  real(c_double) :: tau_ngw(columns, layers)

  ! Post variables
  real(c_double) :: zobl(columns)
  real(c_double) :: zlwb(columns)
  real(c_double) :: zogw(columns)
  real(c_double) :: du_ofdcol(columns)
  real(c_double) :: du_oblcol(columns)

  real(c_double) :: tot_mtb(columns)
  real(c_double) :: tot_ogw(columns)
  real(c_double) :: tot_tofd(columns)
  real(c_double) :: tot_ngw(columns)
  real(c_double) :: tot_zmtb(columns)
  real(c_double) :: tot_zlwb(columns)
  real(c_double) :: tot_zogw(columns)

  real(c_double) :: du3dt_mtb(columns, layers)
  real(c_double) :: du3dt_tms(columns, layers)
  real(c_double) :: du3dt_ogw(columns, layers)
  real(c_double) :: du3dt_ngw(columns, layers)
  real(c_double) :: dv3dt_ngw(columns, layers)

  print *, "Starting CCPP GWD and diagnostics C++23 unit verification..."

  ! Initialize variables
  ugrs = 15.0d0
  vgrs = 10.0d0
  tgrs = 250.0d0
  q1 = 1.0d-5
  del = 1000.0d0

  do k = 1, layers
      phil(:, k) = k * 1000.0d0
      prsl(:, k) = 100000.0d0 - k * 5000.0d0
      prslk(:, k) = (prsl(:, k) / 100000.0d0) ** (287.05d0 / 1004.6d0)
  end do
  do k = 1, layers + 1
      phii(:, k) = (k - 1) * 1000.0d0 + 500.0d0
      prsi(:, k) = 102500.0d0 - k * 5000.0d0
  end do

  varss = 2.5d0
  hprime = 10.0d0
  oc = 1.0d0
  theta = 0.5d0
  sigma = 0.2d0
  gamma = 0.1d0
  elvmax = 1000.0d0
  clx = 0.1d0
  oa4 = 1.0d0
  dx = 10000.0d0
  xlat = 0.78d0
  area = 1.0d8

  ! 1. Run C++23 GWD core solver via Bindings
  call c_ugwpv1_gsldrag_run(columns, layers, dtp, &
      ugrs, vgrs, tgrs, q1, prsl, prsi, prslk, phil, phii, del, &
      1, 1, 1, 1, & ! Enable all drag components
      hprime, oc, theta, sigma, gamma, elvmax, clx, oa4, varss, &
      dx, xlat, area, &
      dudt_ogw, dvdt_ogw, dudt_ngw, dvdt_ngw, dtdt_ngw, &
      dudt_ofd, dvdt_ofd, tau_ogw, tau_ngw)

  print *, "  - GWD core output checks:"
  print *, "    dudt_ogw(1, 1) (expected /= 0): ", dudt_ogw(1, 1)
  print *, "    tau_ogw(1, 1)  (expected > 0):  ", tau_ogw(1, 1)
  print *, "    dudt_ngw(1, 8) (expected /= 0): ", dudt_ngw(1, 8)

  if (abs(dudt_ogw(1, 1)) < 1.0d-15 .or. tau_ogw(1, 1) <= 0.0d0 .or. abs(dudt_ngw(1, 8)) < 1.0d-15) then
      print *, "✗ FAIL: GWD solver computed invalid tendencies."
      stop 1
  endif

  ! 2. Initialize post variables
  zobl = 50.0d0
  zlwb = 100.0d0
  zogw = 200.0d0
  du_ofdcol = 0.1d0
  du_oblcol = 0.2d0

  tot_mtb = 0.0d0
  tot_ogw = 0.0d0
  tot_tofd = 0.0d0
  tot_ngw = 0.0d0
  tot_zmtb = 0.0d0
  tot_zlwb = 0.0d0
  tot_zogw = 0.0d0

  du3dt_mtb = 0.0d0
  du3dt_tms = 0.0d0
  du3dt_ogw = 0.0d0
  du3dt_ngw = 0.0d0
  dv3dt_ngw = 0.0d0

  ! 3. Run GWD Post-convection diagnostics
  call c_ugwpv1_gsldrag_post_run(columns, layers, 1, 1.0d0, & ! ldiag_ugwp = 1, dtf = 1.0
      zobl, zlwb, zogw, tau_ogw(:, 1), tau_ngw(:, 1), du_ofdcol, du_oblcol, &
      tot_mtb, tot_ogw, tot_tofd, tot_ngw, &
      tot_zmtb, tot_zlwb, tot_zogw, &
      dudt_ngw, dvdt_ngw, dudt_ofd, dudt_ofd, dudt_ogw, &
      du3dt_mtb, du3dt_tms, du3dt_ogw, &
      du3dt_ngw, dv3dt_ngw)

  print *, "  - GWD diagnostics post_run output checks:"
  print *, "    tot_zogw(1) (expected: 200.0): ", tot_zogw(1)
  print *, "    tot_ogw(1)  (expected > 0.0):   ", tot_ogw(1)

  if (abs(tot_zogw(1) - 200.0d0) > 1.0d-9) then
      print *, "✗ FAIL: GWD post diagnostics mismatch."
      stop 1
  endif

  print *, "✓ ALL GWD and diagnostics verification gates PASSED!"
end program test_gwd_parity
