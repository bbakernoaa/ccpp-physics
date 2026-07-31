program benchmark_satmedmf
  use, intrinsic :: iso_c_binding
  use mo_satmedmf_cpp_interface
  use satmedmfvdifq, only : satmedmfvdifq_run_fortran
  use omp_lib
  implicit none

  integer(c_size_t), parameter :: layers = 100
  integer(c_size_t), parameter :: columns = 2000 ! 2000 Columns x 100 Layers = 200,000 cells!
  real(c_double), parameter :: dt = 300.0d0

  ! C ABI inputs & outputs
  real(c_double), allocatable :: t_lay_c(:)
  real(c_double), allocatable :: p_lay_c(:)
  real(c_double), allocatable :: rho_c(:)
  real(c_double), allocatable :: u_wind_c(:)
  real(c_double), allocatable :: v_wind_c(:)
  real(c_double), allocatable :: q_vap_c(:)
  real(c_double), allocatable :: te_c(:)

  real(c_double), allocatable :: heat_c(:)
  real(c_double), allocatable :: evap_c(:)
  real(c_double), allocatable :: stress_c(:)
  real(c_double), allocatable :: cfch_c(:)

  real(c_double), allocatable :: c_dt_temp(:)
  real(c_double), allocatable :: c_dt_u(:)
  real(c_double), allocatable :: c_dt_v(:)
  real(c_double), allocatable :: c_dt_te(:)

  ! Fortran Solver Inputs & Outputs
  integer :: im, km, ntrac, ntcw, ntrw, ntiw, ntke, ntqv
  real(c_double) :: my_grav, my_pi, my_rd, my_cp, my_rv, my_hvap, my_hfus, my_fv, my_eps, my_epsm1
  real(c_double) :: delt, xkzm_m, xkzm_h, xkzm_s
  real(c_double) :: dspfac, bl_upfr, bl_dnfr
  real(c_double) :: rlmx, elmx
  logical :: sa3dtke, dspheat, tte_edmf, do_canopy, cplaqm, gen_tend, ldiag3d
  integer :: sfc_rlm, tc_pbl, use_lpt

  real(8), allocatable :: def_1(:,:), def_2(:,:), def_3(:,:)
  real(8), allocatable :: dku3d_h(:,:), dku3d_e(:,:)
  real(8), allocatable :: rtg(:,:,:), tkeh(:,:)
  real(8), allocatable :: u1(:,:), v1(:,:)
  real(8), allocatable :: usfco(:), vsfco(:)
  logical :: use_oceanuv
  real(8), allocatable :: t1(:,:), q1(:,:,:)
  real(8), allocatable :: swh(:,:), hlw(:,:)
  real(8), allocatable :: xmu(:), garea(:), zvfun(:), sigmaf(:)
  real(8), allocatable :: psk(:), rbsoil(:), zorl(:), tsea(:)
  real(8), allocatable :: u10m(:), v10m(:)
  real(8), allocatable :: fm(:), fh(:)
  real(8), allocatable :: evap(:), heat(:), stress(:), spd1(:)
  integer, allocatable :: kpbl(:)
  real(8), allocatable :: prsi(:,:), del(:,:)
  real(8), allocatable :: prsl(:,:), prslk(:,:)
  real(8), allocatable :: phii(:,:), phil(:,:)
  real(8), allocatable :: dusfc(:), dvsfc(:), dtsfc(:), dqsfc(:)
  real(8), allocatable :: hpbl(:)
  real(8), allocatable :: dkt(:,:), dku(:,:)
  integer, allocatable :: kinver(:)
  real(8), allocatable :: claie(:), cfch(:), cfrt(:), cclu(:), cpopu(:)
  real(8), allocatable :: dtend(:,:,:)
  integer, allocatable :: dtidx(:,:)
  integer :: index_of_temperature, index_of_x_wind, index_of_y_wind, index_of_process_pbl
  real(8), allocatable :: ten_t(:,:), ten_u(:,:), ten_v(:,:)
  character(len=128) :: errmsg
  integer :: errflg

  integer(kind=8) :: start_time, end_time, rate
  real(c_double) :: elapsed_cpp, elapsed_fortran, speedup
  integer :: i, k, idx, scenario, t
  real(c_double) :: max_diff, diff_temp, diff_u, diff_v, diff_te

  im = int(columns)
  km = int(layers)
  ntrac = 2
  ntcw = 1
  ntrw = 1
  ntiw = 1
  ntke = 2
  ntqv = 1

  my_grav = 9.80665d0
  my_pi = 3.141592653589793d0
  my_rd = 287.05d0
  my_cp = 1004.6d0
  my_rv = 461.5d0
  my_hvap = 2.5d6
  my_hfus = 3.33d5
  my_fv = 1.0d0
  my_eps = my_rd / my_rv
  my_epsm1 = my_eps - 1.0d0
  delt = dt

  xkzm_m = 1.0d0
  xkzm_h = 1.0d0
  xkzm_s = 1.0d0
  dspfac = 1.0d0
  bl_upfr = 1.0d0
  bl_dnfr = 1.0d0
  rlmx = 300.0d0
  elmx = 300.0d0
  sfc_rlm = 0
  tc_pbl = 1
  use_lpt = 0

  sa3dtke = .false.
  dspheat = .false.
  tte_edmf = .true.
  do_canopy = .true.
  cplaqm = .false.
  gen_tend = .false.
  ldiag3d = .false.
  use_oceanuv = .false.

  index_of_temperature = 1
  index_of_x_wind = 2
  index_of_y_wind = 3
  index_of_process_pbl = 4

  ! Allocate all Fortran parameters
  allocate(def_1(im, km), def_2(im, km), def_3(im, km))
  allocate(dku3d_h(im, km), dku3d_e(im, km))
  allocate(rtg(im, km, ntrac), tkeh(im, km))
  allocate(u1(im, km), v1(im, km))
  allocate(usfco(im), vsfco(im))
  allocate(t1(im, km), q1(im, km, ntrac))
  allocate(swh(im, km), hlw(im, km))
  allocate(xmu(im), garea(im), zvfun(im), sigmaf(im))
  allocate(psk(im), rbsoil(im), zorl(im), tsea(im))
  allocate(u10m(im), v10m(im))
  allocate(fm(im), fh(im))
  allocate(evap(im), heat(im), stress(im), spd1(im))
  allocate(kpbl(im))
  allocate(prsi(im, km+1), del(im, km))
  allocate(prsl(im, km), prslk(im, km))
  allocate(phii(im, km+1), phil(im, km))
  allocate(dusfc(im), dvsfc(im), dtsfc(im), dqsfc(im))
  allocate(hpbl(im))
  allocate(dkt(im, km), dku(im, km))
  allocate(kinver(im))
  allocate(claie(im), cfch(im), cfrt(im), cclu(im), cpopu(im))
  allocate(dtend(im, km, ntrac))
  allocate(dtidx(im, km))
  allocate(ten_t(im, km), ten_u(im, km), ten_v(im, km))

  ! Allocate C ABI arrays
  allocate(t_lay_c(columns * layers))
  allocate(p_lay_c(columns * layers))
  allocate(rho_c(columns * layers))
  allocate(u_wind_c(columns * layers))
  allocate(v_wind_c(columns * layers))
  allocate(q_vap_c(columns * layers))
  allocate(te_c(columns * layers))
  allocate(heat_c(columns))
  allocate(evap_c(columns))
  allocate(stress_c(columns))
  allocate(cfch_c(columns))
  allocate(c_dt_temp(columns * layers))
  allocate(c_dt_u(columns * layers))
  allocate(c_dt_v(columns * layers))
  allocate(c_dt_te(columns * layers))

  call system_clock(count_rate=rate)

  print *, "============================================================="
  print *, "   CCPP SATMEDMF 200,000-CELL HIGH-RESOLUTION BENCHMARK      "
  print *, "           Grid Dimensions: 2000 Columns x 100 Layers        "
  print *, "============================================================="

  heat_c = 0.05d0
  cfch_c = 100.0d0
  heat = heat_c
  cfch = cfch_c
  evap = 0.0001d0
  stress = 0.02d0

  ! Initialize inputs
  do i = 1, int(columns)
     do k = 1, int(layers)
        idx = i + int(columns) * (k - 1)
        t_lay_c(idx) = 288.0d0 - 0.1d0 * k
        te_c(idx) = 1.0d0 - 0.001d0 * k
        p_lay_c(idx) = 100000.0d0 - 100.0d0 * k
        rho_c(idx) = 1.2d0 - 0.001d0 * k
        u_wind_c(idx) = 5.0d0 + 0.02d0 * k
        v_wind_c(idx) = 2.0d0
        q_vap_c(idx) = 0.005d0

        ! Copy to Fortran reference matrices
        t1(i, k) = t_lay_c(idx)
        u1(i, k) = u_wind_c(idx)
        v1(i, k) = v_wind_c(idx)
        q1(i, k, 1) = q_vap_c(idx)
        rtg(i, k, 1) = te_c(idx)
        prsl(i, k) = p_lay_c(idx)
        prslk(i, k) = p_lay_c(idx)
        del(i, k) = 100.0d0
        phii(i, k) = 0.d0
        phil(i, k) = 0.d0
     enddo
     usfco(i) = 0.d0
     vsfco(i) = 0.d0
     xmu(i) = 1.d0
     garea(i) = 1.d0
     zvfun(i) = 1.d0
     sigmaf(i) = 1.d0
     psk(i) = 1.d0
     rbsoil(i) = 1.d0
     zorl(i) = 0.1d0
     u10m(i) = 1.d0
     v10m(i) = 1.d0
     fm(i) = 1.d0
     fh(i) = 1.d0
     spd1(i) = 1.d0
     kinver(i) = 0
  enddo

  ! Profile C++ OpenMP scaling across thread counts (1)
  do t = 1, 1
     if (t /= 1) cycle

     call omp_set_num_threads(t)

     call system_clock(count=start_time)
     call c_satmedmf_run(layers, columns, dt, merge(1, 0, tte_edmf), &
          t_lay_c, p_lay_c, rho_c, u_wind_c, v_wind_c, q_vap_c, te_c, heat_c, evap_c, stress_c, &
          merge(1, 0, do_canopy), cfch_c, c_dt_temp, c_dt_u, c_dt_v, c_dt_te)
     call system_clock(count=end_time)
     elapsed_cpp = real(end_time - start_time, c_double) / real(rate, c_double)

     print '(A, I2, A, F10.6, A)', "  C++ Threads: ", t, " | Execution Time: ", elapsed_cpp, "s"
  enddo

  print *, "-------------------------------------------------------------"

  ! Timing the ACTUAL original Fortran solver subroutine (on standard thread execution)
  print *, "  Running Native Fortran Physics Solver on 1 Thread..."
  call system_clock(count=start_time)
  call satmedmfvdifq_run_fortran(im, km, ntrac, ntcw, ntrw, ntiw, ntke, &
       my_grav, my_pi, my_rd, my_cp, my_rv, my_hvap, my_hfus, my_fv, my_eps, my_epsm1, def_1, def_2, def_3, &
       sa3dtke, dku3d_h, dku3d_e, rtg, u1, v1, t1, q1, usfco, vsfco, use_oceanuv, &
       swh, hlw, xmu, garea, zvfun, sigmaf, psk, rbsoil, zorl, u10m, v10m, fm, fh, &
       tsea, heat, evap, stress, spd1, kpbl, prsi, del, prsl, prslk, phii, phil, &
       delt, tte_edmf, dspheat, dusfc, dvsfc, dtsfc, dqsfc, hpbl, dkt, dku, tkeh, &
       kinver, xkzm_m, xkzm_h, xkzm_s, dspfac, bl_upfr, bl_dnfr, rlmx, elmx, &
       sfc_rlm, tc_pbl, use_lpt, do_canopy, cplaqm, claie, cfch, cfrt, cclu, cpopu, &
       ntqv, dtend, dtidx, index_of_temperature, index_of_x_wind, index_of_y_wind, &
       index_of_process_pbl, gen_tend, ldiag3d, ten_t, ten_u, ten_v, errmsg, errflg)
  call system_clock(count=end_time)
  elapsed_fortran = real(end_time - start_time, c_double) / real(rate, c_double)

  print '(A, F10.6, A)', "  Fortran Execution Time: ", elapsed_fortran, "s"
  print *, "============================================================="

  deallocate(t_lay_c, p_lay_c, rho_c, u_wind_c, v_wind_c, q_vap_c, te_c)
  deallocate(heat_c, evap_c, stress_c, cfch_c)
  deallocate(c_dt_temp, c_dt_u, c_dt_v, c_dt_te)
  deallocate(def_1, def_2, def_3, dku3d_h, dku3d_e)
  deallocate(rtg, tkeh, u1, v1, usfco, vsfco, t1, q1, swh, hlw)
  deallocate(xmu, garea, zvfun, sigmaf, psk, rbsoil, zorl, tsea)
  deallocate(u10m, v10m, fm, fh, evap, heat, stress, spd1, kpbl)
  deallocate(prsi, del, prsl, prslk, phii, phil, dusfc, dvsfc, dtsfc, dqsfc)
  deallocate(hpbl, dkt, dku, kinver)
  deallocate(claie, cfch, cfrt, cclu, cpopu, dtend, dtidx)
  deallocate(ten_t, ten_u, ten_v)

end program benchmark_satmedmf
