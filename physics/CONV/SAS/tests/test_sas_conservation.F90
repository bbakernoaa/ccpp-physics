program test_sas_conservation
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  use sascnvn, only: sascnvn_run
  use shalcnv, only: shalcnv_run
  implicit none

  integer, parameter :: im = 1, km = 10
  real(kind=kind_phys) :: grav_val, cp, hvap, rv, fv, t0c, rgas, cvap, cliq, eps, epsm1, rd
  integer :: jcap, ncloud, mp_phys, mp_phys_mg
  real(kind=kind_phys) :: delt, psp(im), delp(im, km), prslp(im, km), phil(im, km)
  real(kind=kind_phys) :: qlc(im, km), qli(im, km), q1(im, km), t1(im, km), u1(im, km), v1(im, km)
  real(kind=kind_phys) :: dot(im, km)
  real(kind=kind_phys) :: cldwrk(im), rn(im)
  integer :: kbot(im), ktop(im), kcnv(im), islimsk(im)
  real(kind=kind_phys) :: ud_mf(im, km), dd_mf(im, km), dt_mf(im, km), cnvw(im, km), cnvc(im, km)

  ! optional arguments for sascnvn_run
  real(kind=kind_phys) :: qlcn(im,km), qicn(im,km), w_upi(im,km), cf_upi(im,km), cnv_mfd(im,km)
  real(kind=kind_phys) :: cnv_dqldt(im,km), clcn(im,km), cnv_fice(im,km), cnv_ndrop(im,km), cnv_nice(im,km)
  real(kind=kind_phys) :: clam, c0, c1, betal, betas, evfact, evfactl, pgcon

  ! optional arguments for shalcnv_run
  real(kind=kind_phys) :: hpbl(im), heat(im), evap(im)

  character(len=128) :: errmsg
  integer :: errflg
  integer :: i, k
  real(kind=kind_phys) :: m_init_q(im), m_final_q(im)
  logical :: mask(im)

  ! Initialize constants (approximate values)
  grav_val = 9.80665_kind_phys
  cp = 1004.6_kind_phys
  hvap = 2.5e6_kind_phys
  rv = 461.5_kind_phys
  fv = 0.608_kind_phys
  t0c = 273.15_kind_phys
  rgas = 287.05_kind_phys
  rd = 287.05_kind_phys
  cvap = 1846.0_kind_phys
  cliq = 4185.5_kind_phys
  eps = 0.622_kind_phys
  epsm1 = -0.378_kind_phys

  jcap = 0
  delt = 1200.0_kind_phys
  ncloud = 1
  mp_phys = 1
  mp_phys_mg = 0

  clam = 0.1
  c0 = 0.002
  c1 = 0.002
  betal = 0.05
  betas = 0.05
  evfact = 0.5
  evfactl = 0.5
  pgcon = 0.1

  ! Initialize state
  do k = 1, km
    delp(1, k) = 10000.0_kind_phys
    prslp(1, k) = 100000.0_kind_phys - (k-0.5)*10000.0_kind_phys
    phil(1, k) = (k-0.5)*1000.0_kind_phys
    q1(1, k) = 1.0e-3_kind_phys
    t1(1, k) = 300.0_kind_phys - (k-0.5)*5.0_kind_phys
    u1(1, k) = 10.0_kind_phys
    v1(1, k) = 0.0_kind_phys
    qlc(1, k) = 1.0e-5_kind_phys
    qli(1, k) = 1.0e-6_kind_phys
    dot(1, k) = -0.01_kind_phys
  enddo
  ! Force instability at bottom
  t1(1, 1) = 320.0_kind_phys
  q1(1, 1) = 2.0e-2_kind_phys

  psp(1) = 101325.0_kind_phys
  islimsk(1) = 0

  ! Initialize outputs/inouts
  rn = 0.0_kind_phys
  kbot = 0
  ktop = 0
  kcnv = 0
  ud_mf = 0.0_kind_phys
  dd_mf = 0.0_kind_phys
  dt_mf = 0.0_kind_phys
  cnvw = 0.0_kind_phys
  cnvc = 0.0_kind_phys

  qlcn = 0.
  qicn = 0.
  w_upi = 0.
  cf_upi = 0.
  cnv_mfd = 0.
  cnv_dqldt = 0.
  clcn = 0.
  cnv_fice = 0.
  cnv_ndrop = 0.
  cnv_nice = 0.

  hpbl = 2000.
  heat = 0.1
  evap = 1.0e-4

  mask = .true.

  write(*,*) "Testing Deep Convection (SAS)..."
  call calc_column_mass(im, km, mask, q1, delp, m_init_q)

  call sascnvn_run( &
       grav_val,cp,hvap,rv,fv,t0c,rgas,cvap,cliq,eps,epsm1,             &
       im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli,                 &
       q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,islimsk,                &
       dot,ncloud,ud_mf,dd_mf,dt_mf,cnvw,cnvc,                      &
       qlcn,qicn,w_upi,cf_upi,cnv_mfd,                              &
       cnv_dqldt,clcn,cnv_fice,cnv_ndrop,cnv_nice,mp_phys,          &
       mp_phys_mg,clam,c0,c1,betal,betas,evfact,evfactl,pgcon,      &
       errmsg,errflg)

  call calc_column_mass(im, km, mask, q1, delp, m_final_q)
  write(*, '(A, F20.10, A, F20.10)') "q1: Init =", m_init_q(1), " Final =", m_final_q(1)
  write(*, '(A, E15.5)') "Delta+Sink =", m_final_q(1) - m_init_q(1) + rn(1)*1000.0_kind_phys
  write(*, '(A, F10.5, A, I3, A, I3)') "Precip (m) =", rn(1), " kbot =", kbot(1), " ktop =", ktop(1)

  write(*,*) "Testing Shallow Convection (SAS)..."
  ! Reset state for shallow
  do k = 1, km
    q1(1, k) = 1.0e-3_kind_phys
    t1(1, k) = 300.0_kind_phys - (k-0.5)*5.0_kind_phys
  enddo
  t1(1, 1) = 310.0_kind_phys
  q1(1, 1) = 1.5e-2_kind_phys
  qlc = 1.0e-5_kind_phys
  qli = 1.0e-6_kind_phys
  rn = 0.0_kind_phys
  kcnv = 0

  call calc_column_mass(im, km, mask, q1, delp, m_init_q)

  call shalcnv_run( &
       grav_val,cp,hvap,rv,fv,t0c,rd,cvap,cliq,eps,epsm1,             &
       im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli,                 &
       q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk,                       &
       dot,ncloud,hpbl,heat,evap,ud_mf,dt_mf,cnvw,cnvc,             &
       clam,c0,c1,pgcon,                                            &
       errmsg,errflg)

  call calc_column_mass(im, km, mask, q1, delp, m_final_q)
  write(*, '(A, F20.10, A, F20.10, A, E15.5)') "q1: Init =", m_init_q(1), " Final =", m_final_q(1), &
       " Delta =", m_final_q(1) - m_init_q(1)
  write(*, '(A, F10.5, A, I3, A, I3)') "Precip (m) =", rn(1), " kbot =", kbot(1), " ktop =", ktop(1)

end program test_sas_conservation
