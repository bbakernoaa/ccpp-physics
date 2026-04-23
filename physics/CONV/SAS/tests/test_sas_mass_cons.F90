program test_sas_mass_cons
  use machine, only: kind_phys
  use sascnvn, only: sascnvn_run, sascnvn_init
  use shalcnv, only: shalcnv_run, shalcnv_init
  implicit none

  ! Constants (simplified for test)
  real(kind_phys), parameter :: grav = 9.80665_kind_phys
  real(kind_phys), parameter :: cp = 1004.6_kind_phys
  real(kind_phys), parameter :: hvap = 2.501e6_kind_phys
  real(kind_phys), parameter :: rv = 461.5_kind_phys
  real(kind_phys), parameter :: fv = 0.608_kind_phys
  real(kind_phys), parameter :: t0c = 273.15_kind_phys
  real(kind_phys), parameter :: rd = 287.05_kind_phys
  real(kind_phys), parameter :: cvap = 1846.0_kind_phys
  real(kind_phys), parameter :: cliq = 4185.5_kind_phys
  real(kind_phys), parameter :: eps = rd/rv
  real(kind_phys), parameter :: epsm1 = eps - 1.0_kind_phys

  integer, parameter :: im = 1, km = 40
  integer :: jcap = 126
  real(kind_phys) :: delt = 450.0_kind_phys
  real(kind_phys) :: delp(im,km), prslp(im,km), psp(im), phil(im,km)
  real(kind_phys) :: qlc(im,km), qli(im,km), q1(im,km), t1(im,km), u1(im,km), v1(im,km)
  real(kind_phys) :: cldwrk(im), rn(im)
  integer :: kbot(im), ktop(im), kcnv(im)
  integer :: islimsk(im) = 1 ! land
  real(kind_phys) :: dot(im,km)
  integer :: ncloud = 1
  real(kind_phys) :: ud_mf(im,km), dd_mf(im,km), dt_mf(im,km), cnvw(im,km), cnvc(im,km)

  ! Shallow specific
  real(kind_phys) :: hpbl(im), heat(im), evap(im)

  ! Parameters
  real(kind_phys) :: clam = 0.1, c0 = 0.002, c1 = 0.002, pgcon = 0.55
  real(kind_phys) :: betal = 0.15, betas = 0.15, evfact = 0.3, evfactl = 0.3
  character(len=128) :: errmsg
  integer :: errflg

  integer :: i, k
  real(kind_phys) :: m_init, m_final, rain_mass

  ! Initialize
  call sascnvn_init(1, 1, errmsg, errflg)
  call shalcnv_init(.false., .true., 1, 1, errmsg, errflg)

  ! Initialize profiles - make it extremely conducive to convection
  psp(1) = 101325.0_kind_phys
  do k = 1, km
    delp(1,k) = psp(1) / km
    prslp(1,k) = psp(1) - (k-0.5)*delp(1,k)
    phil(1,k) = grav * (k-0.5) * 300.0_kind_phys
    t1(1,k) = 315.0 - (k-1)*6.0 ! Very steep lapse rate
    q1(1,k) = 0.025 * exp(-float(k-1)/3.0) ! Very moist bottom
    u1(1,k) = 10.0
    v1(1,k) = 5.0
    qlc(1,k) = 0.0
    qli(1,k) = 0.0
    dot(1,k) = -5.0 ! Strong upward motion
  enddo

  hpbl(1) = 2000.0
  heat(1) = 0.5
  evap(1) = 0.001

  ! --- TEST DEEP CONVECTION ---
  print *, "Testing Deep Convection (sascnvn)..."
  kcnv(1) = 0
  rn(1) = 0.0
  m_init = 0.0
  do k = 1, km
    m_init = m_init + (q1(1,k) + qlc(1,k) + qli(1,k)) * delp(1,k) / grav
  enddo

  call sascnvn_run( &
     grav,cp,hvap,rv,fv,t0c,rd,cvap,cliq,eps,epsm1, &
     im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli, &
     q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,islimsk, &
     dot,ncloud,ud_mf,dd_mf,dt_mf,cnvw,cnvc, &
     mp_phys=0, mp_phys_mg=1, clam=clam, c0=c0, c1=c1, &
     betal=betal, betas=betas, evfact=evfact, evfactl=evfactl, pgcon=pgcon, &
     errmsg=errmsg, errflg=errflg)

  m_final = 0.0
  do k = 1, km
    m_final = m_final + (q1(1,k) + qlc(1,k) + qli(1,k)) * delp(1,k) / grav
  enddo
  rain_mass = rn(1) * 1000.0

  print *, "  Initial tracer mass:", m_init
  print *, "  Final tracer mass:  ", m_final
  print *, "  Rain mass:         ", rain_mass
  print *, "  Mass Balance Error:", m_final + rain_mass - m_init
  print *, "  Convection Triggered (kcnv):", kcnv(1)
  print *, "  Rain amount (rn):", rn(1)

  ! --- TEST SHALLOW CONVECTION ---
  ! Reset q1, qlc, qli
  do k = 1, km
    q1(1,k) = 0.025 * exp(-float(k-1)/3.0)
    qlc(1,k) = 0.0
    qli(1,k) = 0.0
    kcnv(1) = 0
    rn(1) = 0.0
  enddo

  print *, "Testing Shallow Convection (shalcnv)..."
  m_init = 0.0
  do k = 1, km
    m_init = m_init + (q1(1,k) + qlc(1,k) + qli(1,k)) * delp(1,k) / grav
  enddo

  call shalcnv_run( &
     grav,cp,hvap,rv,fv,t0c,rd,cvap,cliq,eps,epsm1, &
     im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli, &
     q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk, &
     dot,ncloud,hpbl,heat,evap,ud_mf,dt_mf,cnvw,cnvc, &
     clam=clam, c0=c0, c1=c1, pgcon=pgcon, &
     errmsg=errmsg, errflg=errflg)

  m_final = 0.0
  do k = 1, km
    m_final = m_final + (q1(1,k) + qlc(1,k) + qli(1,k)) * delp(1,k) / grav
  enddo
  rain_mass = rn(1) * 1000.0

  print *, "  Initial tracer mass:", m_init
  print *, "  Final tracer mass:  ", m_final
  print *, "  Rain mass:         ", rain_mass
  print *, "  Mass Balance Error:", m_final + rain_mass - m_init
  print *, "  Convection Triggered (kcnv):", kcnv(1)
  print *, "  Rain amount (rn):", rn(1)

end program test_sas_mass_cons
