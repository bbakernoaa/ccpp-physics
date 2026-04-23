program test_sas_comprehensive
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  use sascnvn, only: sascnvn_run
  use shalcnv, only: shalcnv_run
  use funcphys, only: fpvs
  implicit none

  integer, parameter :: im = 1, km = 64
  real(kind=kind_phys) :: grav_val, cp, hvap, rv, fv, t0c, rgas, cvap, cliq, eps, epsm1, rd
  integer :: jcap, ncloud, mp_phys, mp_phys_mg
  real(kind=kind_phys) :: delt, psp(im), delp(im, km), prslp(im, km), phil(im, km)
  real(kind=kind_phys) :: qlc(im, km), qli(im, km), q1(im, km), t1(im, km), u1(im, km), v1(im, km)
  real(kind=kind_phys) :: dot(im, km)
  real(kind=kind_phys) :: cldwrk(im), rn(im)
  integer :: kbot(im), ktop(im), kcnv(im), islimsk(im)
  real(kind=kind_phys) :: ud_mf(im, km), dd_mf(im, km), dt_mf(im, km), cnvw(im, km), cnvc(im, km)

  real(kind=kind_phys) :: qlcn(im,km), qicn(im,km), w_upi(im,km), cf_upi(im,km), cnv_mfd(im,km)
  real(kind=kind_phys) :: cnv_dqldt(im,km), clcn(im,km), cnv_fice(im,km), cnv_ndrop(im,km), cnv_nice(im,km)
  real(kind=kind_phys) :: clam, c0, c1, betal, betas, evfact, evfactl, pgcon
  real(kind=kind_phys) :: hpbl(im), heat(im), evap(im)

  character(len=128) :: errmsg
  integer :: errflg
  integer :: opt, k
  real(kind=kind_phys) :: m_init_q(im), m_final_q(im)
  logical :: mask(im)
  character(len=20) :: method_name
  character(len=10) :: conv_status

  grav_val = 9.80665; cp = 1004.6; hvap = 2.5e6; rv = 461.5; fv = 0.608; t0c = 273.15
  rgas = 287.05; rd = 287.05; cvap = 1846.0; cliq = 4185.5; eps = 0.622; epsm1 = -0.378
  jcap = 0; delt = 1200.0; ncloud = 1; mp_phys = 1; mp_phys_mg = 0
  clam = 0.1; c0 = 0.002; c1 = 0.002; betal = 0.05; betas = 0.05
  evfact = 0.5; evfactl = 0.5; pgcon = 0.1
  mask = .true.

  write(*,*) "| Scheme | Transport Method | Initial Mass | Final Mass | Delta | Status |"
  write(*,*) "|:---|:---|---:|---:|---:|:---|"

  do opt = 1, 3
    if (opt == 1) method_name = "Flux-Form"
    if (opt == 2) method_name = "Hole-Filling"
    if (opt == 3) method_name = "Advective (Leaky)"

    ! --- Deep Convection Test ---
    call set_loaded_gun_state(q1, t1, u1, v1, qlc, qli, delp, prslp, phil, psp, dot, islimsk, im, km, eps, rd)
    call calc_column_mass(im, km, mask, q1, delp, m_init_q)
    kcnv = 0; rn = 0.0
    call sascnvn_run(grav_val,cp,hvap,rv,fv,t0c,rgas,cvap,cliq,eps,epsm1,im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli, &
                     q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,islimsk,dot,ncloud,ud_mf,dd_mf,dt_mf,cnvw,cnvc, &
                     qlcn,qicn,w_upi,cf_upi,cnv_mfd,cnv_dqldt,clcn,cnv_fice,cnv_ndrop,cnv_nice,mp_phys, &
                     mp_phys_mg,clam,c0,c1,betal,betas,evfact,evfactl,pgcon,errmsg,errflg)
    call calc_column_mass(im, km, mask, q1, delp, m_final_q)
    conv_status = "Inactive"
    if (any(rn > 0.0)) conv_status = "Active"
    write(*, '("| Deep | ", A15, " | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", A10, " |")') &
     & method_name, m_init_q(1), m_final_q(1), m_final_q(1) - m_init_q(1) + rn(1)*1000.0_kind_phys, conv_status

    ! --- Shallow Convection Test ---
    call set_loaded_gun_state(q1, t1, u1, v1, qlc, qli, delp, prslp, phil, psp, dot, islimsk, im, km, eps, rd)
    hpbl(1) = 2000.0; heat(1) = 0.5; evap(1) = 1.0e-3
    kcnv = 0; rn = 0.0
    call calc_column_mass(im, km, mask, q1, delp, m_init_q)
    call shalcnv_run(grav_val,cp,hvap,rv,fv,t0c,rd,cvap,cliq,eps,epsm1,im,km,jcap,delt,delp,prslp,psp,phil,qlc,qli, &
                     q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk,dot,ncloud,hpbl,heat,evap,ud_mf,dt_mf,cnvw,cnvc, &
                     clam,c0,c1,pgcon,errmsg,errflg)
    call calc_column_mass(im, km, mask, q1, delp, m_final_q)
    conv_status = "Inactive"
    if (any(rn > 0.0)) conv_status = "Active"
    write(*, '("| Shallow | ", A15, " | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", A10, " |")') &
     & method_name, m_init_q(1), m_final_q(1), m_final_q(1) - m_init_q(1), conv_status
  enddo

contains

  subroutine set_loaded_gun_state(q1, t1, u1, v1, qlc, qli, delp, prslp, phil, psp, dot, islimsk, im, km, eps, rd)
    integer, intent(in) :: im, km
    real(kind=kind_phys), intent(out) :: q1(im, km), t1(im, km), u1(im, km), v1(im, km), qlc(im, km), qli(im, km)
    real(kind=kind_phys), intent(out) :: delp(im, km), prslp(im, km), phil(im, km), psp(im), dot(im, km)
    integer, intent(out) :: islimsk(im)
    real(kind=kind_phys), intent(in) :: eps, rd
    real(kind=kind_phys) :: p_levels(6), t_levels(6), td_levels(6), u_levels(6), v_levels(6)
    real(kind=kind_phys) :: p_k, t_k, td_k, es_k, weight, p_km1
    integer :: k, l
    p_levels =  (/ 1000.0, 850.0, 700.0, 500.0, 300.0, 200.0 /)
    t_levels =  (/ 35.0,   25.0,  12.0,  -12.0, -35.0, -55.0 /) + 273.15
    td_levels = (/ 25.0,   24.0,  20.0,  10.0,  0.0,  -20.0 /) + 273.15
    u_levels = (/ -10.0, 0.0, 35.0, 60.0, 100.0, 110.0 /) * 0.51444
    v_levels = (/ 10.0, 40.0, 35.0, 25.0, 0.0, 0.0 /) * 0.51444
    psp(1) = 101325.0; islimsk(1) = 0; phil(1, 1) = 100.0 * 9.8
    do k = 1, km
      p_k = 1000.0 - (real(k, kind=kind_phys)-0.5) * (800.0 / real(km, kind=kind_phys))
      prslp(1, k) = p_k * 100.0
      delp(1, k) = (800.0 / real(km, kind=kind_phys)) * 100.0
      dot(1, k) = -1.0; qlc(1, k) = 0.0; qli(1, k) = 0.0
      do l = 1, 5
        if (p_k <= p_levels(l) .and. p_k > p_levels(l+1)) then
          weight = (p_k - p_levels(l+1)) / (p_levels(l) - p_levels(l+1))
          t1(1, k) = t_levels(l+1) + weight * (t_levels(l) - t_levels(l+1))
          td_k = td_levels(l+1) + weight * (td_levels(l) - td_levels(l+1))
          u1(1, k) = u_levels(l+1) + weight * (u_levels(l) - u_levels(l+1))
          v1(1, k) = v_levels(l+1) + weight * (v_levels(l) - v_levels(l+1))
          goto 20
        endif
      enddo
      t1(1,k) = t_levels(6); td_k = td_levels(6); u1(1,k) = u_levels(6); v1(1,k) = v_levels(6)
20    continue
      if (k > 1) then
        p_km1 = 1000.0 - (real(k-1, kind=kind_phys)-0.5) * (800.0 / real(km, kind=kind_phys))
        phil(1, k) = phil(1, k-1) + rd * 0.5*(t1(1,k)+t1(1,k-1)) * log(p_km1/p_k)
      endif
      es_k = 611.21 * exp(17.67 * (td_k - 273.15) / (td_k - 29.65)) / 100.0
      q1(1, k) = eps * es_k / (p_k - (1.0-eps)*es_k)
    enddo
  end subroutine
end program test_sas_comprehensive
