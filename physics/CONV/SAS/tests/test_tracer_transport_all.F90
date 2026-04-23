program test_tracer_transport_all
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  implicit none

  integer, parameter :: im = 1, km = 10
  real(kind=kind_phys) :: delt = 600.0_kind_phys
  logical :: cnvflg(im)
  real(kind=kind_phys) :: delp(im, km), mnet(im, km)
  real(kind=kind_phys) :: q_tracer(im, km), q_init(im, km)
  real(kind=kind_phys) :: m_init(im), m_final(im)
  real(kind=kind_phys) :: rn_sink(im)

  cnvflg = .true.
  rn_sink = 0.0_kind_phys

  write(*,*) "### 1. Scenario: Fillable Hole (Local Redistribution)"
  delp = 10000.0_kind_phys
  mnet = 0.0_kind_phys
  mnet(1, 5) = 0.1_kind_phys
  q_init = 1.0e-3_kind_phys
  q_init(1, 5) = -1.0e-4_kind_phys
  call run_comparison(q_init, im, km, delt, cnvflg, delp, mnet, rn_sink)

  write(*,*)
  write(*,*) "### 2. Scenario: The 'Dry Cliff' (Overshoot & Clipping)"
  ! Subsidence pushing dry air into a thin moist layer
  delp = 5000.0_kind_phys
  mnet = -2.0_kind_phys ! Massive subsidence
  q_init = 0.02_kind_phys ! Moist surface
  q_init(1, 6) = 0.016_kind_phys ! Target layer
  q_init(1, 7:km) = 0.0005_kind_phys ! Dry cliff above
  call run_comparison(q_init, im, km, delt, cnvflg, delp, mnet, rn_sink)

  write(*,*)
  write(*,*) "### 3. Scenario: Massive Hole (Total column mass negative)"
  delp = 10000.0_kind_phys
  mnet = 0.0_kind_phys
  q_init = 1.0e-6_kind_phys
  q_init(1, 5) = -1.0e-2_kind_phys
  call run_comparison(q_init, im, km, delt, cnvflg, delp, mnet, rn_sink)

contains

  subroutine run_comparison(q_base, im, km, delt, cnvflg, delp, mnet, rn_sink)
    real(kind=kind_phys), intent(in) :: q_base(im, km)
    integer, intent(in) :: im, km
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(im)
    real(kind=kind_phys), intent(in) :: delp(im, km), mnet(im, km), rn_sink(im)
    real(kind=kind_phys) :: q_work(im, km), m_i(im), m_f(im)

    write(*,*) "| Method | Initial Mass | Final Mass | Delta | Min(q) | Status |"
    write(*,*) "|:---|---:|---:|---:|---:|:---|"

    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_flux_form(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Flux-Form | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " | Conservative |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Hole-Filling | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " | Conservative |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_advective(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Advective | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " | Leaky (Clipped) |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

  end subroutine run_comparison

end program test_tracer_transport_all
