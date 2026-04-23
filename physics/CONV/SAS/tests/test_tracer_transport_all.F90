program test_tracer_transport_all
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  implicit none

  integer, parameter :: im = 1, km = 10
  real(kind=kind_phys) :: delt = 1200.0_kind_phys
  logical :: cnvflg(im)
  real(kind=kind_phys) :: delp(im, km), mnet(im, km)
  real(kind=kind_phys) :: q_tracer(im, km), q_init(im, km)
  real(kind=kind_phys) :: rn_sink(im)

  cnvflg = .true.
  delp = 10000.0_kind_phys
  rn_sink = 0.0_kind_phys
  mnet = 0.0_kind_phys
  mnet(1, 5) = 0.1_kind_phys

  write(*,*) "### Scenario 1: Small Hole (Fillable Locally)"
  q_init = 1.0e-3_kind_phys
  q_init(1, 5) = -1.0e-4_kind_phys
  call run_comparison(q_init, im, km, delt, cnvflg, delp, mnet, rn_sink)

  write(*,*)
  write(*,*) "### Scenario 2: Massive Hole (Total column mass negative)"
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

    write(*,*) "| Method | Initial Mass | Final Mass | Delta | Min(q) |"
    write(*,*) "|:---|---:|---:|---:|---:|"

    ! 1. Flux-Form
    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_flux_form(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Flux-Form | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

    ! 2. Hole-Filling Only
    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Hole-Filling Only | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

    ! 3. Hole-Filling + Mass Fixer
    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_work)
    call apply_mass_fixer(im, km, cnvflg, delp, rn_sink, m_i, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Hole-Fill + Fixer | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

    ! 4. Original Advective
    q_work = q_base
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_i)
    call transport_tracer_advective(im, km, delt, cnvflg, delp, mnet, q_work)
    call calc_column_mass(im, km, cnvflg, q_work, delp, m_f)
    write(*, '("| Advective (Leaky) | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", E12.5, " |")') &
     & m_i(1), m_f(1), m_f(1) - m_i(1), minval(q_work)

  end subroutine run_comparison

end program test_tracer_transport_all
