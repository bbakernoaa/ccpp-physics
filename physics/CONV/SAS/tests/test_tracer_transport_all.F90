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
  real(kind=kind_phys) :: m_init(im), m_final(im)
  real(kind=kind_phys) :: rn_sink(im)
  integer :: opt

  cnvflg = .true.
  delp = 10000.0_kind_phys
  rn_sink = 0.0_kind_phys
  mnet = 0.0_kind_phys
  mnet(1, 5) = 0.1_kind_phys

  ! Consistent initial distribution with a "hole" (negative value)
  ! to show how methods differ.
  q_init = 1.0e-3_kind_phys
  q_init(1, 5) = -1.0e-4_kind_phys

  write(*,*) "| Transport Option | Initial Mass | Final Mass | Delta | Status |"
  write(*,*) "|:---|---:|---:|---:|:---|"

  do opt = 1, 3
    q_tracer = q_init
    call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_init)

    select case(opt)
    case(1)
      call transport_tracer_flux_form(im, km, delt, cnvflg, delp, mnet, q_tracer)
      call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_final)
      write(*, '("| Flux-Form | ", F15.10, " | ", F15.10, " | ", E15.5, " | Conservative |")') &
     &  m_init(1), m_final(1), m_final(1) - m_init(1)
    case(2)
      call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_tracer)
      call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_final)
      write(*, '("| Hole-Filling | ", F15.10, " | ", F15.10, " | ", E15.5, " | Conservative |")') &
     &  m_init(1), m_final(1), m_final(1) - m_init(1)
    case(3)
      call transport_tracer_advective(im, km, delt, cnvflg, delp, mnet, q_tracer)
      call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_final)
      write(*, '("| Advective | ", F15.10, " | ", F15.10, " | ", E15.5, " | Leaky (Clipped) |")') &
     &  m_init(1), m_final(1), m_final(1) - m_init(1)
    end select
  enddo

  ! Mass Fixer test (using the same initial mass but simulated leak)
  q_tracer = q_init
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_init)
  q_tracer = q_tracer * 0.9_kind_phys ! Simulated 10% leak
  call apply_mass_fixer(im, km, cnvflg, delp, rn_sink, m_init, q_tracer)
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_final)
  write(*, '("| Mass-Fixer | ", F15.10, " | ", F15.10, " | ", E15.5, " | Conservative |")') &
     &  m_init(1), m_final(1), m_final(1) - m_init(1)

end program test_tracer_transport_all
