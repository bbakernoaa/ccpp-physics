program test_sas_comprehensive
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
  integer :: opt, fixer_loop
  logical :: use_fixer
  character(len=30) :: test_name

  cnvflg = .true.
  delp = 5000.0_kind_phys
  rn_sink = 0.0_kind_phys
  mnet = -2.0_kind_phys

  q_init = 0.02_kind_phys
  q_init(1, 6) = 0.016_kind_phys
  q_init(1, 7:km) = 0.0005_kind_phys

  write(*,*) "| Test Name | Initial Mass | Final Mass | Delta | Status |"
  write(*,*) "|:---|---:|---:|---:|:---|"

  do opt = 3, 1, -1
    do fixer_loop = 0, 1
      use_fixer = (fixer_loop == 1)
      q_tracer = q_init
      call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_init)

      select case(opt)
      case(1)
        call transport_tracer_flux_form(im, km, delt, cnvflg, delp, mnet, q_tracer)
        if (use_fixer) then
           test_name = "Flux with fixer"
        else
           test_name = "Flux without fixer"
        endif
      case(2)
        call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_tracer)
        if (use_fixer) then
           test_name = "Hole filler with fixer"
        else
           test_name = "Hole filler without fixer"
        endif
      case(3)
        call transport_tracer_advective(im, km, delt, cnvflg, delp, mnet, q_tracer)
        if (use_fixer) then
           test_name = "Legacy (clipping) + fixer"
        else
           test_name = "Legacy"
        endif
      end select

      if (use_fixer) then
        call apply_mass_fixer(im, km, cnvflg, delp, rn_sink, m_init, q_tracer)
      endif

      call calc_column_mass(im, km, cnvflg, q_tracer, delp, m_final)

      write(*, '("| ", A28, " | ", F15.10, " | ", F15.10, " | ", E15.5, " | ", A, " |")') &
       & test_name, m_init(1), m_final(1), m_final(1) - m_init(1), &
       & merge("Conservative    ", "Leaky (Clipped) ", m_final(1)-m_init(1) < 1.0e-10_kind_phys)
    enddo
  enddo

end program test_sas_comprehensive
