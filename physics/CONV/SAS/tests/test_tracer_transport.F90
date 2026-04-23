program test_tracer_transport
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  implicit none

  integer, parameter :: im = 1, km = 10
  real(kind=kind_phys) :: delt = 1200.0_kind_phys
  logical :: cnvflg(im)
  real(kind=kind_phys) :: delp(im, km), mnet(im, km), q_tracer(im, km)
  real(kind=kind_phys) :: col_mass_init(im), col_mass_final(im)
  real(kind=kind_phys) :: rn_sink(im)
  integer :: i, k

  cnvflg = .true.
  delp = 10000.0_kind_phys ! 100 hPa per layer
  rn_sink = 0.0_kind_phys

  ! Initial tracer concentration
  q_tracer = 1.0e-3_kind_phys

  ! Define a net mass flux that would cause some transport
  mnet = 0.0_kind_phys
  mnet(1, 5) = 0.1_kind_phys ! Significant flux at interface 5 (between layer 5 and 6)

  write(*,*) "Testing tracer_transport_mod routines..."

  ! --- 1. Advective (Original leaky method) ---
  q_tracer = 1.0e-3_kind_phys
  ! Force a negative value that will be clipped
  q_tracer(1, 5) = 0.5e-10_kind_phys
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_init)
  call transport_tracer_advective(im, km, delt, cnvflg, delp, mnet, q_tracer)
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_final)
  write(*, '(A, F20.10, A, E15.5)') "Advective:    Init Mass =", col_mass_init(1), " Delta =", col_mass_final(1) - col_mass_init(1)

  ! --- 2. Flux-Form ---
  q_tracer = 1.0e-3_kind_phys
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_init)
  call transport_tracer_flux_form(im, km, delt, cnvflg, delp, mnet, q_tracer)
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_final)
  write(*, '(A, F20.10, A, E15.5)') "Flux-Form:    Init Mass =", col_mass_init(1), " Delta =", col_mass_final(1) - col_mass_init(1)

  ! --- 3. Hole-Filling ---
  q_tracer = 1.0e-3_kind_phys
  ! Force a negative value to trigger hole filling
  q_tracer(1, 5) = -1.0e-4_kind_phys
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_init)
  call transport_tracer_hole_filling(im, km, delt, cnvflg, delp, mnet, q_tracer)
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_final)
  write(*, '(A, F20.10, A, E15.5)') "Hole-Filling: Init Mass =", col_mass_init(1), " Delta =", col_mass_final(1) - col_mass_init(1)

  ! --- 4. Mass Fixer ---
  q_tracer = 1.0e-3_kind_phys
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_init)
  ! simulate a leak
  q_tracer(1, :) = q_tracer(1, :) * 0.9_kind_phys
  ! Fix it
  call apply_mass_fixer(im, km, cnvflg, delp, rn_sink, col_mass_init, q_tracer)
  call calc_column_mass(im, km, cnvflg, q_tracer, delp, col_mass_final)
  write(*, '(A, F20.10, A, E15.5)') "Mass-Fixer:   Init Mass =", col_mass_init(1), " Delta =", col_mass_final(1) - col_mass_init(1)

end program test_tracer_transport
