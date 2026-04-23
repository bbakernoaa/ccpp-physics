program test_multi_tracer
  use machine, only: kind_phys
  use physcons, only: grav => con_g
  use tracer_transport_mod
  implicit none

  integer, parameter :: im = 1, km = 10, ntracers = 3
  real(kind=kind_phys) :: delt = 1200.0_kind_phys
  logical :: cnvflg(im)
  real(kind=kind_phys) :: delp(im, km), mnet(im, km)
  real(kind=kind_phys) :: q_3d(im, km, ntracers), m_init(im, ntracers), m_final(im, ntracers)
  real(kind=kind_phys) :: rn_sink(im, ntracers)
  integer :: n

  cnvflg = .true.
  delp = 10000.0_kind_phys
  mnet = 0.0_kind_phys
  mnet(1, 5) = 0.1_kind_phys
  rn_sink = 0.0_kind_phys

  ! Tracer 1: Uniform
  q_3d(:,:,1) = 1.0e-3_kind_phys
  ! Tracer 2: Hole
  q_3d(:,:,2) = 1.0e-3_kind_phys
  q_3d(1, 5, 2) = -1.0e-4_kind_phys
  ! Tracer 3: Randomish
  do n=1,km
    q_3d(1,n,3) = real(n, kind=kind_phys) * 1.0e-4_kind_phys
  enddo

  write(*,*) "### Testing Multi-Tracer Overloading..."

  call calc_column_mass(im, km, ntracers, cnvflg, q_3d, delp, m_init)
  call transport_tracer_flux_form(im, km, ntracers, delt, cnvflg, delp, mnet, q_3d)
  call calc_column_mass(im, km, ntracers, cnvflg, q_3d, delp, m_final)

  write(*,*) "| Tracer | Initial Mass | Final Mass | Delta |"
  write(*,*) "|---:|---:|---:|---:|"
  do n = 1, ntracers
    write(*, '("| ", I2, " | ", F15.10, " | ", F15.10, " | ", E15.5, " |")') &
     & n, m_init(1,n), m_final(1,n), m_final(1,n) - m_init(1,n)
  enddo

end program test_multi_tracer
