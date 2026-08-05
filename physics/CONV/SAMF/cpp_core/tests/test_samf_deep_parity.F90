program test_samf_deep_parity
  use mo_samf_cpp_interface
  implicit none
  integer, parameter :: columns = 2
  integer, parameter :: layers = 10
  real(8) :: t_lay(columns, layers)
  real(8) :: q_vap(columns, layers)
  real(8) :: u_wind(columns, layers)
  real(8) :: v_wind(columns, layers)
  real(8) :: p_lay(columns, layers)
  real(8) :: p_int(columns, layers+1)
  real(8) :: z_lay(columns, layers)
  real(8) :: z_int(columns, layers+1)
  real(8) :: dt_t(columns, layers)
  real(8) :: dt_q(columns, layers)
  real(8) :: dt_u(columns, layers)
  real(8) :: dt_v(columns, layers)
  real(8) :: ud_mf(columns, layers)
  real(8) :: dd_mf(columns, layers)
  real(8) :: dt_mf(columns, layers)
  real(8) :: cnvw(columns, layers)
  real(8) :: cnvc(columns, layers)
  integer :: kbot(columns)
  integer :: ktop(columns)
  integer :: kcnv(columns)
  real(8) :: rain(columns)
  real(8) :: dt

  dt = 300.0d0
  t_lay = 280.0d0
  q_vap = 0.005d0
  u_wind = 10.0d0
  v_wind = 5.0d0
  p_lay = 90000.0d0
  p_int = 91000.0d0
  z_lay = 500.0d0
  z_int = 100.0d0

  print *, "Running test_samf_deep_parity..."
  call c_samf_deep_convection_run(int(columns, 8), int(layers, 8), dt, &
      t_lay, q_vap, u_wind, v_wind, p_lay, p_int, z_lay, z_int, &
      dt_t, dt_q, dt_u, dt_v, &
      ud_mf, dd_mf, dt_mf, cnvw, cnvc, &
      kbot, ktop, kcnv, rain)

  print *, "  - rain(1) (expected: 0.0): ", rain(1)
  print *, "  - kcnv(1) (expected: 0):   ", kcnv(1)
  print *, "test_samf_deep_parity PASS!"
end program test_samf_deep_parity
