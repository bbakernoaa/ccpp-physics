program test_fixer_conservation
  implicit none
  integer, parameter :: kind_phys = 8
  integer, parameter :: km = 30
  real(kind=kind_phys) :: ctr(km), ctr_orig(km), delp(km), grav, tsumn, tsump, rtnp, tem
  real(kind=kind_phys) :: mass_initial, mass_after_transport, mass_after_fixer, mass_after_clipping
  integer :: k

  grav = 9.80665_kind_phys
  delp(:) = 3000.0_kind_phys

  ! 1. Initial State (all positive)
  ctr(:) = 1.0e-6_kind_phys
  mass_initial = sum(ctr * delp / grav)

  ! 2. After Transport (one layer becomes negative)
  ! Imagine transport created a negative at the top layer km
  ctr(km) = -1.0e-7_kind_phys
  mass_after_transport = sum(ctr * delp / grav)
  ctr_orig = ctr

  print *, '--- BASELINE (km-1) ---'
  ! 3. Mass Fixer (Baseline: km-1)
  tsumn = 0.0
  tsump = 0.0
  do k = 1, km-1
    tem = ctr(k) * delp(k) / grav
    if (ctr(k) < 0.0) tsumn = tsumn + tem
    if (ctr(k) > 0.0) tsump = tsump + tem
  enddo
  rtnp = 1.0
  if (tsump > 0.0 .and. tsumn < 0.0) then
    if (tsump > abs(tsumn)) then
      rtnp = tsumn / tsump
    else
      rtnp = tsump / tsumn
    endif
  endif
  do k = 1, km-1
    if (rtnp < 0.0) then
      if (tsump > abs(tsumn)) then
        if (ctr(k) < 0.0) ctr(k) = 0.0
        if (ctr(k) > 0.0) ctr(k) = (1.0 + rtnp) * ctr(k)
      else
        if (ctr(k) < 0.0) ctr(k) = (1.0 + rtnp) * ctr(k)
        if (ctr(k) > 0.0) ctr(k) = 0.0
      endif
    endif
  enddo
  mass_after_fixer = sum(ctr * delp / grav)

  ! 4. Clipping (what the host model or subsequent physics would do)
  do k = 1, km
    if (ctr(k) < 0.0) ctr(k) = 0.0
  enddo
  mass_after_clipping = sum(ctr * delp / grav)

  print *, 'Initial Mass:        ', mass_initial
  print *, 'Mass after Transport:', mass_after_transport
  print *, 'Mass after Fixer:    ', mass_after_fixer
  print *, 'Mass after Clipping: ', mass_after_clipping
  print *, 'Conservation Error:  ', mass_after_clipping - mass_after_transport
  print *, 'Tracer at KM:        ', ctr(km)

  ! --- FIXED ---
  print *, ''
  print *, '--- FIXED (km) ---'
  ctr = ctr_orig
  tsumn = 0.0
  tsump = 0.0
  do k = 1, km
    tem = ctr(k) * delp(k) / grav
    if (ctr(k) < 0.0) tsumn = tsumn + tem
    if (ctr(k) > 0.0) tsump = tsump + tem
  enddo
  rtnp = 1.0
  if (tsump > 0.0 .and. tsumn < 0.0) then
    if (tsump > abs(tsumn)) then
      rtnp = tsumn / tsump
    else
      rtnp = tsump / tsumn
    endif
  endif
  do k = 1, km
    if (rtnp < 0.0) then
      if (tsump > abs(tsumn)) then
        if (ctr(k) < 0.0) ctr(k) = 0.0
        if (ctr(k) > 0.0) ctr(k) = (1.0 + rtnp) * ctr(k)
      else
        if (ctr(k) < 0.0) ctr(k) = (1.0 + rtnp) * ctr(k)
        if (ctr(k) > 0.0) ctr(k) = 0.0
      endif
    endif
  enddo
  mass_after_fixer = sum(ctr * delp / grav)

  ! Clipping (should do nothing now)
  do k = 1, km
    if (ctr(k) < 0.0) ctr(k) = 0.0
  enddo
  mass_after_clipping = sum(ctr * delp / grav)

  print *, 'Initial Mass:        ', mass_initial
  print *, 'Mass after Transport:', mass_after_transport
  print *, 'Mass after Fixer:    ', mass_after_fixer
  print *, 'Mass after Clipping: ', mass_after_clipping
  print *, 'Conservation Error:  ', mass_after_clipping - mass_after_transport
  print *, 'Tracer at KM:        ', ctr(km)

end program
