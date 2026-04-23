! =========================================================================
! MODULE: Tracer Transport and Mass Conservation for Deep/Shallow Convection
! =========================================================================
module tracer_transport_mod
  use machine,  only: kind_phys
  use physcons, only: grav => con_g
  implicit none
  private

  public :: calc_column_mass, transport_tracer_flux_form, &
            transport_tracer_advective, transport_tracer_hole_filling, &
            apply_mass_fixer

  interface calc_column_mass
    module procedure calc_column_mass_2d
    module procedure calc_column_mass_3d
  end interface

  interface transport_tracer_flux_form
    module procedure transport_tracer_flux_form_2d
    module procedure transport_tracer_flux_form_3d
  end interface

  interface transport_tracer_hole_filling
    module procedure transport_tracer_hole_filling_2d
    module procedure transport_tracer_hole_filling_3d
  end interface

  interface transport_tracer_advective
    module procedure transport_tracer_advective_2d
    module procedure transport_tracer_advective_3d
  end interface

  interface apply_mass_fixer
    module procedure apply_mass_fixer_2d
    module procedure apply_mass_fixer_3d
  end interface

  real(kind=kind_phys), parameter :: rhowater = 1000.0_kind_phys ! kg/m^3

contains

  ! =======================================================================
  ! ROUTINE 1: Calculate Column-Integrated Mass
  ! =======================================================================
  subroutine calc_column_mass_2d(im, km, cnvflg, q_tracer, delp, col_mass)
    integer, intent(in) :: im, km
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: q_tracer(:,:), delp(:,:)
    real(kind=kind_phys), intent(out) :: col_mass(:)
    integer :: i, k

    do i = 1, im
      col_mass(i) = 0.0_kind_phys
    enddo

    do k = 1, km
      do i = 1, im
        if (cnvflg(i)) then
          col_mass(i) = col_mass(i) + (q_tracer(i,k) * delp(i,k) / grav)
        endif
      enddo
    enddo
  end subroutine calc_column_mass_2d

  subroutine calc_column_mass_3d(im, km, ntracers, cnvflg, q_tracer, delp, col_mass)
    integer, intent(in) :: im, km, ntracers
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: q_tracer(:,:,:), delp(:,:)
    real(kind=kind_phys), intent(out) :: col_mass(:,:) ! (im, ntracers)
    integer :: i, k, n

    do n = 1, ntracers
      do i = 1, im
        col_mass(i,n) = 0.0_kind_phys
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            col_mass(i,n) = col_mass(i,n) + (q_tracer(i,k,n) * delp(i,k) / grav)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_column_mass_3d

  ! =======================================================================
  ! ROUTINE 2: Strict Flux-Form Transport (Perfect Conservation)
  ! =======================================================================
  subroutine transport_tracer_flux_form_2d(im, km, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:)

    integer :: i, k
    real(kind=kind_phys) :: layer_mass_rate, layer_mass_above_rate
    real(kind=kind_phys) :: q_upwind, max_flux_up, max_flux_dn
    real(kind=kind_phys), dimension(im, 0:km) :: q_flux

    do i = 1, im
      q_flux(i, 0)  = 0.0_kind_phys
      q_flux(i, km) = 0.0_kind_phys
    enddo

    do k = 1, km - 1
      do i = 1, im
        if (cnvflg(i)) then
          layer_mass_rate       = delp(i,k)   / (grav * delt)
          layer_mass_above_rate = delp(i,k+1) / (grav * delt)

          q_upwind = merge(q_tracer(i,k), q_tracer(i,k+1), mnet(i,k) > 0.0_kind_phys)
          max_flux_up = q_tracer(i,k) * layer_mass_rate
          max_flux_dn = -q_tracer(i,k+1) * layer_mass_above_rate

          q_flux(i,k) = merge( min(mnet(i,k) * q_upwind, max_flux_up), &
                               max(mnet(i,k) * q_upwind, max_flux_dn), &
                               mnet(i,k) > 0.0_kind_phys )
        else
          q_flux(i,k) = 0.0_kind_phys
        endif
      enddo
    enddo

    do k = 1, km
      do i = 1, im
        if (cnvflg(i)) then
          q_tracer(i,k) = q_tracer(i,k) + (q_flux(i,k-1) - q_flux(i,k)) * delt * grav / delp(i,k)
        endif
      enddo
    enddo
  end subroutine transport_tracer_flux_form_2d

  subroutine transport_tracer_flux_form_3d(im, km, ntracers, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km, ntracers
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:,:)

    integer :: i, k, n
    real(kind=kind_phys) :: layer_mass_rate, layer_mass_above_rate
    real(kind=kind_phys) :: q_upwind, max_flux_up, max_flux_dn
    real(kind=kind_phys), dimension(im, 0:km, ntracers) :: q_flux

    do n = 1, ntracers
      do i = 1, im
        q_flux(i, 0, n)  = 0.0_kind_phys
        q_flux(i, km, n) = 0.0_kind_phys
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km - 1
        do i = 1, im
          if (cnvflg(i)) then
            layer_mass_rate       = delp(i,k)   / (grav * delt)
            layer_mass_above_rate = delp(i,k+1) / (grav * delt)

            q_upwind = merge(q_tracer(i,k,n), q_tracer(i,k+1,n), mnet(i,k) > 0.0_kind_phys)
            max_flux_up = q_tracer(i,k,n) * layer_mass_rate
            max_flux_dn = -q_tracer(i,k+1,n) * layer_mass_above_rate

            q_flux(i,k,n) = merge( min(mnet(i,k) * q_upwind, max_flux_up), &
                                   max(mnet(i,k) * q_upwind, max_flux_dn), &
                                   mnet(i,k) > 0.0_kind_phys )
          else
            q_flux(i,k,n) = 0.0_kind_phys
          endif
        enddo
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            q_tracer(i,k,n) = q_tracer(i,k,n) + (q_flux(i,k-1,n) - q_flux(i,k,n)) * delt * grav / delp(i,k)
          endif
        enddo
      enddo
    enddo
  end subroutine transport_tracer_flux_form_3d

  ! =======================================================================
  ! ROUTINE 3: Conservative Advection + Two-Pass Hole-Filling
  ! =======================================================================
  subroutine transport_tracer_hole_filling_2d(im, km, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:)

    integer :: i, k
    real(kind=kind_phys) :: mass_deficit, q_upwind
    real(kind=kind_phys), dimension(im, 0:km) :: q_flux

    do i = 1, im
      q_flux(i, 0)  = 0.0_kind_phys
      q_flux(i, km) = 0.0_kind_phys
    enddo

    do k = 1, km - 1
      do i = 1, im
        if (cnvflg(i)) then
          q_upwind = merge(q_tracer(i,k), q_tracer(i,k+1), mnet(i,k) > 0.0_kind_phys)
          q_flux(i,k) = mnet(i,k) * q_upwind
        else
          q_flux(i,k) = 0.0_kind_phys
        endif
      enddo
    enddo

    do k = 1, km
      do i = 1, im
        if (cnvflg(i)) then
          q_tracer(i,k) = q_tracer(i,k) + (q_flux(i,k-1) - q_flux(i,k)) * delt * grav / delp(i,k)
        endif
      enddo
    enddo

    do k = km, 2, -1
      do i = 1, im
        if (cnvflg(i) .and. q_tracer(i,k) < 0.0_kind_phys) then
          mass_deficit = abs(q_tracer(i,k)) * (delp(i,k) / grav)
          q_tracer(i,k-1) = q_tracer(i,k-1) - (mass_deficit / (delp(i,k-1) / grav))
          q_tracer(i,k) = 0.0_kind_phys
        endif
      enddo
    enddo

    do k = 1, km - 1
      do i = 1, im
        if (cnvflg(i) .and. q_tracer(i,k) < 0.0_kind_phys) then
          mass_deficit = abs(q_tracer(i,k)) * (delp(i,k) / grav)
          q_tracer(i,k+1) = q_tracer(i,k+1) - (mass_deficit / (delp(i,k+1) / grav))
          q_tracer(i,k) = 0.0_kind_phys
        endif
      enddo
    enddo

    do i = 1, im
      if (cnvflg(i) .and. q_tracer(i,km) < 0.0_kind_phys) then
        q_tracer(i,km) = 0.0_kind_phys
      endif
    enddo
  end subroutine transport_tracer_hole_filling_2d

  subroutine transport_tracer_hole_filling_3d(im, km, ntracers, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km, ntracers
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:,:)

    integer :: i, k, n
    real(kind=kind_phys) :: mass_deficit, q_upwind
    real(kind=kind_phys), dimension(im, 0:km, ntracers) :: q_flux

    do n = 1, ntracers
      do i = 1, im
        q_flux(i, 0, n)  = 0.0_kind_phys
        q_flux(i, km, n) = 0.0_kind_phys
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km - 1
        do i = 1, im
          if (cnvflg(i)) then
            q_upwind = merge(q_tracer(i,k,n), q_tracer(i,k+1,n), mnet(i,k) > 0.0_kind_phys)
            q_flux(i,k,n) = mnet(i,k) * q_upwind
          else
            q_flux(i,k,n) = 0.0_kind_phys
          endif
        enddo
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            q_tracer(i,k,n) = q_tracer(i,k,n) + (q_flux(i,k-1,n) - q_flux(i,k,n)) * delt * grav / delp(i,k)
          endif
        enddo
      enddo
    enddo

    do n = 1, ntracers
      do k = km, 2, -1
        do i = 1, im
          if (cnvflg(i) .and. q_tracer(i,k,n) < 0.0_kind_phys) then
            mass_deficit = abs(q_tracer(i,k,n)) * (delp(i,k) / grav)
            q_tracer(i,k-1,n) = q_tracer(i,k-1,n) - (mass_deficit / (delp(i,k-1) / grav))
            q_tracer(i,k,n) = 0.0_kind_phys
          endif
        enddo
      enddo

      do k = 1, km - 1
        do i = 1, im
          if (cnvflg(i) .and. q_tracer(i,k,n) < 0.0_kind_phys) then
            mass_deficit = abs(q_tracer(i,k,n)) * (delp(i,k) / grav)
            q_tracer(i,k+1,n) = q_tracer(i,k+1,n) - (mass_deficit / (delp(i,k+1) / grav))
            q_tracer(i,k,n) = 0.0_kind_phys
          endif
        enddo
      enddo

      do i = 1, im
        if (cnvflg(i) .and. q_tracer(i,km,n) < 0.0_kind_phys) then
          q_tracer(i,km,n) = 0.0_kind_phys
        endif
      enddo
    enddo
  end subroutine transport_tracer_hole_filling_3d

  ! =======================================================================
  ! ROUTINE 4: Traditional Advection (Original Leaky Method for Baseline)
  ! =======================================================================
  subroutine transport_tracer_advective_2d(im, km, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:)

    integer :: i, k
    real(kind=kind_phys) :: q_tendency
    real(kind=kind_phys), parameter :: qmin = 1.0e-10_kind_phys

    do k = 1, km
      do i = 1, im
        if (cnvflg(i)) then
          if (k == 1) then
             q_tendency = merge(mnet(i,k) * (q_tracer(i,k+1) - q_tracer(i,k)) / delp(i,k), &
                                mnet(i,k) * q_tracer(i,k) / delp(i,k), &
                                mnet(i,k) <= 0.0_kind_phys)
          else if (k == km) then
             q_tendency = merge(-mnet(i,k-1) * q_tracer(i,k) / delp(i,k), &
                                mnet(i,k-1) * (q_tracer(i,k) - q_tracer(i,k-1)) / delp(i,k), &
                                mnet(i,k-1) <= 0.0_kind_phys)
          else
             if (mnet(i,k) <= 0.0_kind_phys) then
                q_tendency = mnet(i,k) * (q_tracer(i,k+1) - q_tracer(i,k)) / delp(i,k)
             else
                q_tendency = mnet(i,k-1) * (q_tracer(i,k) - q_tracer(i,k-1)) / delp(i,k)
             endif
          endif
          q_tracer(i,k) = q_tracer(i,k) + (q_tendency * grav * delt)
          q_tracer(i,k) = max(q_tracer(i,k), qmin)
        endif
      enddo
    enddo
  end subroutine transport_tracer_advective_2d

  subroutine transport_tracer_advective_3d(im, km, ntracers, delt, cnvflg, delp, mnet, q_tracer)
    integer, intent(in) :: im, km, ntracers
    real(kind=kind_phys), intent(in) :: delt
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), mnet(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:,:)

    integer :: i, k, n
    real(kind=kind_phys) :: q_tendency
    real(kind=kind_phys), parameter :: qmin = 1.0e-10_kind_phys

    do n = 1, ntracers
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k == 1) then
               q_tendency = merge(mnet(i,k) * (q_tracer(i,k+1,n) - q_tracer(i,k,n)) / delp(i,k), &
                                  mnet(i,k) * q_tracer(i,k,n) / delp(i,k), &
                                  mnet(i,k) <= 0.0_kind_phys)
            else if (k == km) then
               q_tendency = merge(-mnet(i,k-1) * q_tracer(i,k,n) / delp(i,k), &
                                  mnet(i,k-1) * (q_tracer(i,k,n) - q_tracer(i,k-1,n)) / delp(i,k), &
                                  mnet(i,k-1) <= 0.0_kind_phys)
            else
               if (mnet(i,k) <= 0.0_kind_phys) then
                  q_tendency = mnet(i,k) * (q_tracer(i,k+1,n) - q_tracer(i,k,n)) / delp(i,k)
               else
                  q_tendency = mnet(i,k-1) * (q_tracer(i,k,n) - q_tracer(i,k-1,n)) / delp(i,k)
               endif
            endif
            q_tracer(i,k,n) = q_tracer(i,k,n) + (q_tendency * grav * delt)
            q_tracer(i,k,n) = max(q_tracer(i,k,n), qmin)
          endif
        enddo
      enddo
    enddo
  end subroutine transport_tracer_advective_3d

  ! =======================================================================
  ! ROUTINE 5: Proportional Mass Fixer
  ! =======================================================================
  subroutine apply_mass_fixer_2d(im, km, cnvflg, delp, rn_sink, col_mass_init, q_tracer)
    integer, intent(in) :: im, km
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), rn_sink(:), col_mass_init(:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:)

    integer :: i, k
    real(kind=kind_phys), dimension(im) :: col_mass_final, mass_error, factor

    call calc_column_mass_2d(im, km, cnvflg, q_tracer, delp, col_mass_final)

    do i = 1, im
      if (cnvflg(i)) then
        mass_error(i) = col_mass_init(i) - col_mass_final(i) - (rn_sink(i) * rhowater)
        factor(i) = merge(1.0_kind_phys + (mass_error(i) / max(col_mass_final(i), 1.0e-12_kind_phys)), &
                          1.0_kind_phys, &
                          abs(mass_error(i)) > 1.0e-12_kind_phys .and. col_mass_final(i) > 0.0_kind_phys)
      else
        factor(i) = 1.0_kind_phys
      endif
    enddo

    do k = 1, km
      do i = 1, im
        if (cnvflg(i)) then
          q_tracer(i,k) = q_tracer(i,k) * factor(i)
        endif
      enddo
    enddo
  end subroutine apply_mass_fixer_2d

  subroutine apply_mass_fixer_3d(im, km, ntracers, cnvflg, delp, rn_sink, col_mass_init, q_tracer)
    integer, intent(in) :: im, km, ntracers
    logical, intent(in) :: cnvflg(:)
    real(kind=kind_phys), intent(in) :: delp(:,:), rn_sink(:,:), col_mass_init(:,:)
    real(kind=kind_phys), intent(inout) :: q_tracer(:,:,:)

    integer :: i, k, n
    real(kind=kind_phys), dimension(im, ntracers) :: col_mass_final, mass_error, factor

    call calc_column_mass_3d(im, km, ntracers, cnvflg, q_tracer, delp, col_mass_final)

    do n = 1, ntracers
      do i = 1, im
        if (cnvflg(i)) then
          mass_error(i,n) = col_mass_init(i,n) - col_mass_final(i,n) - (rn_sink(i,n) * rhowater)
          factor(i,n) = merge(1.0_kind_phys + (mass_error(i,n) / max(col_mass_final(i,n), 1.0e-12_kind_phys)), &
                            1.0_kind_phys, &
                            abs(mass_error(i,n)) > 1.0e-12_kind_phys .and. col_mass_final(i,n) > 0.0_kind_phys)
        else
          factor(i,n) = 1.0_kind_phys
        endif
      enddo
    enddo

    do n = 1, ntracers
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            q_tracer(i,k,n) = q_tracer(i,k,n) * factor(i,n)
          endif
        enddo
      enddo
    enddo
  end subroutine apply_mass_fixer_3d

end module tracer_transport_mod
