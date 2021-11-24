!======================================================================
!  Nudged Elastic Band algorithm
!======================================================================
! f2py compilation example:
! python3 -m numpy.f2py -c neb_fortran.f90 -m neb_fortran --opt='-Ofast'
!
!  Functions
!  ---------
!
!    neb
!

subroutine neb(mep_guess, grid, Z, dZ, max_step, max_iter, grid_d, &
               spring, span, rescale, direction, n, nknots, mep_path)

  !--------------------------------------------------------------------
  ! Nudged Elastic Band over grid points w/ LOWESS interpolation
  !--------------------------------------------------------------------

  implicit none

  ! Variable definition
  integer, intent(in)  :: n
  integer, intent(in)  :: nknots

  real(8), intent(in)  :: mep_guess(nknots,2)
  real(8), intent(in)  :: grid(n,2), Z(n), dZ(n,2)
  real(8), intent(in)  :: max_step
  integer, intent(in)  :: max_iter
  real(8), intent(in)  :: grid_d
  real(8), intent(in)  :: spring
  real(8), intent(in)  :: span
  real(8), intent(in)  :: rescale(2)                   ! re-scale forces limits (min,max)
  integer, intent(in)  :: direction                    ! '1' downhill or '-1' uphill

  real(8), intent(out) :: mep_path(nknots,2)

  integer              :: i, iter, maxpoint(1)
  real(8)              :: mep_last(nknots,2)            ! previous iteration mep_path
  real(8)              :: E_path(nknots)                ! energy of every point
  real(8)              :: T_path(nknots,2)              ! normal tangent-to-the-path vector
  real(8)              :: P_path(nknots,2)              ! perpendicular-to-the-path force vector
  real(8)              :: S_path(nknots,2)              ! NEB spring force vector
  real(8)              :: P_mod(nknots), S_mod(nknots)  ! modulus of every point
  real(8)              :: F_path(nknots,2)              ! final force vector
  real(8)              :: gradx, grady, grad(2)
  real(8)              :: step, displacement
  real(8)              :: E_min, E_max, E_diff
  logical              :: p_in

  ! initialize variables
  mep_path = mep_guess
  mep_last = mep_path
  E_path   = 0.D0
  T_path   = 0.D0
  P_path   = 0.D0
  S_path   = 0.D0
  F_path   = 0.D0
  P_mod    = 0.D0
  S_mod    = 0.D0
  step     = 1.D9

  ! main NEB loop
  do iter=1,max_iter

    ! energy of the MEP points
    do i=1,nknots
      CALL grid2point_lowess(mep_path(i,:), grid, Z, span, n, E_path(i))
    enddo

    do i=2,nknots-1

      ! tangent vector by bisection
      T_path(i,:) = vnorm2( vnorm2( mep_path(i,:) - mep_path(i-1,:) ) &
                          + vnorm2( mep_path(i+1,:) - mep_path(i,:) ) )

      ! gradient
      CALL grid2point_lowess(mep_path(i,:), grid, dZ(:,1), span, n, gradx)
      CALL grid2point_lowess(mep_path(i,:), grid, dZ(:,2), span, n, grady)
      grad = - direction * (/ gradx, grady /)

      ! perpendicular force
      P_path(i,:) = grad - vproj2(grad,T_path(i,:))

      ! spring force
      S_path(i,:) = spring * ( (mep_path(i+1,:) - mep_path(i,:)) + (mep_path(i-1,:) - mep_path(i,:)) )

      ! modulus
      P_mod(i) = NORM2(P_path(i,:))
      S_mod(i) = NORM2(S_path(i,:))

    enddo

    ! rescale forces respect the energy difference
    E_min = MINVAL(E_path(2:nknots-1))
    E_max = MAXVAL(E_path(2:nknots-1))
    E_diff = E_max - E_min
    do i=2,nknots-1
      if (ABS(P_mod(i)) < 1.D-05) then
        P_mod(i) = 0.D0
      else
        P_path(i,:) = P_path(i,:) / P_mod(i) * ( (E_path(i) - E_min)/E_diff * (rescale(2)-rescale(1)) + rescale(1) )
      endif
      if (ABS(S_mod(i)) < 1.D-05) then
        S_mod(i) = 0.D0
      else
        S_path(i,:) = S_path(i,:) / S_mod(i)
      endif
    enddo

    ! adapt step size
    maxpoint = MAXLOC(E_path(2:nknots-1))
    step = MIN(max_step, ABS(P_mod(maxpoint(1))/E_max))

    ! final force vector
    F_path = P_path/2.D0 + S_path/2.D0

    ! move MEP points
    mep_path = mep_path + F_path * step

    ! check if any point is out of the grid
    do i=2,nknots-1
      CALL point_in(mep_path(i,:), grid, grid_d/SQRT(2.D0), 1, n, p_in)
      if ( .NOT. p_in) then
        mep_path(i,:) = mep_last(i,:)
      endif
    enddo

    ! calculate convergence
    displacement = 0.D0
    do i=2,nknots-1
      displacement = displacement + NORM2(mep_path(i,:) - mep_last(i,:))
    enddo
    displacement = displacement / nknots

    ! save mep path to compare
    mep_last = mep_path

    ! PRINT*, iter, displacement, P_mod(maxpoint(1)), step

    ! check convergence
    if ( displacement < 1.D-5 .OR. P_mod(maxpoint(1)) < 1.D-3 ) EXIT

  enddo


  contains

    ! normalize vector - 2D -------------------------------------------
    function vnorm2(vec)

      implicit none
      real(8), intent(in)  :: vec(2)
      real(8)              :: vnorm2(2)

      vnorm2 = vec / NORM2(vec)

    end function

    ! vector projection - 2D ------------------------------------------
    function vproj2(vec1, vec2)

      implicit none
      real(8), intent(in)  :: vec1(2), vec2(2)
      real(8)              :: vproj2(2)

      if (ALL(vec2==0.D0)) then
        vproj2 = 0.D0
      else
        vproj2 = vec2 * DOT_PRODUCT(vec1,vec2) / DOT_PRODUCT(vec2,vec2)
      endif

    end function

end subroutine

include '../interpolation_fortran.f90'