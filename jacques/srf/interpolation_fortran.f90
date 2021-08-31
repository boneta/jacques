!======================================================================
!  Interpolation Operations
!======================================================================
! f2py compilation example:
! python3 -m numpy.f2py -c interpolation_fortran.f90 -m interpolation_fortran --opt='-Ofast'
!
!  Functions
!  ---------
!
!    point_in
!    smooth_sma
!    grid2point_lowess
!    grid2grid_lowess
!

subroutine point_in(coord, grid, thr, n, m, result)

  !--------------------------------------------------------------------
  ! Evaluate if some coordinates are inside a grid (close to any grid point)
  !--------------------------------------------------------------------

  implicit none

  integer, intent(in)  :: n, m
  real(8), intent(in)  :: coord(n), grid(m,n), thr

  logical, intent(out) :: result

  integer              :: i, j
  real(8)              :: dist(m)

  do i=1,m
    dist(i) = 0.D0
    do j=1,n
      dist(i) = dist(i) + (grid(i,j) - coord(j))**2
    enddo
  enddo

  ! result true if less than thr, else false
  result = SQRT(MINVAL(dist)) <= thr

end subroutine

subroutine smooth_sma(grid, neighbours, n, grid_smooth)

  !--------------------------------------------------------------------
  ! Smoothing of an array: Simple Moving central-Average (SMA)
  !--------------------------------------------------------------------

  implicit none

  integer, intent(in)  :: n
  real(8), intent(in)  :: grid(n)
  integer, intent(in)  :: neighbours

  real(8), intent(out) :: grid_smooth(n)

  integer              :: i
  real(8)              :: divisor


  ! main averages along array
  divisor = 1.D0 / (neighbours*2+1)
  do i=neighbours+1,n-neighbours
    grid_smooth(i) = SUM(grid(i-neighbours:i+neighbours)) * divisor
  enddo
  ! tails
  do i=1,neighbours-1
    divisor = 1.D0 / (i*2+1)
    grid_smooth(i+1) = SUM(grid(:i*2+1)) * divisor
    grid_smooth(n-i) = SUM(grid(n-i*2:)) * divisor
  enddo
  ! extremes
  grid_smooth(1) = grid(1)
  grid_smooth(n) = grid(n)

end subroutine


! ======================= 'lowess' interpolation ======================
! local weighted scatterplot smoothing with gaussian weighting
! probably based on code by Jean-Pierre Moreau
! http://jean-pierre.moreau.pagesperso-orange.fr/f_function.html

subroutine grid2point_lowess(coord, grid, Z, span, n, point)

  !--------------------------------------------------------------------
  ! Local regression from grid to a single point
  !--------------------------------------------------------------------

  implicit none

  integer, intent(in)  :: n
  real(8), intent(in)  :: coord(2), grid(n,2), Z(n), span

  real(8), intent(out) :: point

  integer              :: i
  real(8), parameter   :: thr = 1.D-5     ! consider zero if lower
  real(8)              :: x, y
  real(8)              :: span_inv, w(n)


  ! invert smoothing parameter number
  span_inv = 1.D0 / span

  ! iterate over all the grid points
  do i=1, n
    ! calculate absolute diferences of distances
    x = ABS((grid(i,1)-coord(1))*span_inv)
    y = ABS((grid(i,2)-coord(2))*span_inv)
    ! calculate weights
    if (x>y) then
      w(i) = EXP( - ( x*SQRT(1.D0+y**2/(x**2) )**2 ))
    elseif (y<thr) then
      w(i) = 0.D0
    else
      w(i) = EXP( - ( y*SQRT(1.D0+x**2/(y**2) )**2 ))
    endif
  enddo

  ! calculate point
  point = SUM(Z*w) / SUM(w)

end subroutine

subroutine grid2grid_lowess(grid, Z, grid2, span, n, m, Zf)

  !--------------------------------------------------------------------
  ! Local regression from grid to grid
  !--------------------------------------------------------------------

  implicit none

  integer, intent(in)  :: n, m
  real(8), intent(in)  :: grid(n,2), Z(n), grid2(m,2), span

  real(8), intent(out) :: Zf(m)

  integer              :: i,j
  real(8), parameter   :: thr = 1.D-5     ! consider zero if lower
  real(8)              :: x, y
  real(8)              :: span_inv, w(n)


  ! invert smoothing parameter number
  span_inv = 1.D0 / span

  ! interpolate every point of the final grid
  do j=1, m
    ! iterate over all the grid points
    do i=1, n
      ! calculate absolute diferences of distances
      x = ABS((grid(i,1)-grid2(j,1))*span_inv)
      y = ABS((grid(i,2)-grid2(j,2))*span_inv)
      ! calculate weights
      if (x>y) then
        w(i) = EXP( - ( x*SQRT(1.D0+y**2/(x**2) )**2 ))
      elseif (y<thr) then
        w(i) = 0.D0
      else
        w(i) = EXP( - ( y*SQRT(1.D0+x**2/(y**2) )**2 ))
      endif
    enddo

    ! calculate point
    Zf(j) = SUM(Z*w) / SUM(w)

  enddo

end subroutine