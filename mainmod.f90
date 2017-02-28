module evolution
  implicit none
  contains
!-------------------------------------------------------------------------------
! Randomly choose two parents
  subroutine choose_parents(p,N,np,parent1,parent2,T,U)
    implicit none
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: parent1, parent2, U
    real(8)                            :: randnum, T
    integer                            :: i, N, np, index1, index2
    logical                            :: p1, p2

    p1 = .FALSE.
    p2 = .FALSE.
    do while(.not. p1)
      call random_number(randnum)
        index1 = ceiling(randnum * np)
      call random_number(randnum)
      if(randnum < exp(-U(index1)/T)) then
        do i = 1, N
          parent1(i) = p(index1,i)
        end do
        p1 = .TRUE.
      end if
    end do

    do while(.not. p2)
      call random_number(randnum)
      index2 = ceiling(randnum*np)
      if(index2 == index1) cycle
      call random_number(randnum)
      if(randnum < exp(-U(index2)/T)) then
        do i = 1, N
          parent2(i) = p(index2,i)
        end do
        p2 = .TRUE.
      end if
    end do

  end subroutine choose_parents
!-------------------------------------------------------------------------------
  ! Produce child with the randomly chosen parents
  subroutine produce_child(parent1,parent2,child,N)
    implicit none
    real(8),allocatable,dimension(:) :: parent1, parent2, child
    real(8)                          :: randnum
    integer                          :: N, index, i

    ! Determine point to split parents
    call random_number(randnum)
    index = ceiling(randnum*N)
    if(mod(index,3)/=0) then
      if(mod(index,3) == 1) then
        index = index + 2
      else
        index = index + 1
      end if
    end if

    ! Ensure the cut isn't at the last atom
    if(index == N) index = index - 3

    do i = 1, index + 1
      child(i) = parent1(i)
    end do

    do i = index + 1, N
      child(i) = parent2(i)
    end do

  end subroutine produce_child
!-------------------------------------------------------------------------------
! Produces a child with a mutation by moving it randomly between 5 and 25 times
! a random distance between 0 and 1 in a random direction
  subroutine mutate_child(child,N,edge)
    implicit none
    real(8),allocatable,dimension(:) :: child
    real(8)                          :: randnum, edge
    integer                          :: index, i, j, N

    call random_number(randnum)
    index = ceiling(randnum*21) + 4
    do i = 1, index
      do j = 1, N
        call random_number(randnum)
        child(j) = child(j) + (1.d0/3.d0 * randnum - 1.d0/6.d0)
        call pbc(child,edge,j)
      end do
    end do

  end subroutine mutate_child
!-------------------------------------------------------------------------------
! Decide whether or not to keep the child and add it to the population or
! discard it. Population size is constant
  subroutine keep_child(p,child,U,U_child,N,np,dE)
    implicit none
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: child, U
    real(8)                            :: U_child, dE
    integer                            :: N, np, i
    logical                            :: diff

    diff = .TRUE.
    if(U_child < maxval(U) - dE) then
      do i = 1, np
        if(abs(U_child - U(i)) < dE) diff = .FALSE.
      end do
      if(diff) then
        do i = 1, N
          p(maxloc(U),i) = child(i)
        end do
        print*, 'Child kept'
      end if
      U(maxloc(U)) = U_child
    end if

  end subroutine keep_child
!-------------------------------------------------------------------------------
! Moves clusters back to the origin
  subroutine pbc(child,edge,i)
    implicit none
    real(8),allocatable,dimension(:) :: child
    real(8)                          :: edge
    integer                          :: i

    if(child(i) < 0) child(i) = child(i) + edge
    if(child(i) > edge) child(i) = child(i) - edge

  end subroutine pbc
!-------------------------------------------------------------------------------
! Calculates the cluster energy of the population with the Lennard Jones
! potential
  double precision function U_func(p,N,i)
    implicit none
    real(8),allocatable,dimension(:,:) :: p
    real(8)                            :: distance
    integer                            :: N, i, j, k

    !Reset energy before every calculation
    U_func = 0.d0
    do j = 1, (N/3)-1
      do k = j+1, N/3
        distance = sqrt((p(i,1+(j-1)*3)-p(i,1+(k-1)*3))**2 &
              + (p(i,2+(j-1)*3) - p(i,2+(k-1)*3))**2 &
              + (p(i,3+(j-1)*3) - p(i,3+(k-1)*3))**2)
        U_func = U_func + (distance ** (-12) - distance ** (-6))
      end do
    end do
    U_func = 4 * U_func
  end function U_func
!-------------------------------------------------------------------------------
! Calculates the cluster energy of a gene with the Lennard Jones
! potential
  double precision function U_child_func(child,N)
    implicit none
    real(8),allocatable,dimension(:)   :: child
    real(8)                            :: distance
    integer                            :: N, j, k

    !Reset energy before every calculation
    U_child_func = 0.d0
    do j = 1, (N/3)-1
      do k = j+1, N/3
        distance = sqrt((child(1+(j-1)*3)-child(1+(k-1)*3))**2 &
              + (child(2+(j-1)*3) - child(2+(k-1)*3))**2 &
              + (child(3+(j-1)*3) - child(3+(k-1)*3))**2)
        U_child_func = U_child_func + (distance ** (-12) - distance ** (-6))
      end do
    end do
    U_child_func = 4 * U_child_func
  end function U_child_func
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Seed generation subroutines for random number generator
subroutine init_random_seed()
    implicit none
    integer                          :: i, n, clock
    integer,allocatable,dimension(:) :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i-1, i = 1, n) /)
    call random_seed(put = seed)

end subroutine

!-------------------------------------------------------------------------------
subroutine generate_seed(seed)
    implicit none
    integer :: seed, clock

    call system_clock(count=clock)
    seed = clock + 37

end subroutine

subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
      icount, numres, ifault )

!*****************************************************************************80
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, f )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
      real ( kind = 8 ) del
      real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
      real ( kind = 8 ), parameter :: eps = 0.001D+00
      real ( kind = 8 ), external :: fn
      integer ( kind = 4 ) i
      integer ( kind = 4 ) icount
      integer ( kind = 4 ) ifault
      integer ( kind = 4 ) ihi
      integer ( kind = 4 ) ilo
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jcount
      integer ( kind = 4 ) kcount
      integer ( kind = 4 ) konvge
      integer ( kind = 4 ) l
      integer ( kind = 4 ) numres
      real ( kind = 8 ) p(n,n+1)
      real ( kind = 8 ) p2star(n)
      real ( kind = 8 ) pbar(n)
      real ( kind = 8 ) pstar(n)
      real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
      real ( kind = 8 ) reqmin
      real ( kind = 8 ) rq
      real ( kind = 8 ) start(n)
      real ( kind = 8 ) step(n)
      real ( kind = 8 ) x
      real ( kind = 8 ) xmin(n)
      real ( kind = 8 ) y(n+1)
      real ( kind = 8 ) y2star
      real ( kind = 8 ) ylo
      real ( kind = 8 ) ynewlo
      real ( kind = 8 ) ystar
      real ( kind = 8 ) z
    !
    !  Check the input parameters.
    !
      if ( reqmin <= 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( n < 1 ) then
        ifault = 1
        return
      end if

      if ( konvge < 1 ) then
        ifault = 1
        return
      end if
    !
    !  Initialization.
    !
      icount = 0
      numres = 0
      jcount = konvge
      del = 1.0D+00
      rq = reqmin * real ( n, kind = 8 )
    !
    !  Initial or restarted loop.
    !
      do

        p(1:n,n+1) = start(1:n)
        y(n+1) = fn ( start )
        icount = icount + 1
    !
    !  Define the initial simplex.
    !
        do j = 1, n
          x = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          y(j) = fn ( start )
          icount = icount + 1
          start(j) = x
        end do
    !
    !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
    !  the vertex of the simplex to be replaced.
    !
        ilo = minloc ( y(1:n+1), 1 )
        ylo = y(ilo)
    !
    !  Inner loop.
    !
        do while ( icount < kcount )
    !
    !  YNEWLO is, of course, the HIGHEST value???
    !
          ihi = maxloc ( y(1:n+1), 1 )
          ynewlo = y(ihi)
    !
    !  Calculate PBAR, the centroid of the simplex vertices
    !  excepting the vertex with Y value YNEWLO.
    !
          do i = 1, n
            pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
          end do
    !
    !  Reflection through the centroid.
    !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar = fn ( pstar )
          icount = icount + 1
    !
    !  Successful reflection, so extension.
    !
          if ( ystar < ylo ) then

            p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
            y2star = fn ( p2star )
            icount = icount + 1
    !
    !  Retain extension or contraction.
    !
            if ( ystar < y2star ) then
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
            else
              p(1:n,ihi) = p2star(1:n)
              y(ihi) = y2star
            end if
    !
    !  No extension.
    !
          else

            l = 0
            do i = 1, n + 1
              if ( ystar < y(i) ) then
                l = l + 1
              end if
            end do

            if ( 1 < l ) then

              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
    !
    !  Contraction on the Y(IHI) side of the centroid.
    !
            else if ( l == 0 ) then

              p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Contract the whole simplex.
    !
              if ( y(ihi) < y2star ) then

                do j = 1, n + 1
                  p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                  xmin(1:n) = p(1:n,j)
                  y(j) = fn ( xmin )
                  icount = icount + 1
                end do

                ilo = minloc ( y(1:n+1), 1 )
                ylo = y(ilo)

                cycle
    !
    !  Retain contraction.
    !
              else
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              end if
    !
    !  Contraction on the reflection side of the centroid.
    !
            else if ( l == 1 ) then

              p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Retain reflection?
    !
              if ( y2star <= ystar ) then
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              else
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
              end if

            end if

          end if
    !
    !  Check if YLO improved.
    !
          if ( y(ihi) < ylo ) then
            ylo = y(ihi)
            ilo = ihi
          end if

          jcount = jcount - 1

          if ( 0 < jcount ) then
            cycle
          end if
    !
    !  Check to see if minimum reached.
    !
          if ( icount <= kcount ) then

            jcount = konvge

            x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
            z = sum ( ( y(1:n+1) - x )**2 )

            if ( z <= rq ) then
              exit
            end if

          end if

        end do
    !
    !  Factorial tests to check that YNEWLO is a local minimum.
    !
        xmin(1:n) = p(1:n,ilo)
        ynewlo = y(ilo)

        if ( kcount < icount ) then
          ifault = 2
          exit
        end if

        ifault = 0

        do i = 1, n
          del = step(i) * eps
          xmin(i) = xmin(i) + del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) - del - del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) + del
        end do

        if ( ifault == 0 ) then
          exit
        end if
    !
    !  Restart the procedure.
    !
        start(1:n) = xmin(1:n)
        del = eps
        numres = numres + 1

      end do

      return
    end

end module
