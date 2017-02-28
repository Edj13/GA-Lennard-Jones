program genetic
  implicit none
  real(8),allocatable,dimension(:,:) :: p, p_child
  real(8),allocatable,dimension(:)   :: U, parent1, parent2, child, fitness
  real(8),allocatable,dimension(:)   :: U_child, wa, H, G, W, prevE
  real(8)                            :: randnum, mu, potential, epsil
  integer                            :: np, N, i, j, steps
  integer                            :: maxsteps, nmat, f, ier, equalprev
  character(len=25)                  :: filename, tempN

  call init_random_seed()
  np = 10
  N = 7 * 3
  nmat = int(8d-1*np)
  mu = 2.d-1
  maxsteps = 1000000
  epsil = 1.d-3

  open(80,file='maxpotential.csv')
  open(90,file='minpotential.csv')

  allocate(p(np,N))
  allocate(p_child(nmat,N))
  allocate(U(np))
  allocate(U_child(nmat))
  allocate(parent1(N))
  allocate(parent2(N))
  allocate(child(N))
  allocate(fitness(np))
  allocate(wa(N))
  allocate(H(N*(N+1)/2))
  allocate(G(N))
  allocate(W(3*N))
  allocate(prevE(maxsteps))

  ! Randomly place each gene in computational box
  do i = 1, np
    do j = 1, N
      call random_number(randnum)
      p(i,j) = 3.5d0 * randnum
    end do
  end do

  ! Calculate initial cluster energies and relax to nearest local minimum
  do i = 1, np
    call LJpop(p,U,N,i)
    wa(:) = p(i,:)
    call zxmin(N,3,1000,0,wa,H,G,potential,W,ier)
    p(i,:) = wa(:)
    U(i) = potential
  end do

  do steps = 1, maxsteps
    call determine_fitness(np,fitness,U)
    do f = 1, nmat
      call choose_parents(p,N,np,parent1,parent2,fitness)
      call produce_child(parent1,parent2,child,N)
      p_child(f,:) = child
    end do

    do f = 1, nmat
      call random_number(randnum)
      if(randnum < mu) call mutate_child(p_child,N,f)
    end do


    ! Calculate cluster energy of childen and relax to nearest local minimum
    do i = 1, nmat
      call LJpop(p_child,U_child,N,i)
      wa(:) = p_child(i,:)
      call zxmin(N,5,10000000,0,wa,H,G,potential,W,ier)
      p_child(i,:) = wa(:)
      U_child(:) = potential
    end do

    ! Choose the lowest energies in the current generation and the children
    ! for inclusion into the next generation
    call keep_child(p,p_child,U,U_child,N,np,nmat)
    print*, "Max: ", maxval(U), "Min: ", minval(U), steps
    write(80,*) steps, maxval(U)
    write(90,*) steps, minval(U)
    prevE(steps) = minval(U)
    if (steps /= 1 .and. abs(prevE(steps-1) - minval(U)) < epsil) then
      equalprev = equalprev + 1
    else
      equalprev = 0
    end if

    !if (equalprev == 1000) exit
    if(abs(maxval(U) - minval(U)) < epsil) exit
  end do

  write(tempN,*) N/3
  filename = trim(adjustl(tempN)) // '.xyz'
  open(60, file=filename)
  write(60,*) n / 3.0
  write(60,*) 'Comment line'
  do i = 1, (n/3)
    write(60,"(A1,3F9.2)") 'C', p(minloc(U),1+(i-1)*3), &
                                 p(minloc(U),2+(i-1)*3), &
                                 p(minloc(U),3+(i-1)*3)
  end do
  close(60)

contains
!-------------------------------------------------------------------------------
! Determine fitness of each gene in the population
  subroutine determine_fitness(np,fitness,U)
    implicit none
    real(8),allocatable,dimension(:)   :: fitness, U
    real(8)                            :: rho
    integer                            :: np, i

    rho = 3.d0

    do i = 1, np
      fitness(i) = exp(-rho*((U(i) - minval(U))/(maxval(U) - minval(U))))
    end do
  end subroutine determine_fitness
!-------------------------------------------------------------------------------
! Randomly choose two parents
  subroutine choose_parents(p,N,np,parent1,parent2,fitness)
    implicit none
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: parent1, parent2, fitness
    real(8)                            :: randnum
    integer                            :: i, N, np, index1, index2
    logical                            :: p1, p2

    p1 = .FALSE.
    p2 = .FALSE.
    do while(.not. p1)
      call random_number(randnum)
      index1 = ceiling(randnum * np)
      call random_number(randnum)
      if(randnum < fitness(index1)) then
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
      if(randnum < fitness(index2)) then
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
    real(8),allocatable,dimension(:) :: parent1, parent2, child, Lp1, Lp2
    real(8),allocatable,dimension(:) :: distance, wa
    real(8)                          :: randnum, x, tempd, operator
    integer                          :: N, index, S, i, j, k
    integer                          :: atomsleft

    ! Randomly choose an operator to create a child
    call random_number(operator)

    if (operator < (1.d0/3.d0)) then
      allocate(Lp1(N))
      allocate(Lp2(N))
      allocate(wa(N/3))
      allocate(distance(N/3))

      ! Determine cut point
      call random_number(randnum)
      index = int(randnum*N)
      if(mod(index-1,3)/=0) then
        if(mod(index-1,3) == 1) then
          index = index + 2
        else
          index = index + 1
        end if
      end if

      ! Calculate distance between each atom in parent 1 and the cut point.
      do i = 1, (N/3)
        distance(i) = sqrt((parent1(1+(i-1)*3)-parent1(index))**2 &
              + (parent1(2+(i-1)*3) - parent1(index+1))**2 &
              + (parent1(3+(i-1)*3) - parent1(index+2))**2)
      end do

      ! Order distance and Lp1 from smallest to largest
      do i = 1, (N/3)
        index = minloc(distance,dim=1)
        Lp1(1+(i-1)*3) = parent1(index)
        Lp1(2+(i-1)*3) = parent1(index+1)
        Lp1(3+(i-1)*3) = parent1(index+2)
        distance(index) = huge(x)
      end do

      call random_number(randnum)
      S = int(1 + randnum * (real(N)/3.d0-1))

      do i = 1, S
        child(1+(i-1)*3) = Lp1(1+(i-1)*3)
        child(2+(i-1)*3) = Lp1(2+(i-1)*3)
        child(3+(i-1)*3) = Lp1(3+(i-1)*3)
      end do

      ! Calculate distance between each atom in parent 2 and the cut point.
      do i = 1, (N/3)
        distance(i) = sqrt((parent2(1+(i-1)*3)-parent1(index))**2 &
              + (parent2(2+(i-1)*3) - parent1(index+1))**2 &
              + (parent2(3+(i-1)*3) - parent1(index+2))**2)
      end do

      ! Order distance and Lp2 from smallest to largest
      do i = 1, (N/3)
        index = minloc(distance,dim=1)
        Lp2(1+(i-1)*3) = parent1(index)
        Lp2(2+(i-1)*3) = parent1(index+1)
        Lp2(3+(i-1)*3) = parent1(index+2)
        distance(index) = huge(x)
      end do

      ! Remove from Lp2 atoms that are too close to particles already copied to
      ! the child
      atomsleft = (N/3)
      do i = 1, S
        do j = 1, (N/3)
          tempd = sqrt((child(1+(i-1)*3) - parent2(1+(j-1)*3))**2 + &
                       (child(2+(i-1)*3) - parent2(2+(j-1)*3))**2 + &
                       (child(3+(i-1)*3) - parent2(3+(j-1)*3))**2)
          if (tempd < 0.5) then
            do k = 1, (N/3) - 1
              Lp2(1+(k-1)*3) = Lp2(4+(k-1)*3)
              Lp2(2+(k-1)*3) = Lp2(5+(k-1)*3)
              Lp2(3+(k-1)*3) = Lp2(6+(k-1)*3)
            end do
            atomsleft = atomsleft - 1
          end if
        end do
      end do

      do i = 1, (N/3-S)
        if (atomsleft <= 0) then
          call random_number(randnum)
          child(1+(i+S)*3-3) = (N/3) ** (1.d0/3.d0) * randnum
          call random_number(randnum)
          child(2+(i+S)*3-3) = (N/3) ** (1.d0/3.d0) * randnum
          call random_number(randnum)
          child(3+(i+S)*3-3) = (N/3) ** (1.d0/3.d0) * randnum
          cycle
        end if
        child(1+(i+S)*3-3) = Lp2(1+(i-1)*3)
        child(2+(i+S)*3-3) = Lp2(2+(i-1)*3)
        child(3+(i+S)*3-3) = Lp2(3+(i-1)*3)
        atomsleft = atomsleft - 1
      end do
    else if (operator >= (1.d0/3.d0) .and. operator < (2.d0/3.d0)) then
      do i = 1, N
        child(i) = (parent1(i) + parent2(i)) / 2.d0
      end do
    else
      do i = 1, N
        child(i) = sqrt(parent1(i)*parent2(i))
      end do
    end if

  end subroutine produce_child
!-------------------------------------------------------------------------------
! Produces a child with a mutation by moving it randomly between 5 and 25 times
! a random distance between 0 and 1 in a random direction
  subroutine mutate_child(p_child,N,f)
    implicit none
    real(8),allocatable,dimension(:,:) :: p_child
    real(8)                            :: randnum
    integer,allocatable,dimension(:)   :: used
    integer                            :: j, k, N, f

    allocate(used(ceiling((real(N)/3.d0))))
    used = 0
    do k = 1, N/3
      45 call random_number(randnum)
      j = ceiling(randnum * (real(N)/3.d0))
      if(any(used == j)) go to 45
      used(k) = j
      call random_number(randnum)
      p_child(f,1+(j-1)*3) = 3.5d0 * randnum
      call random_number(randnum)
      p_child(f,2+(j-1)*3) = 3.5d0 * randnum
      call random_number(randnum)
      p_child(f,3+(j-1)*3) = 3.5d0 * randnum
    end do
    deallocate(used)
  end subroutine mutate_child
!-------------------------------------------------------------------------------
! Decide whether or not to keep the child and add it to the population or
! discard it. Population size is constant
  subroutine keep_child(p,p_child,U,U_child,N,np,nmat)
    implicit none
    real(8),allocatable,dimension(:,:) :: p, p_child, wa
    real(8),allocatable,dimension(:)   :: U, U_child, wb
    real(8)                            :: x
    integer                            :: N, np, i, index, nmat

    allocate(wa(np+nmat,N))
    allocate(wb(np+nmat))

    wb = huge(x)
    do i = 1, np
      do j = 1, N
        wa(i,j) = p(i,j)
      end do
      wb(i) = U(i)
    end do

    do i = np+1, nmat + np
      do j = 1, N
        wa(i,j) = p_child(i-np,j)
      end do
      wb(i) = U_child(i-np)
    end do
    !print*, 'All energies:', wb
    do i = 1, np
      index = minloc(wb,dim=1)
      p(i,:) = wa(index,:)
      U(i) = wb(index)
      wb(index) = huge(x)
    end do
    !print*, 'Kept energies', U

  end subroutine keep_child
!-------------------------------------------------------------------------------
! Calculates the cluster energy of the population with the Lennard Jones
! potential
  subroutine LJpop(p,U,N,i)
    implicit none
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: U
    real(8)                            :: distance
    integer                            :: N, i, j, k

    !Reset energy before every calculation
    U(i) = 0.d0
    do j = 1, (N/3)-1
      do k = j+1, N/3
        distance = sqrt((p(i,1+(j-1)*3)-p(i,1+(k-1)*3))**2 &
              + (p(i,2+(j-1)*3) - p(i,2+(k-1)*3))**2 &
              + (p(i,3+(j-1)*3) - p(i,3+(k-1)*3))**2)
        U(i) = U(i) + (distance ** (-12) - distance ** (-6))
      end do
    end do
    U(i) = 4 * U(i)
  end subroutine LJpop

!-------------------------------------------------------------------------------
end program genetic


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

!-------------------------------------------------------------------------------
! Calculates the cluster energy of a gene with the Lennard Jones
! potential
  subroutine funct(N,child,U_child)
    implicit none
    real(8),dimension(7*3)           :: child
    real(8)                            :: U_child, distance
    integer                            :: N, j, k

    !Reset energy before every calculation
    U_child = 0.d0
    do j = 1, (N/3)-1
      do k = j+1, N/3
        distance = sqrt((child(1+(j-1)*3)-child(1+(k-1)*3))**2 &
              + (child(2+(j-1)*3) - child(2+(k-1)*3))**2 &
              + (child(3+(j-1)*3) - child(3+(k-1)*3))**2)
        U_child = U_child + (distance ** (-12) - distance ** (-6))
      end do
    end do
    U_child = 4 * U_child
  end subroutine funct

!**************************************************************
       Subroutine ZXMIN (N,NSIG,MAXFN,IOPT,X,H,G,F,W,IER)
!   IMSL ROUTINE NAME   - ZXMIN
!
!--  ---------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - MINIMUM OF A FUNCTION OF N VARIABLES USING
!                           A QUASI-NEWTON METHOD
!
!   USAGE               - CALL ZXMIN (FUNCT,N,NSIG,MAXFN,IOPT,X,H,G,F,
!                           W,IER)
!
!   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
!                           THE FUNCTION F FOR GIVEN PARAMETER VALUES
!                           X(1),X(2),...,X(N).
!                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
!                           CALL FUNCT(N,X,F)
!                           WHERE X IS A VECTOR OF LENGTH N.
!                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT
!                           IN THE CALLING PROGRAM. FUNCT MUST NOT
!                           ALTER THE VALUES OF X(I),I=1,...,N OR N.
!                N      - THE NUMBER OF PARAMETERS (I.E., THE LENGTH
!                           OF X) (INPUT)
!                NSIG   - CONVERGENCE CRITERION. (INPUT). THE NUMBER
!                           OF DIGITS OF ACCURACY REQUIRED IN THE
!                           PARAMETER ESTIMATES.
!                           THIS CONVERGENCE CONDITION IS SATISIFIED IF
!                           ON TWO SUCCESSIVE ITERATIONS, THE PARAMETER
!                           ESTIMATES (I.E.,X(I), I=1,...,N) AGREE,
!                           COMPONENT BY COMPONENT, TO NSIG DIGITS.
!                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
!                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT)
!                IOPT   - OPTIONS SELECTOR. (INPUT)
!                         IOPT = 0 CAUSES ZXMIN TO INITIALIZE THE
!                           HESSIAN MATRIX H TO THE IDENTITY MATRIX.
!                         IOPT = 1 INDICATES THAT H HAS BEEN INITIALIZED
!                           BY THE USER TO A POSITIVE DEFINITE MATRIX.
!                         IOPT = 2 CAUSES ZXMIN TO COMPUTE THE DIAGONAL
!                           VALUES OF THE HESSIAN MATRIX AND SET H TO
!                           A DIAGONAL MATRIX CONTAINING THESE VALUES.
!                         IOPT = 3 CAUSES ZXMIN TO COMPUTE AN ESTIMATE
!                           OF THE HESSIAN IN H.
!                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
!                           VALUES.
!                         ON INPUT, X MUST CONTAIN THE INITIAL
!                           PARAMETER ESTIMATES.
!                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER
!                           ESTIMATES AS DETERMINED BY ZXMIN.
!                H      - VECTOR OF LENGTH N*(N+1)/2 CONTAINING AN
!                           ESTIMATE OF THE HESSIAN MATRIX
!                           D**2F/(DX(I)DX(J)), I,J=1,...,N.
!                           H IS STORED IN SYMMETRIC STORAGE MODE.
!                         ON INPUT, IF IOPT = 0, 2, OR 3 ZXMIN INITIA-
!                           LIZES H. AN INITIAL SETTING OF H BY THE
!                           USER IS INDICATED BY IOPT=1.
!                           H MUST BE POSITIVE DEFINITE. IF IT IS NOT,
!                           A TERMINAL ERROR OCCURS.
!                         ON OUTPUT, H CONTAINS AN ESTIMATE OF THE
!                           HESSIAN AT THE FINAL PARAMETER ESTIMATES
!                           (I.E., AT X(1),X(2),...,X(N))
!                G      - A VECTOR OF LENGTH N CONTAINING AN ESTIMATE
!                           OF THE GRADIENT DF/DX(I),I=1,...,N AT THE
!                           FINAL PARAMETER ESTIMATES. (OUTPUT)
!                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
!                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
!                W      - A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
!                         ON OUTPUT, WORK(I), CONTAINS FOR
!                           I = 1, THE NORM OF THE GRADIENT (I.E.,
!                             SQRT(G(1)**2+G(2)**2+...+G(N)**2))
!                           I = 2, THE NUMBER OF FUNCTION EVALUATIONS
!                             PERFORMED.
!                           I = 3, AN ESTIMATE OF THE NUMBER OF
!                             SIGNIFICANT DIGITS IN THE FINAL
!                             PARAMETER ESTIMATES.
!                IER    - ERROR PARAMETER (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 IMPLIES THAT THE INITIAL HESSIAN
!                             USED BY ZXMIN IS NOT POSITIVE DEFINITE,
!                             EVEN AFTER ADDING A MULTIPLE OF THE
!                             IDENTITY TO MAKE ALL DIAGONAL ELEMENTS
!                             POSITIVE.
!                           IER = 130 IMPLIES THAT THE ITERATION WAS
!                             TERMINATED DUE TO ROUNDING ERRORS
!                             BECOMING DOMINANT. THE PARAMETER
!                             ESTIMATES HAVE NOT BEEN DETERMINED TO
!                             NSIG DIGITS.
!                           IER = 131 IMPLIES THAT THE ITERATION WAS
!                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO,ZXMJN
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!--  ---------------------------------------------------------------------
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,MAXFN,IOPT,IER
      DOUBLE PRECISION   X(N),G(N),H(1),F,W(1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER IG,IGG,IS,IDIFF,IR,IJ,I,J,NM1,JJ,JP1,L,KJ,K,IFN,LINK,ITN,II,IM1,JNT,NP1,JB,NJ
      DOUBLE PRECISION   REPS,AX,ZERO,ONE,HALF,SEVEN,FIVE,TWELVE,TEN,HH
      DOUBLE PRECISION   EPS,HJJ,V,DF,RELX,GS0,DIFF,AEPS,ALPHA,FF,TOT
      DOUBLE PRECISION   F1,F2,Z,GYS,DGS,SIG,ZZ,GNRM,P1,HHH,GHH,H2,F11
      DOUBLE PRECISION   F12,F21,F22,HMAX,HMIN
      DATA               REPS/2.775557562D-17/,AX/0.1D0/
      DATA ZERO/0.0D0/,ONE/1.0D0/,HALF/0.5D0/,SEVEN/7.0D0/,FIVE/5.0D0/,TWELVE/12.0D0/,TEN/10.0D0/,P1/0.1D0/
!                                  INITIALIZATION
!                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      HH = DSQRT(REPS)
      H2 = DSQRT(HH)
      EPS = TEN**(-NSIG)
      IG = N
      IGG = N+N
      IS = IGG
      IDIFF = 1
      IR = N
      W(1) = -ONE
      W(2) = ZERO
      W(3) = ZERO
!                                  EVALUATE FUNCTION AT STARTING POINT
      DO 5 I=1,N
         G(I) = X(I)
    5 CONTINUE
      CALL FUNCT(N,G,F)
      IFN = 1
      IF (IOPT.EQ.1) GO TO 50
!                                  SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
      IF (N.EQ.1) GO TO 20
      IJ = 2
      DO 15 I=2,N
         DO 10 J=2,I
            H(IJ) = ZERO
            IJ = IJ+1
   10    CONTINUE
         IJ = IJ+1
   15 CONTINUE
   20 IF (IOPT.NE.0) GO TO 30
!                                  SET DIAGONAL ELEMENTS OF H TO ONE
      IJ = 0
      DO 25 I=1,N
         IJ = IJ+I
         H(IJ) = ONE
   25 CONTINUE
      GO TO 95
!                                  GET DIAGONAL ELEMENTS OF HESSIAN
   30 IM1 = 1
      NM1 = 1
      NP1 = N+1
      DO 35 I=2,NP1
         HHH = H2*DMAX1(DABS(X(IM1)),AX)
         G(IM1) = X(IM1)+HHH
         CALL FUNCT(N,G,F2)
         G(IM1) = X(IM1)-HHH
         CALL FUNCT(N,G,FF)
         H(NM1) = (FF-F+F2-F)/(HHH*HHH)
         G(IM1) = X(IM1)
         IM1 = I
         NM1 = I+NM1
   35 CONTINUE
      IFN = IFN+N+N
      IF (IOPT.NE.3 .OR. N.EQ.1) GO TO 50
!                                  GET THE REST OF THE HESSIAN
      JJ = 1
      II = 2
      DO 45 I=2,N
         GHH = H2*DMAX1(DABS(X(I)),AX)
         DO 40 J=1,JJ
            HHH = H2*DMAX1(DABS(X(J)),AX)
            G(I) = X(I)+GHH
            G(J) = X(J)+HHH
            CALL FUNCT(N,G,F22)
            G(I) = X(I)-GHH
            CALL FUNCT(N,G,F12)
            G(J) = X(J)-HHH
            CALL FUNCT(N,G,F11)
            G(I) = X(I)+GHH
            CALL FUNCT(N,G,F21)
            H(II) = (F22-F21-F12+F11)/(4.D0*HHH*GHH)
            G(J) = X(J)
            II = II+1
   40    CONTINUE
         G(I) = X(I)
         JJ = JJ+1
         II = II+1
   45 CONTINUE
      IFN = IFN+((N*N-N)*2)
!                                  ADD MULTIPLE OF IDENTITY TO
!                                  MAKE DIAGONAL ELEMENTS POSITIVE
   50 HMIN = H(1)
      HMAX = H(1)
      NM1 = 1
      DO 55 I=1,N
         HMIN = DMIN1(HMIN,H(NM1))
         HMAX = DMAX1(HMAX,H(NM1))
         NM1 = NM1+I+1
   55 CONTINUE
      HMIN = DMAX1(0.01D0*(DABS(HMAX)+DABS(HMIN))-HMIN,0.0D0)
      NM1 = 1
      DO 60 I=1,N
         H(NM1) = H(NM1)+HMIN
         NM1 = NM1+I+1
   60 CONTINUE
!                                  FACTOR H TO L*D*L-TRANSPOSE
      IR = N
      IF (N.GT.1) GO TO 65
      IF (H(1).GT.ZERO) GO TO 95
      H(1) = ZERO
      IR = 0
      GO TO 90
   65 NM1 = N-1
      JJ = 0
      DO 85 J=1,N
         JP1 = J+1
         JJ = JJ+J
         HJJ = H(JJ)
         IF (HJJ.GT.ZERO) GO TO 70
         H(JJ) = ZERO
         IR = IR-1
         GO TO 85
   70    IF (J.EQ.N) GO TO 85
         IJ = JJ
         L = 0
         DO 80 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)/HJJ
            KJ = IJ
            DO 75 K=I,N
               H(KJ+L) = H(KJ+L)-H(KJ)*V
               KJ = KJ+K
   75       CONTINUE
            H(IJ) = V
   80    CONTINUE
   85 CONTINUE
   90 IF (IR.EQ.N) GO TO 95
      IER = 129
      GO TO 9000
   95 ITN = 0
      DF = -ONE
!                                  EVALUATE GRADIENT W(IG+I),I=1,...,N
  100 LINK = 1
      GO TO 280
  105 CONTINUE
!                                  BEGIN ITERATION LOOP
      IF (IFN.GE.MAXFN) GO TO 240
      ITN = ITN+1
      DO 110 I=1,N
         W(I) = -W(IG+I)
  110 CONTINUE
!                                  DETERMINE SEARCH DIRECTION W
!                                    BY SOLVING H*W = -G WHERE
!                                    H = L*D*L-TRANSPOSE
      IF (IR.LT.N) GO TO 140
!                                  N .EQ. 1
      G(1) = W(1)
      IF (N.GT.1) GO TO 115
      W(1) = W(1)/H(1)
      GO TO 140
!                                  N .GT. 1
  115 II = 1
!                                  SOLVE L*W = -G
      DO 125 I=2,N
         IJ = II
         II = II+I
         V = W(I)
         IM1 = I-1
         DO 120 J=1,IM1
            IJ = IJ+1
            V = V-H(IJ)*W(J)
  120    CONTINUE
         G(I) = V
         W(I) = V
  125 CONTINUE
!                                  SOLVE (D*LT)*Z = W WHERE
!                                  LT = L-TRANSPOSE
      W(N) = W(N)/H(II)
      JJ = II
      NM1 = N-1
      DO 135 NJ=1,NM1
!                                  J = N-1,N-2,...,1
         J = N-NJ
         JP1 = J+1
         JJ = JJ-JP1
         V = W(J)/H(JJ)
         IJ = JJ
         DO 130 I=JP1,N
            IJ = IJ+I-1
            V = V-H(IJ)*W(I)
  130    CONTINUE
         W(J) = V
  135 CONTINUE
!                                  DETERMINE STEP LENGTH ALPHA
  140 RELX = ZERO
      GS0 = ZERO
      DO 145 I=1,N
         W(IS+I) = W(I)
         DIFF = DABS(W(I))/DMAX1(DABS(X(I)),AX)
         RELX = DMAX1(RELX,DIFF)
         GS0 = GS0+W(IG+I)*W(I)
  145 CONTINUE
      IF (RELX.EQ.ZERO) GO TO 245
      AEPS = EPS/RELX
      IER = 130
      IF (GS0.GE.ZERO) GO TO 245
      IF (DF.EQ.ZERO) GO TO 245
      IER = 0
      ALPHA = (-DF-DF)/GS0
      IF (ALPHA.LE.ZERO) ALPHA = ONE
      ALPHA = DMIN1(ALPHA,ONE)
      IF (IDIFF.EQ.2) ALPHA = DMAX1(P1,ALPHA)
      FF = F
      TOT = ZERO
      JNT = 0
!                                  SEARCH ALONG  X+ALPHA*W
  150 IF (IFN.GE.MAXFN) GO TO 240
      DO 155 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
  155 CONTINUE
      CALL FUNCT(N,W,F1)
      IFN = IFN+1
      IF (F1.GE.F) GO TO 180
      F2 = F
      TOT = TOT+ALPHA
  160 IER = 0
      F = F1
      DO 165 I=1,N
         X(I) = W(I)
  165 CONTINUE
      IF (JNT-1) 170, 200, 205
  170 IF (IFN.GE.MAXFN) GO TO 240
      DO 175 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
  175 CONTINUE
      CALL FUNCT(N,W,F1)
      IFN = IFN+1
      IF (F1.GE.F) GO TO 205
      IF (F1+F2.GE.F+F .AND. SEVEN*F1+FIVE*F2.GT.TWELVE*F) JNT = 2
      TOT = TOT+ALPHA
      ALPHA = ALPHA+ALPHA
      GO TO 160
  180 CONTINUE
      IF (F.EQ.FF .AND. IDIFF.EQ.2 .AND. RELX.GT.EPS) IER = 130
      IF (ALPHA.LT.AEPS) GO TO 245
      IF (IFN.GE.MAXFN) GO TO 240
      ALPHA = HALF*ALPHA
      DO 185 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
  185 CONTINUE
      CALL FUNCT(N,W,F2)
      IFN = IFN+1
      IF (F2.GE.F) GO TO 195
      TOT = TOT+ALPHA
      IER = 0
      F = F2
      DO 190 I=1,N
         X(I) = W(I)
  190 CONTINUE
      GO TO 200
  195 Z = P1
      IF (F1+F.GT.F2+F2) Z = ONE+HALF*(F-F1)/(F+F1-F2-F2)
      Z = DMAX1(P1,Z)
      ALPHA = Z*ALPHA
      JNT = 1
      GO TO 150
  200 IF (TOT.LT.AEPS) GO TO 245
  205 ALPHA = TOT
!                                  SAVE OLD GRADIENT
      DO 210 I=1,N
         W(I) = W(IG+I)
  210 CONTINUE
!                                  EVALUATE GRADIENT W(IG+I), I=1,...,N
      LINK = 2
      GO TO 280
  215 IF (IFN.GE.MAXFN) GO TO 240
      GYS = ZERO
      DO 220 I=1,N
         GYS = GYS+W(IG+I)*W(IS+I)
         W(IGG+I) = W(I)
  220 CONTINUE
      DF = FF-F
      DGS = GYS-GS0
      IF (DGS.LE.ZERO) GO TO 105
      IF (DGS+ALPHA*GS0.GT.ZERO) GO TO 230
!                                  UPDATE HESSIAN H USING
!                                    COMPLEMENTARY DFP FORMULA
      SIG = ONE/GS0
      IR = -IR
      CALL ZXMJN(H,N,W,SIG,G,IR,0,ZERO)
      DO 225 I=1,N
         G(I) = W(IG+I)-W(IGG+I)
  225 CONTINUE
      SIG = ONE/(ALPHA*DGS)
      IR = -IR
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      GO TO 105
!                                  UPDATE HESSIAN USING
!                                    DFP FORMULA
  230 ZZ = ALPHA/(DGS-ALPHA*GS0)
      SIG = -ZZ
      CALL ZXMJN(H,N,W,SIG,G,IR,0,REPS)
      Z = DGS*ZZ-ONE
      DO 235 I=1,N
         G(I) = W(IG+I)+Z*W(IGG+I)
  235 CONTINUE
      SIG = ONE/(ZZ*DGS*DGS)
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      GO TO 105
  240 IER = 131
!                                  MAXFN FUNCTION EVALUATIONS
      GO TO 250
  245 IF (IDIFF.EQ.2) GO TO 250
!                                  CHANGE TO CENTRAL DIFFERENCES
      IDIFF = 2
      GO TO 100
  250 IF (IER.NE.0) GO TO 255
      IF (RELX.LE.EPS) GO TO 255
      GO TO 100
!                                  MOVE GRADIENT TO G AND RETURN
  255 GNRM = ZERO
      DO 260 I=1,N
         G(I) = W(IG+I)
         GNRM = GNRM+G(I)*G(I)
  260 CONTINUE
      GNRM = DSQRT(GNRM)
      W(1) = GNRM
      W(2) = IFN
      W(3) = -DLOG10(DMAX1(REPS,RELX))
!                                  COMPUTE H = L*D*L-TRANSPOSE
      IF (N.EQ.1) GO TO 9000
      NP1 = N+1
      NM1 = N-1
      JJ = (N*(NP1))/2
      DO 275 JB=1,NM1
         JP1 = NP1-JB
         JJ = JJ-JP1
         HJJ = H(JJ)
         IJ = JJ
         L = 0
         DO 270 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)*HJJ
            KJ = IJ
            DO 265 K=I,N
               H(KJ+L) = H(KJ+L)+H(KJ)*V
               KJ = KJ+K
  265       CONTINUE
            H(IJ) = V
  270    CONTINUE
         HJJ = H(JJ)
  275 CONTINUE
      GO TO 9000
!                                  EVALUATE GRADIENT
  280 IF (IDIFF.EQ.2) GO TO 290
!                                  FORWARD DIFFERENCES
!                                    GRADIENT = W(IG+I), I=1,...,N
      DO 285 I=1,N
         Z = HH*DMAX1(DABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         CALL FUNCT(N,X,F1)
         W(IG+I) = (F1-F)/Z
         X(I) = ZZ
  285 CONTINUE
      IFN = IFN+N
      GO TO (105, 215), LINK
!                                  CENTRAL DIFFERENCES
!                                    GRADIENT = W(IG+I), I=1,...,N
  290 DO 295 I=1,N
         Z = HH*DMAX1(DABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         CALL FUNCT(N,X,F1)
         X(I) = ZZ-Z
         CALL FUNCT(N,X,F2)
         W(IG+I) = (F1-F2)/(Z+Z)
         X(I) = ZZ
  295 CONTINUE
      IFN = IFN+N+N
      GO TO (105, 215), LINK
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HZXMIN )
 9005 RETURN
      END
      SUBROUTINE UERTST (IER,NAME)
!   IMSL ROUTINE NAME   - UERTST
!
!--  ---------------------------------------------------------------------
!
!   COMPUTER            - VAX/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
!
!   USAGE               - CALL UERTST (IER,NAME)
!
!   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
!                           IER = I+J WHERE
!                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
!                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
!                             I =  32 IMPLIES WARNING MESSAGE.
!                             J = ERROR CODE RELEVANT TO CALLING
!                                 ROUTINE.
!                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
!                           THE NAME OF THE CALLING ROUTINE. (INPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UGETIO,USPKD
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
!                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
!                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
!                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
!                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
!                UGETIO AS FOLLOWS..
!                                NIN = 0
!                                NOUT = NEW OUTPUT UNIT NUMBER
!                                CALL UGETIO(3,NIN,NOUT)
!                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!--  ---------------------------------------------------------------------
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER            NAME(1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER I,IEQ,IEQDF,IOUNIT,LEVEL,LEVOLD,NAMEQ(6),NAMSET(6),NAMUPK(6),NIN,NMTB
      DATA               NAMSET/1HU,1HE,1HR,1HS,1HE,1HT/
      DATA               NAMEQ/6*1H /
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
      DATA               NAMUPK/1HZ,1HX,1HM,1HI,1HN,1H /
!                                  UNPACK NAME INTO NAMUPK
!                                  FIRST EXECUTABLE STATEMENT
!     CALL USPKD (NAME,6,NAMUPK,NMTB)
!                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
!                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
!                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
!                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
!                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
!                                  CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
!                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
!
!                                  SAVE P FOR P = R CASE
!                                    P IS THE PAGE NAMUPK
!                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 RETURN
      END
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
!   IMSL ROUTINE NAME   - UGETIO
!
!--  ---------------------------------------------------------------------
!
!   COMPUTER            - VAX/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
!                           VALUES FOR INPUT AND OUTPUT UNIT
!                           IDENTIFIERS.
!
!   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
!
!   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
!                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
!                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
!                           AND NOUT, RESPECTIVELY.
!                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
!                           RESET FOR SUBSEQUENT USE.
!                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
!                           RESET FOR SUBSEQUENT USE.
!                NIN    - INPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
!                NOUT   - OUTPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
!                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
!                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
!                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
!                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!--  ---------------------------------------------------------------------
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/35/,NOUTD/36/
!                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END
      SUBROUTINE ZXMJN  (A,N,Z,SIG,W,IR,MK,EPS)
!   IMSL ROUTINE NAME   - ZXMJN
!
!--  ---------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES ZXMIN AND
!                           ZXMWD
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - DOUBLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!--  ---------------------------------------------------------------------
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR,MK
      DOUBLE PRECISION   A(1),Z(N),SIG,W(N),EPS
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,JJ,IJ,JP1,I,II,MM
      DOUBLE PRECISION   ZERO,ONE,FOUR,TI,V,TIM,AL,R,B,GM,Y
      DATA               ZERO/0.0D0/,ONE/1.0D0/,FOUR/4.0D0/
!                                  UPDATE FACTORS GIVEN IN A
!                                    SIG*Z*Z-TRANSPOSE IS ADDED
!                                  FIRST EXECUTABLE STATEMENT
      IF (N.GT.1) GO TO 5
!                                  N .EQ. 1
      A(1) = A(1)+SIG*Z(1)*Z(1)
      IR = 1
      IF (A(1).GT.ZERO) GO TO 9005
      A(1) = ZERO
      IR = 0
      GO TO 9005
!                                  N .GT. 1
    5 IF (SIG.GT.ZERO) GO TO 65
      IF (SIG.EQ.ZERO.OR.IR.EQ.0) GO TO 9005
      TI = ONE/SIG
      JJ = 0
      IF (MK.EQ.0) GO TO 15
!                                  L*W = Z ON INPUT
      DO 10 J=1,N
         JJ = JJ+J
         IF (A(JJ).NE.ZERO) TI = TI+(W(J)*W(J))/A(JJ)
   10 CONTINUE
      GO TO 40
!                                  SOLVE L*W = Z
   15 DO 20 J=1,N
         W(J) = Z(J)
   20 CONTINUE
      DO 35 J=1,N
         JJ = JJ+J
         V = W(J)
         IF (A(JJ).GT.ZERO) GO TO 25
         W(J) = ZERO
         GO TO 35
   25    TI = TI+(V*V)/A(JJ)
         IF (J.EQ.N) GO TO 35
         IJ = JJ
         JP1 = J+1
         DO 30 I=JP1,N
            IJ = IJ+I-1
            W(I) = W(I)-V*A(IJ)
   30    CONTINUE
   35 CONTINUE
!                                  SET TI, TIM AND W
   40 IF (IR.LE.0) GO TO 45
      IF (TI.GT.ZERO) GO TO 50
      IF (MK-1) 65,65,55
   45 TI = ZERO
      IR = -IR-1
      GO TO 55
   50 TI = EPS/SIG
      IF (EPS.EQ.ZERO) IR = IR-1
   55 TIM = TI
      II = JJ
      I = N
      DO 60 J=1,N
         IF (A(II).NE.ZERO) TIM = TI-(W(I)*W(I))/A(II)
         W(I) = TI
         TI = TIM
         II = II-I
         I = I-1
   60 CONTINUE
      MM = 1
      GO TO 70
   65 MM = 0
      TIM = ONE/SIG
   70 JJ = 0
!                                  UPDATE A
      DO 110 J=1,N
         JJ = JJ+J
         IJ = JJ
         JP1 = J+1
!                                  UPDATE A(J,J)
         V = Z(J)
         IF (A(JJ).GT.ZERO) GO TO 85
!                                  A(J,J) .EQ. ZERO
         IF (IR.GT.0.OR.SIG.LT.ZERO.OR.V.EQ.ZERO) GO TO 80
         IR = 1-IR
         A(JJ) = (V*V)/TIM
         IF (J.EQ.N) GO TO 9005
         DO 75 I=JP1,N
            IJ = IJ+I-1
            A(IJ) = Z(I)/V
   75    CONTINUE
         GO TO 9005
   80    TI = TIM
         GO TO 110
!                                  A(J,J) .GT. ZERO
   85    AL = V/A(JJ)
         TI = W(J)
         IF (MM.EQ.0) TI = TIM+V*AL
         R = TI/TIM
         A(JJ) = R*A(JJ)
         IF (R.EQ.ZERO) GO TO 115
         IF (J.EQ.N) GO TO 115
!                                  UPDATE REMAINDER OF COLUMN J
         B = AL/TI
         IF (R.GT.FOUR) GO TO 95
         DO 90 I=JP1,N
            IJ = IJ+I-1
            Z(I) = Z(I)-V*A(IJ)
            A(IJ) = A(IJ)+B*Z(I)
   90    CONTINUE
         GO TO 105
   95    GM = TIM/TI
         DO 100 I=JP1,N
            IJ = IJ+I-1
            Y = A(IJ)
            A(IJ) = B*Z(I)+Y*GM
            Z(I) = Z(I)-V*Y
  100    CONTINUE
  105    TIM = TI
  110 CONTINUE
  115 IF (IR.LT.0) IR = -IR
 9005 CONTINUE
      RETURN
      END
