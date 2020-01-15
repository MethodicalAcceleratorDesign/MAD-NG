subroutine micit(a, xin, cin, res, nx, rms, im, ic, iter, ny, ax, cinx, &
                 xinx, resx, rho, ptop, rmss, xrms, xptp, xiter, ifail)
    implicit none
    ! ****************************************************
    !                                                    *
    !    Driving routine for MICADO correction           *
    !                                                    *
    !     Author: WFH  05.02.02                          *
    !                                                    *
    ! ****************************************************
    ! RMS = value of tolerance for correction
    integer :: im, ic, iter, ifail
    double precision :: a(im,ic), xin(im), cin(ic), res(im)
    integer :: nx(ic)
    double precision :: ax(im,ic), cinx(ic), xinx(im), resx(im), rho(3*ic), ptop(ic)
    double precision :: rmss(ic), xrms(ic), xptp(ic), xiter(ic)

!    integer :: j

    integer :: i, ny(ic)
    double precision :: rms, rm
    double precision, parameter :: zero=0d0

    write (*,'(/A, I5, A/)') 'start MICADO correction with ',iter,' correctors'

    AX(:im,:ic) = A(:im,:ic)
    CINX(:ic) = zero

    NY(1:ic) = (/ (i, i = 1, ic) /) ! NY(i) = i

    XINX(:im) = XIN(:im)
    RESX(:im) = zero

!  write(*,'("A[",i2,"x",i2,"]=")') im, ic
!  do i=1,im ; write(*,'(1x,10f9.5)') (a(i,j),j=1,ic) ; enddo
!
!  write(*,'("B[",i2,"x",i2,"]=")') im, 1
!  write(*,'(1x,10f9.5)') (xin(j),j=1,im)

    rm = sqrt(dot_product(XINX,XINX)/float(im))
    if(rm.le.rms) then
       write(*,*) '++++++ WARNING: RMS already smaller than desired '
       write(*,*) '++++++ WARNING: no correction is done            '
       rms = rm
       iter = 0
       ifail = -2
    else
       !open(61,file='fort.61')
       call htls(ax, xinx, im, ic, cinx, ny, resx, rms, 3, &
                 iter, rho, ptop, rmss, xrms, xptp, xiter, ifail)
       !close(61)
    endif

    CIN(:ic) = CINX(:ic)
    RES(:im) = RESX(:im)
    NX(NY(1:ic)) = (/ (i, i = 1, ic) /) ! NX(NY(i)) = i

!  write(*,'("X[",i2,"x",i2,"]=")') ic, 1
!  write(*,'(1x,10f9.5)') (cin(j),j=1,ic)
!
!  write(*,'("NY[",i2,"x",i2,"]=")') ic, 1
!  write(*,'(1x,10i3)') (ny(i),i=1,ic)
!
!  write(*,'("NX[",i2,"x",i2,"]=")') ic, 1
!  write(*,'(1x,10i3)') (nx(i),i=1,ic)

    return
end subroutine micit

subroutine htls(a, b, m, n, x, ipiv, r, rms, prtlev, iter, rho, ptop, &
                rmss, xrms, xptp, xiter, ifail)
     implicit none
     !*********************************************************************
     !     Subroutine HTLS to make Householder transform                  *
     !                                                                    *
     !     Authors:     many                Date:  17.09.1989             *
     !                                                                    *
     !*********************************************************************
     !     dimension of array RHO should be 3*N
     !     M  = NMTOT nr available monitors
     !     N  = NCTOT nr available independent correctors

     integer, intent(IN)  :: m, n, prtlev
     integer, intent(OUT) :: iter, ifail
     double precision :: a(m,n), b(m), x(n), r(m)
     integer :: ipiv(n)
     double precision :: rms
     double precision :: rho(3*n), ptop(n), rmss(n), xrms(n), xptp(n), xiter(n)

     integer :: i, j, j1, k, k2, k3, kpiv, kkpiv, kk, ki
     integer :: kn, kl
     double precision :: ptp, g, h, h2, h3, sig, beta, piv, pivt, kvpiv, rm, pt

     double precision, parameter :: reps7=1.d-7, zero=0d0, one=1d0
     character(len=4) :: units='mrad'

     ifail = 0
     ptp = zero

     rm = sqrt(dot_product(B(:m),B(:m))/float(m))
     pt = MAXVAL(B(:m))-MINVAL(B(:m))

     ! --- calculate first pivot
     !==========================

     RHO(:3*n) = zero

     k2 = n + 1
     piv = zero
     kpiv = 1

     do k = 1, n ! initialisation loop over correctors
        h = dot_product(A(1:m,k),A(1:m,k)) ; rho(k)  = h
        g = dot_product(A(1:m,k),B(1:m))   ; rho(k2) = g
        if (h .ne. zero) then
           pivt = g*g/h
        else
           pivt = zero
        endif
        if(pivt-piv .gt. zero) then
           piv = pivt
           kpiv = k
        endif
        k2 = k2 + 1
     enddo

!  write(*,'("A[",i2,"x",i2,"]=")') m, n
!  do i=1,m ; write(*,'(1x,10f9.5)') (a(i,j),j=1,n) ; enddo
!
!  write(*,'("B[",i2,"x",i2,"]=")') m, 1
!  write(*,'(1x,10f9.5)') (b(j),j=1,m)
!
!  write(*,'("P[",i2,"x",i2,"]=")') n, 1
!  write(*,'(1x,10i3)') (ipiv(i),i=1,n)

     ! --- boucle pour chaque iteration
     do k = 1, iter

        if (kpiv .ne. k) then
           kkpiv = kpiv
           ! --- on echange les K et KPIV si KPIV plus grand que K
           call swapreal(rho(k),rho(kpiv))
           k2 = n + k;   k3 = n + kpiv
           call swapreal(rho(k2), rho(k3))
           do i = 1, m ! swap a(i,k) and a(i,kpiv) for all monitors
              call swapreal(a(i,k),a(i,kpiv))
           enddo
          call swapint(ipiv(k),ipiv(kpiv))
        else
          kkpiv = -1
        endif
        kvpiv = piv

!  write(*,'(1x,72("-"))')
!  write(*,'(" k=",i4,", kpiv=",i4,", vpiv=",f9.5)') k, kkpiv, kvpiv
!
!  write(*,'("P[",i2,"x",i2,"]=")') n, 1
!  write(*,'(1x,10i3)') (ipiv(i),i=1,n)

        ! --- calcul de beta,sigma et uk
        sig = sqrt(dot_product(A(k:m,k),A(k:m,k)))
        h2 = a(k,k)
        sig = sign(sig,h2)
        beta = h2 + sig
        a(k,k) = beta
        beta = one/(sig*beta)

!  write(*,'(" sigma=",f9.4,",  beta=",f9.4)') sig, beta

        ! --- on garde SIGMA dans RHO(N+K)
        j = n + k
        rho(j) = -sig

        if (k .ne. n) then
           ! --- transformation de A
           do i = 1, n-k
              h = beta * dot_product(A(k:m,k),A(k:m,k+i))
              A(k:m,k+i) = A(k:m,k+i) - A(k:m,k)*h
           enddo
        endif

!  write(*,'("A[",i2,"x",i2,"]=")') m, n
!  do i=1,m ; write(*,'(1x,10f9.5)') (a(i,j),j=1,n) ; enddo

        ! --- transformation de B dans HTBL
        h3 = beta * dot_product(A(k:m,k),B(k:m))
        B(k:m) = B(k:m) - A(k:m,k)*h3

!  write(*,'("B[",i2,"x",i2,"]=")') m, 1
!  write(*,'(1x,10f9.5)') (b(j),j=1,m)

        ! --- recherche du pivot (K+1)
        !=============================
        rho(k) = sqrt(piv)

        if (k .eq. n) go to 11

        piv = zero
        kpiv = k + 1
        j1 = kpiv
        k2 = n + j1

        do j = j1, n
           h = rho(j) - a(k,j)*a(k,j)

           if(abs(h) .lt. reps7) then
              write(*,*) 'Correction process aborted'
              write(*,*) 'during ',k,'th iteration'
              write(*,*) 'Division by zero expected'
              write(*,*) 'Probably two kickers too close: ',h
              write(*,*) 'SUSPECTED KICKER: ',J,'  '
              write(*,*) 'magnet no: ', ipiv(j)
              ifail = -1
              return
           endif

           rho(j) = h
           g = rho(k2) - a(k,j)*b(k)
           rho(k2) = g
           if (h .ne. zero) then
              pivt = g*g/h
           else
              pivt = zero
           endif

           if (pivt .ge. piv) then
              kpiv=j
              piv=pivt
           endif

           k2 = k2 + 1
        enddo

!  write(*,'("RHO[",i2,"x",i2,"]=")') 1, 2*n
!  write(*,'(1x,20f9.5)') (rho(j),j=1,2*n)

        ! --- calcul des x
11      x(k) = b(k)/rho(n+k)

        if (k .ne. 1) then
           do i = 2,k
              kk = k - i + 1
              x(kk) = b(kk)
              ki = kk + 1
!  write(*,'(" x(kk)=",f9.4,", kk=",i2,", ki=",i2,", k=",i2)') x(kk), kk, ki, k
              do j = ki, k
                 x(kk) = x(kk) - a(kk,j)*x(j)
!  write(*,'(" x(kk)=",f9.4,", kk=",i2,", j=",i2,", a(kk,j)=",f9.4,", x(j)=",f9.4)') x(kk), kk, j, a(kk,j), x(j)

              enddo
              x(kk) = x(kk) / rho(n+kk)
!  write(*,'("X[",i2,"x",i2,"]=")') 1, n
!  write(*,'(1x,10f9.5)') (x(j),j=1,n)
           enddo
        endif

!  write(*,'("X[",i2,"x",i2,"]=")') 1, n
!  write(*,'(1x,10f9.5)') (x(j),j=1,n)

        ! --- save residual orbit and inverse sign of corrections (convention!)
        R(:m) = B(:m)
        X(:k) = -X(:k)

        ! --- calcul du vecteur residuel
        !===============================
        !     transform orbit R back to "normal space"
        R(1:k) = zero
        do i = 1, k
           kl = k - i + 1
           kn = k - i + 1 + n
           beta = -one / (rho(kn)*a(kl,kl))
           h = beta * dot_product(A(kl:m,kl),R(kl:m))
           R(kl:m) = R(kl:m) - A(kl:m,kl)*h
        enddo

!  write(*,'("R[",i2,"x",i2,"]=")') 1, m
!  write(*,'(1x,10f9.5)') (r(j),j=1,m)

        rmss(k) = sqrt(dot_product(R,R)/float(m))
        ptop(k) = MAXVAL(R) - MINVAL(R)

        if (k .lt. n) then
           xiter(k+1) = k
           xrms(k+1)  = rmss(k)
           xptp(k+1)  = ptop(k)
        endif

        ! --- write intermediate results to fort.61 file
        if (prtlev .ge. 2) then
          if (k .eq. 1) then
             write(61,'(/" ***********    start MICADO    ***********"/)')
             write(61,'(" iter",5X,"corrector",13X,A4,6X,"mrad",7X,"rms",11X,"ptop",/)') units
             write(61,'(4x,"0",42x,f12.8,f15.8)') rm, pt
          endif

          write(61,'(/,1x,72("-"),/)')
          write(61,'(1x,i4,3x,38(" "),1x,f12.8,f15.8)') k,rmss(k),ptop(k)

          write(61,'(I3,1X,"magnet ",I9,9X,F18.14)') (i, ipiv(i), x(i), i=1,k)

          write(61,'(/," residual orbit after iteration ",i4,":")') k
          write(61,'(1x,8f9.4)') (r(i),i=1,m)

          write(61,'(/," permutations after iteration   ",i4,": kpiv=",i4," vpiv=",f9.5)') k, kkpiv, kvpiv
          write(61,'(1x,10i6)') (ipiv(i),i=1,n)

          if (k .eq. iter) &
            write(61,'(/" ***********    end   MICADO    ***********"/)')

        endif

        if (ptop(k).le.ptp .or. rmss(k).le.rms ) then
           ! --- correction is already good enough:
           !=======================================
           ptp = ptop(k)
           rms = rmss(k)
           iter = k
           return
        endif

     enddo

     return
end subroutine htls

subroutine swapint(a, b)
  integer, intent(in out) :: a, b
  integer :: temp
  temp = a ; a = b ; b = temp
end subroutine swapint

subroutine swapreal(a, b)
  double precision, intent(in out) :: a, b
  double precision :: temp
  temp = a ; a = b ; b = temp
end subroutine swapreal

