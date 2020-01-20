!******************************************************************************!
!**                                                                          **!
!**  -- MICADO --                                                            **!
!**                                                                          **!
!******************************************************************************!

subroutine micit(cin, res, nx, rms, im, ic, iter, ny, ax, cinx, &
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

    integer :: j

    integer :: i, ny(ic)
    double precision :: rms, rm
    double precision, parameter :: zero=0d0

    write (*,'(/A, I5, A/)') 'start MICADO correction with ',iter,' correctors'

    AX(:im,:ic) = A(:im,:ic)
    CINX(:ic) = zero

    NY(1:ic) = (/ (i, i = 1, ic) /) ! NY(i) = i

    XINX(:im) = XIN(:im)
    RESX(:im) = zero

  write(*,'("A[",i2,"x",i2,"]=")') im, ic
  do i=1,im ; write(*,'(1x,10f9.5)') (a(i,j),j=1,ic) ; enddo

  write(*,'("B[",i2,"x",i2,"]=")') im, 1
  write(*,'(1x,10f9.5)') (xin(j),j=1,im)

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

  write(*,'("X[",i2,"x",i2,"]=")') ic, 1
  write(*,'(1x,10f9.5)') (cin(j),j=1,ic)

  write(*,'("NY[",i2,"x",i2,"]=")') ic, 1
  write(*,'(1x,10i3)') (ny(i),i=1,ic)

  write(*,'("NX[",i2,"x",i2,"]=")') ic, 1
  write(*,'(1x,10i3)') (nx(i),i=1,ic)

    return

!******************************************************************************!
  CONTAINS

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

end subroutine micit

!******************************************************************************!
!**                                                                          **!
!**  -- SVD COND --                                                          **!
!**                                                                          **!
!******************************************************************************!

subroutine svddec(svdmat, umat, vmat, ws, wvec, sortw, &
                  sngcut, sngval, im, ic, iflag, sing)
    implicit none
    ! ****************************************************
    !                                                    *
    !    Performs SVD analysis for any response matrix   *
    !                                                    *
    !     Author: GJR, 2015-07-01                        *
    !             based on work by WFH  12.09.02         *
    !                                                    *
    ! ****************************************************
    integer :: im, ic
    double precision :: svdmat(im,ic) !, a(im,ic)
    double precision :: umat(im,ic)
    double precision :: vmat(ic,ic)
    double precision :: ws(ic), wvec(ic)
    double precision :: sngval, sngcut
    integer :: sortw(ic), iflag, sing(2,ic)

    integer :: i, j, jj, ii, errflag
    integer, parameter :: nsing=5
    double precision :: rat

    ! SVDMAT(:im,:ic) = A(:im,:ic)

    call prepsvd(im, ic, svdmat, wvec, umat, vmat, errflag, ws)
    if (errflag .ne. 0) then
      write(*,*) 'end SVD with error code: ',errflag
      iflag = -1
      return
    endif

    call rvord(wvec,sortw,ws,ic) ! useless, SV are already sorted...

  write(*,'("V[",i2,"x",i2,"]=")') ic, ic
  do i=1,ic ; write(*,'(1x,10f9.5)') (vmat(i,j),j=1,ic) ; enddo

  write(*,'("S[",i2,"x",i2,"]=")') ic, 1
  write(*,'(1x,10f9.5)') (wvec(i),i=1,ic)

  write(*,'("W[",i2,"x",i2,"]=")') ic, 1
  write(*,'(1x,10i3)') (sortw(i),i=1,ic)

  write(*,'(" nc=",i3,", sngcut=",f9.5,", sngval=",f9.5)') min(nsing,ic), sngcut, sngval

    iflag = 0
    do  ii = 1, min(nsing,ic)
       i = sortw(ii)
!  write(*,'(" wvec(i)=",f9.5,", i=",i3)') wvec(i), i
       if ( abs(wvec(i)) .lt. sngcut) then
           do  j = 1, ic-1
              do  jj = j+1, ic

!  write(*,'(" vmat(j,i)=",f9.5,", j=",i3)') vmat(j,i), j
                 if( abs(vmat(j,i)) .gt. 1.0d-4) then
                    rat = abs(vmat(j,i)) + abs(vmat(jj,i))
                    rat = rat/abs(abs(vmat(j,i)) - abs(vmat(jj,i)))

!  write(*,'(" rat=",f9.5, ", 1/rat=",f9.5)') rat, 1/rat
                    if (rat .gt. sngval) then
                       if (iflag .lt. ic) then
                          iflag = iflag + 1 ! recorded indexes are C indexes!
                          sing(1,iflag) =  j - 1
                          sing(2,iflag) = jj - 1
                       endif
                    endif

                 endif

              enddo
           enddo
        endif
     enddo

  write(*,'("C[",i2,"x",i2,"]=")') iflag, 1
  write(*,'(1x,10i3)') (sing(1,i),i=1,iflag)

!******************************************************************************!
  CONTAINS

subroutine prepsvd(im, ic, svdmat, wvec, umat, vmat, iflag, ws)
    implicit none
    ! ******************************************************
    !                                                      *
    !   Prepare SVD and correction for any response matrix *
    !   interface routine to the SVD subroutine that takes *
    !   extended matrices such that the number of lines is *
    !   always the largest of im and ic
    !                                                      *
    !     Author: GJR, 2015-07-01                          *
    !                                                      *
    ! ******************************************************
    integer :: im, ic, iflag
    double precision :: svdmat(im,ic), umat(im,ic), vmat(ic,ic)
    double precision :: ws(ic), wvec(ic)

    double precision :: vmat2(im,ic)
    double precision :: umat2(ic,ic), svdmat2(ic,ic)

    double precision, parameter :: zero=0d0

    if (im .ge. ic) then
       VMAT2 = zero ; UMAT = zero
       call svd(im, im, ic, svdmat, wvec, .true., umat, .true., vmat2, iflag, ws)
       VMAT(:ic,:ic) = VMAT2(:ic,:ic)

    elseif (im .lt. ic) then
       SVDMAT2 = zero ;
       SVDMAT2(:im,:ic) = SVDMAT(:im,:ic)
       UMAT2 = zero ;  VMAT = zero
       call svd(ic, im, ic, svdmat2, wvec, .true., umat2, .true., vmat, iflag, ws)
       UMAT(:im,:ic) = UMAT2(:im,:ic)
       SVDMAT(:im,:ic) = SVDMAT2(:im,:ic)
    endif

    return
end subroutine prepsvd

subroutine rvord(inv,outv,ws,n)
     implicit   none
     ! 2015-Apr-28  15:45:16  ghislain: analysis
     ! subroutine to sort the indexes of elements stored in input vector INV
     ! in reverse order in output vector OUTV.
     ! WS is a workspace vector, N is the dimension of the vectors.
     ! Supposition that INV contains positive numbers only! LD: ok since S_i > 0
     integer  :: n
     double precision  ::  inv(n), ws(n)
     integer  :: outv(n)

     integer  :: i, j, jmax

     WS(:n) = INV(:n)

     do j = 1,n
        jmax = 1
        do i = 1,n
           if (ws(i).gt.ws(jmax)) jmax = i
        enddo
        outv(n-j+1) = jmax
        ws(jmax) = 0.0
     enddo

     return
end subroutine rvord

subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
    implicit  none
    ! ------------------------------------------------------------------
    !     this subroutine is a translation of the algol procedure svd,
    !     NUM. MATH. 14, 403-420(1970) by golub and reinsch.
    !     handbook for auto. comp., vol 2 -linear algebra, 134-151(1971).
    !
    !     this subroutine determines the singular value decomposition
    !          t
    !     a=usv  of a real m by n rectangular matrix.  householder
    !     bidiagonalization and a variant of the qr algorithm are used.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.  note that nm must be at least
    !          as large as the maximum of m and n.
    !
    !        m is the number of rows of a (and u).
    !
    !        n is the number of columns of a (and u) and the order of v.
    !
    !        a contains the rectangular input matrix to be decomposed.
    !
    !        matu should be set to .true. if the u matrix in the
    !          decomposition is desired, and to .false. otherwise.
    !
    !        matv should be set to .true. if the v matrix in the
    !          decomposition is desired, and to .false. otherwise.
    !
    !     on output
    !
    !        a is unaltered (unless overwritten by u or v).
    !
    !        w contains the n (non-negative) singular values of a (the
    !          diagonal elements of s).  they are unordered.  if an
    !          error exit is made, the singular values should be correct
    !          for indices ierr+1,ierr+2,...,n.
    !
    !        u contains the matrix u (orthogonal column vectors) of the
    !          decomposition if matu has been set to .true.  otherwise
    !          u is used as a temporary array.  u may coincide with a.
    !          if an error exit is made, the columns of u corresponding
    !          to indices of correct singular values should be correct.
    !
    !        v contains the matrix v (orthogonal) of the decomposition if
    !          matv has been set to .true.  otherwise v is not referenced.
    !          v may also coincide with a if u is not needed.  if an error
    !          exit is made, the columns of v corresponding to indices of
    !          correct singular values should be correct.
    !
    !        ierr is set to
    !          zero       for normal return,
    !          k          if the k-th singular value has not been
    !                     determined after 30 iterations.
    !
    !        rv1 is a temporary storage array.
    !
    !     calls pythag for  dsqrt(a*a + b*b) .
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    ! ------------------------------------------------------------------
    integer :: nm, m, n
    double precision :: a(nm,n), w(n), u(nm,n), v(nm,n), rv1(n)
    logical :: matu, matv

    integer :: i, j, k, l, ii, i1, kk, k1, ll, l1, mn, its, ierr
    double precision :: c, f, g, h, s, x, y, z, tst1, tst2, scale
    double precision :: pythag

    ierr = 0
    l1 = 0

    U(1:m,1:n) = A(1:m,1:n)
    !     .......... householder reduction to bidiagonal form ..........
    g = 0.0d0
    scale = 0.0d0
    x = 0.0d0

    do i = 1, n
        l = i + 1
        rv1(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m) go to 210
        !
        do k = i, m
            scale = scale + dabs(u(k,i))
        enddo

        if (scale .eq. 0.0d0) go to 210

        U(i:m,i) = U(i:m,i) / scale
        s = dot_product(U(i:m,i),U(i:m,i))

        f = u(i,i)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,i) = f - g
        if (i .eq. n) go to 190

        do j = l, n
           s = dot_product(U(i:m,i),U(i:m,j))
           f = s / h
           U(i:m,j) = U(i:m,j) + f*U(i:m,i)
        enddo

190     continue
        U(i:m,i) = scale * U(i:m,i)

210     w(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m .or. i .eq. n) go to 290

        do k = l, n
            scale = scale + dabs(u(i,k))
        enddo

        if (scale .eq. 0.0d0) go to 290

        U(i,l:n) = U(i,l:n) / scale
        s = dot_product(U(i,l:n),U(i,l:n))

        f = u(i,l)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,l) = f - g

        RV1(l:n) = U(i,l:n) / h

        if (i .eq. m) go to 270

        do j = l, m
           s = dot_product(U(j,l:n),U(i,l:n))
           U(j,l:n) = U(j,l:n) + s * RV1(l:n)
        enddo

270     continue
        U(i,l:n) = scale * U(i,l:n)

290     x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
     enddo
    !     .......... accumulation of right-hand transformations ..........
    if (.not. matv) go to 410
    !     .......... for i=n step -1 until 1 do -- ..........
    do ii = 1, n
        i = n + 1 - ii
        if (i .eq. n) go to 390
        if (g .eq. 0.0d0) go to 360

        do j = l, n
            !     .......... double division avoids possible underflow ..........
            v(j,i) = (u(i,j) / u(i,l)) / g
        enddo

        do j = l, n
           s = dot_product(U(i,l:n), V(l:n,j))
           V(l:n,j) = V(l:n,j) + s*V(l:n,i)
        enddo

360     continue
        V(i,l:n) = 0.0d0
        V(l:n,i) = 0.d0

390     continue
        v(i,i) = 1.0d0
        g = rv1(i)
        l = i
    enddo
    !     .......... accumulation of left-hand transformations ..........
410 if (.not. matu) go to 510
    !     ..........for i=min(m,n) step -1 until 1 do -- ..........
    mn = n
    if (m .lt. n) mn = m

    do ii = 1, mn
        i = mn + 1 - ii
        l = i + 1
        g = w(i)
        if (i .eq. n) go to 430

        U(i,l:n) = 0.0d0

430     if (g .eq. 0.0d0) go to 475
        if (i .eq. mn) go to 460

        do j = l, n
           s = dot_product(U(l:m,i),U(l:m,j))
           !     .......... double division avoids possible underflow ..........
           f = (s / u(i,i)) / g

           U(i:m,j) = U(i:m,j) + f*U(i:m,i)
        enddo

460     continue
        U(i:m,i) = U(i:m,i) / g

        go to 490

475     continue
        U(i:m,i) = 0.0d0

490     u(i,i) = u(i,i) + 1.0d0
    enddo
    !     .......... diagonalization of the bidiagonal form ..........
510 tst1 = x
    !     .......... for k=n step -1 until 1 do -- ..........
    do kk = 1, n
        k1 = n - kk
        k = k1 + 1
        its = 0
        !     .......... test for splitting.
        !                for l=k step -1 until 1 do -- ..........
520     do ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
            !     .......... rv1(1) is always zero, so there is no exit
            !                through the bottom of the loop ..........
            tst2 = tst1 + dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
        enddo
        !     .......... cancellation of rv1(l) if l greater than 1 ..........
540     c = 0.0d0
        s = 1.0d0

        do i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560

            do j = 1, m
                y = u(j,l1)
                z = u(j,i)
                u(j,l1) = y * c + z * s
                u(j,i) = -y * s + z * c
            enddo

560     continue
        enddo

        !     .......... test for convergence ..........
565     z = w(k)
        if (l .eq. k) go to 650
        !     .......... shift from bottom 2 by 2 minor ..........
        if (its .eq. 30) go to 1000
        its = its + 1
        x = w(l)
        y = w(k1)
        g = rv1(k1)
        h = rv1(k)
        f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
        g = pythag(f,1.0d0)
        f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
        !     .......... next qr transformation ..........
        c = 1.0d0
        s = 1.0d0

        do i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575

            do j = 1, n
                x = v(j,i1)
                z = v(j,i)
                v(j,i1) = x * c + z * s
                v(j,i) = -x * s + z * c
            enddo

575         z = pythag(f,h)
            w(i1) = z
            !     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
580         f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600

            do j = 1, m
                y = u(j,i1)
                z = u(j,i)
                u(j,i1) = y * c + z * s
                u(j,i) = -y * s + z * c
            enddo

600     continue
        enddo

        rv1(l) = 0.0d0
        rv1(k) = f
        w(k) = x
        go to 520
        !     .......... convergence ..........
650     if (z .ge. 0.0d0) go to 700
        !     .......... w(k) is made non-negative ..........
        w(k) = -z
        if (.not. matv) go to 700

        V(1:n,k) = -V(1:n,k)

700 continue
    enddo

    go to 1001
    !     .......... set error -- no convergence to a
    !                singular value after 30 iterations ..........
1000 ierr = k
1001 return
end subroutine svd

end subroutine svddec

double precision function pythag(a,b)
    implicit  none
    !
    !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
    !
    double precision :: a, b
    double precision :: p, r, s, t, u

    p = max(dabs(a),dabs(b))
    if (p .eq. 0.0d0) go to 20

    r = (min(dabs(a),dabs(b))/p)**2

10  continue
    t = 4.0d0 + r
    if (t .eq. 4.0d0) go to 20
    s = r/t
    u = 1.0d0 + 2.0d0*s
    p = u*p
    r = (s/u)**2 * r
    go to 10

20  pythag = p

   return
end function pythag

