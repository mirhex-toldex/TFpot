!####################################################################### 
! Program to calculate the time propagation of the 1D Schrodinger 
! equation in a finite difference uniform grid. 
!
! Dr. Remigio Cabrera-Trujillo		June, 9th, 2005
!#######################################################################

      PROGRAM SE_1D
!
      IMPLICIT NONE
      INTEGER it, tSkip, tsk, npt                ! Time iteration
      INTEGER wrtindex, InitFlag
      INTEGER xNumPoints, xp, tp, ix             ! # point for x and t
      CHARACTER*64 InitFile
      DOUBLE PRECISION xMin, xMax, dx, x         ! x grid boundary and Delta x
      DOUBLE PRECISION tmin, tmax, dt, t         ! t grid boundary and Delta t
      DOUBLE PRECISION Mass, Vo
      DOUBLE PRECISION Ve, g0, ke, Pi, R0, width, vel, norma, norma1
      DOUBLE PRECISION KineticE, Norm, PotentialV
      DOUBLE PRECISION, ALLOCATABLE :: V(:)
      DOUBLE COMPLEX, ALLOCATABLE :: Psi(:), xDL(:), xD(:), xDU(:), &
                                     xDU2(:)
      INTEGER, ALLOCATABLE :: xIPIV(:)

!     Reads intial time, final time, steps to skip, Delta t.
      read(5,*)
      read(5,*) tmin, tmax, tSkip, dt
      read(5,*)
      read(5,*)
!     Reads space grid as xMin, xMax and number of points
      read(5,*) xMin,xMax,xNumPoints
      read(5,*)
      read(5,*)
!     For the Schroedinger with a square potential reads
!     Mass, Potential height/deep Vo, center position R0,
!     Initial gaussian wavepacket center g0, width and velocity.
!
      read(5,*) Mass, Vo, R0, g0, Width, Vel
      read(5,*)
      read(5,*)
!     In case of providing with a file with an initial wavefunction
!     0 is no, 1 is yes
      read(5,*) InitFlag
      read(5,'(A)') InitFile

      Pi = 3.14159265358979323846d0

      xp = xNumPoints+1
      dx = (xmax - xmin)/xNumPoints 
      npt =  (tmax-tmin)/dt + 1

      ALLOCATE(Psi(xp),xDL(xp), xD(xp), xDU(xp), xDU2(xp), &
               xIPIV(xp), V(xp))

      ! #########################################################
      ! * Initialize the wavefunction to a Gaussian wavepacket  *
      ! #########################################################

      IF (InitFlag .eq. 1) THEN
         open(unit=7,file=InitFile)
         DO it=1, xp
            READ(7,*) x, Psi(it)
         ENDDO
         close(7)
      ELSE
         CALL Psi1s( xp, xmin, dx, Psi, g0, width, vel)
         CALL NormPsi(xp, dx, Psi)
      ENDIF

      !#######################################################################
      ! Get the LU decomposition for the x component of the
      ! big matrix (1-v d**2) f
      !#######################################################################

      CALL LUdecomp(xp, Dx, Dt, xDL, xD, xDU, xDU2, xIPIV, mass)

      !#########################################################
      !* Loop in time starts here
      !#########################################################

      tsk=0
      wrtindex=100
      DO it=0, npt

         t=tmin + it*dt

         !#########################################################
         ! *  Error in grid to include normalization              *
         ! #########################################################

         norma1= Norma(xp, dx, Psi)

         !######################################################### 
         ! *  Now, let's calculate the Kinetic energy of this wavefunction
         ! #########################################################

         ke = KineticE(xp, dx, Psi, mass)

         !#########################################################
         ! *  Now, let's calculate the Potential energy of this wavefunction
         ! #########################################################
 
         Ve = PotentialV(xp, xmin, dx, Psi, R0, V, Vo, t)

         !#########################################################
         !*  Here we calculate e^{i(t-t0)V/2} Psi(x,t) and return 
         !*  it in the same array Psi
         !#########################################################
 
         CALL ExpVPsi(xp, Dt, Psi, V)
         
         !#########################################################
         !*  Here we apply  (1-i Dt T) Psi(x,t) and return
         !*  it in the same array Psi via the Crank-Nicholson method
         !#########################################################

         CALL xCNimplicit(xp, dx, dt, xDL, xD, xDU, xDU2, xIPIV, Psi, &
                          mass) 

         !#########################################################
         !*  Here we calculate e^{i(t-t0)V/2} Psi(x,y,z,t) and return
         !*  it in the same array Psi
         !#########################################################

         CALL ExpVPsi(xp, Dt, Psi, V)

         !#########################################################
         !*  If within the time step to write, do so! 
         !#########################################################

         IF (it - tsk .EQ. 0 ) THEN
             WRITE(6,'(F10.5,4F16.10)') t, norma1, ke, Ve, ke+Ve
             tsk = tsk + tSkip
             wrtindex = wrtindex + 1

             open(unit=wrtindex)
             DO ix=1, xp
                x = xmin + (ix-1)*dx
                WRITE(wrtindex,'(2F12.5,2F12.6)') t, x, REAL(Psi(ix)), &
                                               AIMAG(Psi(ix))
             ENDDO
             CALL FLUSH(wrtindex)
             close(wrtindex)
         ENDIF
      ENDDO

      CALL rtcoef(xmin, xp, dx, R0, Psi)

      DEALLOCATE(Psi,xDL, xD, xDU, xDU2, xIPIV, V)

      STOP
      END

!/*****************************************************************
! *  Psi0 function initialization to array Psi of Gaussian wavepacket 
! *****************************************************************/
      SUBROUTINE Psi1s(nx, xm, dx, Psi, R0, w, vel)
      IMPLICIT NONE
      INTEGER j, nx
      DOUBLE PRECISION xm, dx, x, R0, w, vel
      DOUBLE COMPLEX Psi(nx), I

      I=(0.0d0, 1.0d0)

      DO j=1, nx
         x = xm + (j-1)*dx
         Psi(j) =  EXP(-((x-R0)/w)**2 + I*vel*x)
      ENDDO

      RETURN
      END

!/*****************************************************************
! * Calculates the norm and renormalizes an array g by using the
! * grid function based on the trapezoidal scheme.
! *****************************************************************/
                                                                        
      SUBROUTINE NormPsi(nx, dx, Psi)
      IMPLICIT NONE
      INTEGER i, nx
      DOUBLE PRECISION dx, IntTrapz, norm1
      DOUBLE PRECISION f(nx)
      DOUBLE COMPLEX Psi(nx)

      DO i=1, nx
          f(i) = ABS(Psi(i))**2
      ENDDO

      Norm1 = IntTrapz(nx, dx, f)
      Norm1 = 1.0d0/DSQRT(Norm1)

      Psi = Norm1*Psi
      
      RETURN
      END

!###################################################################### 
! This subrutine calculates the LU decomposition for the Crank-Nicholson
! model of the implict time propagation method.
!
! R. Cabrera-Trujillo, January 9th, 2005
!
!###################################################################### 

      SUBROUTINE LUdecomp(n, dx, dt, DL, D, DU, DU2, IPIV, m)
      IMPLICIT NONE

      INTEGER n 
      INTEGER ii, INFO 
      INTEGER IPIV(n)
      DOUBLE PRECISION dx, dt, m
      DOUBLE COMPLEX I, nu, DL(n), D(n), DU(n), DU2(n)

      I=(0.0d0, 1.0d0)
      
      nu = -I*dt/(4.0d0*m*dx**2)
      
      DO ii = 1, n
        DL(ii) = nu
        D(ii) = 1.0d0 - 2.0d0*nu
        DU(ii) = nu
      ENDDO

      !#######################################################################
      ! ZGTTRF -compute an LU factorization of a complex tridiagonal matrix A
      ! LAPACK   using elimination with partial pivoting and row interchanges
      !#######################################################################

      CALL ZGTTRF(n, DL, D, DU, DU2, IPIV, INFO)

      IF ( INFO .NE. 0 ) THEN
           WRITE(6,*) " Error in ZGTTRF diagonalization (LU) "
           STOP
      ENDIF

      RETURN
      END


!/*****************************************************************
! * Calculates the norm i.e <g|g> of an array g by using the
! * Integrate function based on the trapezoidal scheme.
! *****************************************************************/
                                                                        
      DOUBLE PRECISION FUNCTION Norma(nx, dx, Psi)
      IMPLICIT NONE
      INTEGER i, nx
      DOUBLE PRECISION dx, IntTrapz
      DOUBLE PRECISION f(nx)
      DOUBLE COMPLEX Psi(nx)

      DO i=1, nx
          f(i) = REAL(Psi(i)*CONJG(Psi(i)))
      ENDDO

      Norma = IntTrapz(nx, dx, f)

      RETURN
      END

!/**********************************************************************
! * Calculates the Kinetic Energy 0.5*<f|-Nabla^2 |f> for the 1D array f
! *********************************************************************/

      DOUBLE PRECISION FUNCTION  KineticE(nx, dx, ff, m)
      IMPLICIT NONE
      INTEGER i, nx
      DOUBLE PRECISION dx, dx2, IntTrapz, m
      DOUBLE PRECISION, ALLOCATABLE :: tmp(:)
      DOUBLE COMPLEX D2rx, ff(nx)

      ALLOCATE(tmp(nx))

      dx2 = dx*dx

      DO i=1,  nx  
         IF ( i==1 ) THEN
             D2rx = (- 2.0d0*ff(i)+ ff(i+1))/dx2 
         ELSE IF ( i==nx ) THEN
             D2rx =(ff(i-1)- 2.0d0*ff(i))/dx2
         ELSE
             D2rx = (ff(i-1)- 2.0d0*ff(i)+ ff(i+1))/dx2
         ENDIF
         tmp(i) = REAL(ff(i)*CONJG(D2rx))
      ENDDO
      KineticE = -0.5d0*IntTrapz(nx, dx, tmp)/m

      RETURN 
      END

!/**********************************************************
!* Caluclates the averaged Potential Energy <f|V|f> using the
!* trapezoidal integration Integral.
!**********************************************************/

      DOUBLE PRECISION FUNCTION PotentialV( nx, xm, dx, func, Ro, &
                                            V, Vo, t)
      IMPLICIT NONE
      INTEGER i, nx
      DOUBLE PRECISION xm, dx, x, Vpot, IntTrapz, Ro, Vo, V(nx)
      DOUBLE PRECISION xt(nx), t
      DOUBLE COMPLEX func(nx)

      DO  i=1, nx
          x=xm + (i-1)*dx
          V(i) = Vpot(x, Ro, Vo, t)
          xt(i) = REAL(func(i)*CONJG(func(i))*V(i))
      ENDDO

      PotentialV = IntTrapz(nx, dx, xt)
      
      RETURN 
      END


      SUBROUTINE ExpVPsi(nx, dt, Psi, V)
      IMPLICIT NONE

      INTEGER k, nx
      DOUBLE PRECISION dt, V(nx)
      DOUBLE COMPLEX I, Psi(nx)
      
      I=(0.0d0,1.0d0)

      DO k=1, nx
          Psi(k) = EXP(-I*dt*V(k)/2.0d0)*Psi(k)
      ENDDO

      RETURN 
      END

!####################################################################### 
! Solves the (1-nu*dx^2)f matrix by invoking the ZGTTRS subrutine
! in the Crank-Nicholson scheme.
!
! R. Cabrera-Trujillo, January 9th, 2005
!####################################################################### 

      SUBROUTINE xCNimplicit(nx, dx, dt, xDL, xD, xDU, xDU2, xIPIV, &
                             Psi,m)
      IMPLICIT NONE

      INTEGER nx, ii, INFO, xIPIV(nx)
      DOUBLE PRECISION dx, dt,m
      DOUBLE COMPLEX nu, I, xDL(nx), xD(nx), xDU(nx), xDU2(nx), Psi(nx)
      DOUBLE COMPLEX NewPsi(nx)

      I=(0.0d0, 1.0d0)
 
      nu = -I*dt/(4.0d0*m*dx**2)

      NewPsi(1) = (1.0d0+nu)*Psi(1) - nu*Psi(2)
      DO ii = 2, nx-1
         NewPsi(ii) = (1.0d0+2.0d0*nu)*Psi(ii)- nu*(Psi(ii-1)+Psi(ii+1))
      ENDDO
      NewPsi(nx) = (1.0d0+2.0d0*nu)*Psi(nx)- nu*Psi(nx-1)

      !################################################################ 
      !  ZGTTRS - solve the systems of equations A * X = B 
      !################################################################
 
      CALL ZGTTRS('N',nx,1,xDL,xD,xDU,xDU2,xIPIV,NewPsi,nx,INFO)

      IF ( INFO .NE. 0 ) THEN
           WRITE(6,*) " xCNimplicit failed in solve ZGTTRS "
           STOP
      ENDIF

      DO ii = 1, nx
         Psi(ii) = NewPsi(ii)
      ENDDO

      RETURN
      END

!/**************************************************
! * Integration by using the trapezoidal method    *
! **************************************************/
      DOUBLE PRECISION FUNCTION IntTrapz(nx, dx, f)
      IMPLICIT NONE
      INTEGER i, nx
      DOUBLE PRECISION dx, f(nx), tmp

      tmp=0.0d0
      DO i=2, nx-1
          tmp = tmp + f(i)
      ENDDO

      IntTrapz= (tmp+0.5d0*(f(1)+f(nx)))*dx

      RETURN
      END

!/******************************************************
!* Potential for the interaction. In this case it is a
!* Constant potential 
!*******************************************************/

      DOUBLE PRECISION FUNCTION Vpot(x, Ro, Vo, t)
      IMPLICIT NONE
      DOUBLE PRECISION x, Ro, Vo, t
    
      IF (x.LT.-Ro.OR.x.GT.Ro) THEN
          Vpot = 0.0d0
      ELSE
          Vpot = Vo
      ENDIF

      RETURN
      END

!/******************************************************
!* Reflection and Transmittion coefficients at the end of
!* the dynamics: Trapezoidal integration
!*******************************************************/
!
      SUBROUTINE rtcoef(xm, nx, dx, Ro, f)
      IMPLICIT NONE
      INTEGER nx, i
      DOUBLE PRECISION x, Ro, xm, dx, rtmp, ttmp
      DOUBLE COMPLEX f(nx)

      rtmp=0.0d0
      ttmp=0.0d0
      DO i=2, nx-1
          x=xm + (i-1)*dx
          IF (x.LT.-Ro) THEN
             rtmp = rtmp + REAL(f(i)*CONJG(f(i)))
          ENDIF
          IF (x.GT.Ro) THEN
             ttmp = ttmp + REAL(f(i)*CONJG(f(i)))
          ENDIF
      ENDDO

      rtmp= (rtmp+0.5d0*f(1))*dx
      ttmp= (ttmp+0.5d0*f(nx))*dx

      WRITE(6,*) '# Reflection and transmition coeff=', rtmp, ttmp

      RETURN
      END

