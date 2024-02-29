!234567890
      PROGRAM CLASS_TRAJ_RECT

      IMPLICIT NONE
      INTEGER :: TSKIP,INITFLAG, WRITEYN
      INTEGER :: XN,YN, ZN, TPOINTS
      DOUBLE PRECISION :: TDELT, TNORM, XAVG, YAVG, ZAVG
      DOUBLE PRECISION :: XMIN,XMAX,YMIN,YMAX, ZMIN,ZMAX
      DOUBLE PRECISION :: xDelta, yDelta, zDelta
      DOUBLE PRECISION :: MASS1, MASS(2), MASS12, ZT, ZP, X0
      DOUBLE PRECISION, ALLOCATABLE :: R(:,:), P(:,:)

      INTEGER :: IX,IY, IZ, IT, K,XPNUMPOINTS,MATRIXDIM
      DOUBLE PRECISION :: NORM,TAVG,VAVG,TOTENERGY, VRAVG, ET
      DOUBLE PRECISION :: T, TMIN, TMAX, TMPX, TMPY, TMPZ, GETTDELT
      DOUBLE PRECISION :: PI,ERROR, OLDE, ETOL, TTMP, DTTMP
      DOUBLE PRECISION, ALLOCATABLE :: X(:), HXDIAG(:), HXOFFDIAG(:)
      DOUBLE PRECISION, ALLOCATABLE :: Y(:), HYDIAG(:),HYOFFDIAG(:)
      DOUBLE PRECISION, ALLOCATABLE :: Z(:), HZDIAG(:), HZOFFDIAG(:)
      DOUBLE PRECISION, ALLOCATABLE :: DENZ(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: XDL(:), YDL(:), ZDL(:)      ! ARRAYS FOR
      DOUBLE COMPLEX, ALLOCATABLE :: XD(:), YD(:), ZD(:)         ! THE LU
      DOUBLE COMPLEX, ALLOCATABLE :: XDU(:), YDU(:), ZDU(:)      ! DECOMP
      DOUBLE COMPLEX, ALLOCATABLE :: XDU2(:), YDU2(:), ZDU2(:)
      INTEGER, ALLOCATABLE :: XIPIV(:), YIPIV(:), ZIPIV(:)

      DOUBLE PRECISION, ALLOCATABLE :: V(:,:,:)
      DOUBLE PRECISION :: TDELT2,TDELT4
      DOUBLE COMPLEX, ALLOCATABLE :: PSI(:,:,:)
      DOUBLE COMPLEX :: CTDELT2

      CHARACTER*64 :: INITFILE

      ALLOCATE(R(2,3), P(2,3))

      PI = 3.14159265358979323846D0

      READ(5,*)
      READ(5,*) TMIN,TMAX,TDELT,TSKIP
      READ(5,*)
      READ(5,*)
      READ(5,*) XMIN,XMAX,XN
      READ(5,*)
      READ(5,*)
      READ(5,*) YMIN,YMAX,YN
      READ(5,*)
      READ(5,*)
      READ(5,*) ZMIN,ZMAX,ZN
      READ(5,*)
      READ(5,*)
      READ(5,*) MASS1, MASS(1), MASS(2), ZT, ZP
      READ(5,*)
      READ(5,*)
      READ(5,*) R(1,1), R(1,2), R(1,3)
      READ(5,*)
      READ(5,*)
      READ(5,*) R(2,1), R(2,2), R(2,3)
      READ(5,*)
      READ(5,*)
      READ(5,*) P(1,1), P(1,2), P(1,3)
      READ(5,*)
      READ(5,*)
      READ(5,*) P(2,1), P(2,2), P(2,3)
      READ(5,*)
      READ(5,*)
      READ(5,*) INITFLAG, WRITEYN
      READ(5,1001) INITFILE

      MATRIXDIM = XN*YN*ZN

      TDELT2 = TDELT/2.0D0            ! FACTOR OF 0.5 FOR CN SPLITTING
      CTDELT2 = (0.0D0,1.0D0)*TDELT2
      xDelta =(xMax-xMin)/(xN+1)
      yDelta =(yMax-yMin)/(yN+1)
      zDelta =(zMax-zMin)/(zN+1)

      ALLOCATE(X(0:XN+1), HXDIAG(XN),HXOFFDIAG(XN-1))
      ALLOCATE(Y(0:YN+1), HYDIAG(YN),HYOFFDIAG(YN-1))
      ALLOCATE(Z(0:ZN+1), HZDIAG(ZN),HZOFFDIAG(ZN-1))
      ALLOCATE(PSI(XN,YN,ZN), DENZ(XN, ZN),V(XN,YN,ZN))
      ALLOCATE(XDL(XN-1), XD(XN), XDU(XN-1), XDU2(XN-2), XIPIV(XN))
      ALLOCATE(YDL(YN-1), YD(YN), YDU(YN-1), YDU2(YN-2), YIPIV(YN))
      ALLOCATE(ZDL(ZN-1), ZD(ZN), ZDU(ZN-1), ZDU2(ZN-2), ZIPIV(ZN))

      CALL GRIDGEN(XN,xDelta,xMin, YN, yDelta,yMin, ZN, zDelta, zMin, X,Y,Z)

      CALL CALC1DHAMILTONIAN(MASS1,XN,X,HXDIAG,HXOFFDIAG, xDelta)
      CALL CALC1DHAMILTONIAN(MASS1,YN,Y,HYDIAG,HYOFFDIAG, yDelta)
      CALL CALC1DHAMILTONIAN(MASS1,ZN,Z,HZDIAG,HZOFFDIAG, zDelta)

      CALL LUDECOMP(XN, CTDELT2, HXDIAG,HXOFFDIAG, XDL, XD, XDU, XDU2, XIPIV)
      CALL LUDECOMP(YN, CTDELT2, HYDIAG,HYOFFDIAG, YDL, YD, YDU, YDU2, YIPIV)
      CALL LUDECOMP(ZN, CTDELT2, HZDIAG,HZOFFDIAG, ZDL, ZD, ZDU, ZDU2, ZIPIV)

      CALL INITPSI(XN,YN,ZN,X, Y, Z,xDelta, yDelta, zDelta,PSI,INITFLAG,&
                   INITFILE )

      WRITE(6,'(6A)') "#     TIME     NORM          TAVG",&
                     "             VAVG           ET", &
                     "           R1X	        R1Y	       R1Z", &
                     "           R2X	        R2Y	       R2Z", &
                     "           P1X            P1Y            P1Z", &
                     "           P2X            P2Y            P2Z"

      TNORM = NORM(MATRIXDIM,PSI)*xDelta*yDelta*zDelta

      CALL POTENTIAL(XN,YN,ZN,X,Y,Z,V,ZT,ZP, R)
 
      CALL CALCTAVG(XN, YN, ZN,HXDIAG, HXOFFDIAG,HYDIAG,HYOFFDIAG,&
                    HZDIAG, HZOFFDIAG,PSI,TAVG)

      TAVG=TAVG*xDelta*yDelta*zDelta
 
      CALL CALCVAVG(XN, YN, ZN, V,PSI,VAVG, R, P, MASS, ZT, ZP, ET)

      VAVG=VAVG*xDelta*yDelta*zDelta

      T=TMIN
 
      WRITE(6,10) T,TNORM,TAVG,VAVG,ET, R(1,1), R(1,2), &
                  R(1,3), R(2,1), R(2,2), R(2,3), P(1,1), P(1,2), &
                  P(1,3), P(2,1), P(2,2), P(2,3)

      IF (WRITEYN.EQ.1) THEN
         DO IX = 1,XN
           DO IZ = 1,ZN
              TMPY=0.0D0
              DO IY = 1,YN
                 TMPY = TMPY + CONJG(PSI(IX,IY, IZ))*PSI(IX,IY, IZ)
              ENDDO
              DENZ(IX, IZ) = TMPY
              WRITE(K,*) X(IX), Z(IZ), DENZ(IX,IZ)
           ENDDO
         ENDDO
      ENDIF

      TPOINTS=(TMAX-TMIN)/(TDELT*TSKIP)+1

      K=0
      DO WHILE (T.LT. TMAX)

        K=K+1

        DO IT = 1,TSKIP

          T = TDELT + T

          CALL CALCR(MASS, TDELT, R, P, ZT, ZP)

          CALL EXPV(XN,YN,ZN,CTDELT2,V,PSI)

          CALL ZIMPLICIT(XN,YN,ZN,CTDELT2,ZDL, ZD, ZDU, ZDU2,ZIPIV, &
                    HZDIAG,HZOFFDIAG,PSI)

          CALL YIMPLICIT(XN,YN,ZN,CTDELT2,YDL, YD, YDU, YDU2,YIPIV, &
                    HYDIAG,HYOFFDIAG,PSI)

          CALL XIMPLICIT(XN,YN,ZN,CTDELT2, XDL, XD, XDU,XDU2,XIPIV, &
                    HXDIAG,HXOFFDIAG,PSI)

          CALL EXPV(XN,YN,ZN,CTDELT2,V,PSI)

          CALL POTENTIAL(XN,YN,ZN,X,Y,Z,V,ZT,ZP, R)

        ENDDO

        CALL CALCTAVG(XN, YN, ZN,HXDIAG, HXOFFDIAG,HYDIAG,HYOFFDIAG,&
                      HZDIAG, HZOFFDIAG,PSI,TAVG)

        TAVG=TAVG*xDelta*yDelta*zDelta

        CALL CALCVAVG(XN, YN, ZN, V,PSI,VAVG, R, P, MASS, ZT, ZP, ET)

        VAVG=VAVG*xDelta*yDelta*zDelta

        TNORM = NORM(MATRIXDIM,PSI)*xDelta*yDelta*zDelta

        WRITE(6,10) T,TNORM,TAVG,VAVG,ET, R(1,1), R(1,2), &
                    R(1,3), R(2,1), R(2,2), R(2,3), P(1,1), P(1,2), &
                    P(1,3), P(2,1), P(2,2), P(2,3)

        call FLUSH(6)

        IF (WRITEYN.EQ.1) THEN
            DO IX = 1,XN
               DO IZ = 1,ZN
                  TMPY=0.0D0
                  DO IY = 1,YN
                     TMPY = TMPY + CONJG(PSI(IX,IY, IZ))*PSI(IX,IY, IZ)
                  ENDDO
                  DENZ(IX, IZ) = TMPY*yDelta
                  WRITE(K,*) X(IX), Z(IZ), DENZ(IX,IZ)
               ENDDO
            ENDDO
        ENDIF
      ENDDO

      OPEN(UNIT=1000)
      DO IZ = 1,ZN
         TMPZ=0.0D0
         DO IX = 1,XN
            DO IY = 1,YN
               TMPZ = TMPZ + CONJG(PSI(IX,IY, IZ))*PSI(IX,IY, IZ)
            ENDDO
         ENDDO
         WRITE(1000,20) Z(IZ), TMPZ*xDelta*yDelta
      ENDDO
      CLOSE(UNIT=1000)

      DEALLOCATE(X,HXDIAG,HXOFFDIAG)
      DEALLOCATE(Y,HYDIAG,HYOFFDIAG)
      DEALLOCATE(Z,HZDIAG,HZOFFDIAG)
      DEALLOCATE(PSI,V)
      DEALLOCATE(XDL, XD, XDU, XDU2, XIPIV)
      DEALLOCATE(YDL, YD, YDU, YDU2, YIPIV)

   10 FORMAT(F10.3,F10.5,20E16.7)
   20 FORMAT(2E16.7)
   30 FORMAT(2F10.4,I6,1P,E15.7)
   40 FORMAT('INTENSITY: ',1P,E15.7)
   50 FORMAT('******************************************************')
!'
 1001 FORMAT(A64)

      STOP
      END

!     GENERATION OF THE NUMERICAL GRID
!234567890
      subroutine GridGen(xNumPoints,xDelt,xMin, yNumPoints, yDelt, yMin, &
                         zNumPoints,zDelt, zMin, x, y, z)

      IMPLICIT NONE
      integer xNumPoints,yNumPoints, zNumPoints
      double precision xDelt,yDelt,zDelt, xMin, yMin, zMin
      double precision x(0:xNumPoints+1),y(0:yNumPoints+1),z(0:zNumPoints+1)

      integer i,ix, iy, iz

!  x -----------------------------

       do i = 0, xNumPoints+1
        x(i) = i*xDelt + xMin
       enddo

!  y -----------------------------

       do i = 0, yNumPoints+1
        y(i) = i*yDelt + yMin
       enddo

!  z -----------------------------

       do i = 0, zNumPoints+1
        z(i) = i*zDelt + zMin
       enddo

      return
      end

!     Calculates the nu in FD
!234567890
      subroutine Calc1DHamiltonian(Mass,Nx,x,HDiag,HOffDiag,Delta)

      IMPLICIT NONE
      integer Nx, i
      double precision Mass,x(0:Nx+1), HDiag(Nx), HOffDiag(Nx-1)
      double precision m, Delta

      m = 0.5d0/Mass

      do i = 1,Nx-1
         HOffDiag(i) = -m/Delta**2
      enddo

      do i = 1,Nx
         HDiag(i) = 2.0d0*m/Delta**2 
      enddo

      return
      end

!###################################################################### 
! This subrutine calculates the LU decomposition for the Crank-Nicholson
! model of the implict time propagation method.
!
! R. Cabrera-Trujillo, January 9th, 2005
!
!###################################################################### 

      SUBROUTINE LUdecomp(n, tDeltPar,HDiag, HOffDiag, DL, D, DU, DU2, &
                          IPIV)
      IMPLICIT NONE

      INTEGER n, INFO, IPIV(n)
      DOUBLE COMPLEX DL(n-1), D(n), DU(n-1), DU2(n-2), tDeltPar
      DOUBLE PRECISION HDiag(n), HOffDiag(n)

      DL = tDeltPar*HOffDiag
      D = 1.0d0+tDeltPar*HDiag
      DU = DL

!#########################################################################
! ZGTTRF  - compute an LU factorization of a complex tridiagonal matrix A
!           using elimination with partial pivoting and row interchanges
!#########################################################################
      CALL ZGTTRF(n, DL, D, DU, DU2, IPIV, INFO)

      IF ( INFO .NE. 0 ) THEN
           WRITE(6,*) " Error in ZGTTRF diagonalization (LU) "
           STOP
      ENDIF

      RETURN
      END

!     Calculates the norm of a vector
!234567890
      double precision function Norm(MatrixDim,Psi)

      IMPLICIT NONE
      integer i, MatrixDim
      double complex Psi(MatrixDim)

      Norm = 0.0d0

      do i=1,MatrixDim
         norm = norm + REAL(Psi(i)*CONJG(Psi(i)))
      enddo

      return
      end

!     Normalizes a vector
!234567890
      subroutine Normalize(MatrixDim,Psi)

      IMPLICIT NONE
      integer MatrixDim
      double complex Psi(MatrixDim)

      double precision Norm,NormFac

      NormFac = Norm(MatrixDim,Psi)
      NormFac = 1.0d0/dsqrt(NormFac)

      Psi = NormFac*Psi

      return
      end

!   Calculates the potential generated by the nuclei with an electron
!234567890
      SUBROUTINE POTENTIAL(XN,YN,ZN,X,Y,Z,V,ZT, ZP, R)

      IMPLICIT NONE
      INTEGER XN,YN, ZN, IX,IY,IZ
      DOUBLE PRECISION X(0:XN+1),Y(0:YN+1),Z(0:ZN+1)
      DOUBLE PRECISION V(XN,YN,ZN), ZT, ZP 
      DOUBLE PRECISION RE, R(2,3), RR, PI

      PI = 3.14159265358979323846D0

!
      DO IX = 1,XN
         DO IY = 1,YN
            DO IZ = 1,ZN
               RE = DSQRT((R(1,1)-X(IX))**2 + (R(1,2)-Y(IY))**2 + &
                          (R(1,3)-Z(IZ))**2)
               RR = DSQRT((R(2,1)-X(IX))**2 + (R(2,2)-Y(IY))**2 + &
                          (R(2,3)-Z(IZ))**2)
               V(IX,IY,IZ) = -ZT/RE - ZP/RR
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

!     Calculates the average of the potential energy and returns 
!     the partial total Energy ET
!234567890
      SUBROUTINE CALCVAVG(NX, NY, NZ, V,PSI,VAVG, R, P, MASS, ZT, ZP,ET)

      IMPLICIT NONE
      INTEGER NX, NY, NZ
      DOUBLE COMPLEX PSI(NX, NY, NZ)
      DOUBLE PRECISION VAVG, V(NX, NY, NZ), R(2,3), P(2,3), MASS(2), &
                       ET, ZT, ZP

      INTEGER I, J, K

      VAVG = 0.0D0
      DO I = 1,NX
         DO J=1, NY
            DO K=1, NZ
               VAVG = VAVG + CONJG(PSI(I,J,K))*V(I,J,K)*PSI(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      
      VAVG=REAL(VAVG) 

      ET = ZT*ZP/DSQRT((R(1,1)-R(2,1))**2+ (R(1,2)-R(2,2))**2+ &
                       (R(1,3)-R(2,3))**2) + &
           0.50D0*(P(1,1)**2+P(1,2)**2+P(1,3)**2)/MASS(1) + &
           0.50D0*(P(2,1)**2+P(2,2)**2+P(2,3)**2)/MASS(2)

      RETURN
      END

!    Calculates the averag of the electronic Kinetic energy
!234567890
      SUBROUTINE CALCTAVG(NX,NY,NZ,HXDIAG,HXOFFDIAG,HYDIAG,HYOFFDIAG,&
                          HZDIAG,HZOFFDIAG,PSI,TAVG)

      IMPLICIT NONE
      INTEGER NX,NY,NZ
      DOUBLE PRECISION HXDIAG(NX),HXOFFDIAG(NX-1)
      DOUBLE PRECISION HYDIAG(NY),HYOFFDIAG(NY-1)
      DOUBLE PRECISION HZDIAG(NZ),HZOFFDIAG(NZ-1)
      DOUBLE COMPLEX PSI(NX,NY,NZ), TMPX(NX), TMPY(NY), TMPZ(NZ)
      DOUBLE PRECISION TAVG, X(0:NX+1),Y(0:NY+1), Z(0:NZ+1)

      INTEGER IX,IY, IZ
      DOUBLE PRECISION EX, EY, EZ, TAVG1D

      EX = 0.0D0

!  X
      DO IY = 1,NY
         DO IZ=1,NZ
            DO IX=1,NX
               TMPX(IX) = PSI(IX,IY,IZ)
            ENDDO
            EX = EX + TAVG1D(NX,HXDIAG,HXOFFDIAG,TMPX)
         ENDDO
      ENDDO

! Y
      EY = 0.0D0

      DO IX = 1,NX
         DO IZ = 1,NZ
            DO IY=1,NY
               TMPY(IY) = PSI(IX,IY,IZ)
            ENDDO
            EY = EY + TAVG1D(NY,HYDIAG,HYOFFDIAG,TMPY)
         ENDDO
      ENDDO

! Z
      EZ = 0.0D0

      DO IX = 1,NX
         DO IY = 1,NY
            DO IZ=1,NZ
               TMPZ(IZ) = PSI(IX,IY,IZ)
            ENDDO
            EZ = EZ + TAVG1D(NZ,HZDIAG,HZOFFDIAG,TMPZ)
         ENDDO
      ENDDO

      TAVG = EX + EY + EZ

      RETURN
      END

! Calculates the new position by Euler method frpm the Momentum
!234567890
      SUBROUTINE CALCR(MASS, DT, R, P, ZT, ZP)

      IMPLICIT NONE
      INTEGER :: I, J
      DOUBLE PRECISION :: R(2,3), P(2,3), ZT, ZP, MASS(2), DT

      DO I=1,2
         DO J=1,3
            R(I,J) = R(I,J) + P(I,J)*DT/MASS(I)
         ENDDO
      ENDDO

      RETURN
      END

!    Does the E^{-iV dt}Psi
!234567890
      subroutine ExpV(Nx, Ny, Nz, ctDeltPar,V,Psi)

      IMPLICIT NONE
      integer Nx, Ny, Nz, i,j,k
      double precision V(Nx, Ny, Nz)
      double complex Psi(Nx, Ny, Nz), ctDeltPar

      do i = 1,Nx
         do j=1, Ny
            do k=1, Nz
               Psi(i,j,k)=CDEXP(-ctDeltPar*V(i,j,k))*Psi(i,j,k)
            enddo
         enddo
      enddo

      return
      end


!     This is the CN for Z
!234567890
      subroutine zImplicit(Nx,Ny,Nz,tDeltPar,DL, D, DU, DU2, IPIV,&
                           HzDiag,HzOffDiag,Psi)

      IMPLICIT NONE
      integer Nx,Ny,Nz,INFO, IPIV(Nz)
      integer ix,iy,iz

      double precision HzDiag(Nz),HzOffDiag(Nz-1)
      double complex Psi(Nx,Ny,Nz), tDeltPar, tmp(Nz)
      double complex DL(Nz-1),D(Nz),DU(Nz-1),DU2(Nz-2),NewPsi(Nz)

      do ix = 1,Nx
         do iy = 1,Ny
            do iz = 1, Nz
               tmp(iz) = Psi(ix,iy,iz)
            enddo

            call RHS(Nz,tDeltPar,HzDiag,HzOffDiag,tmp, NewPsi)

            CALL ZGTTRS('N',Nz,1,DL,D,DU,DU2,IPIV,NewPsi,Nz,INFO)

            IF ( INFO .NE. 0 ) THEN
               WRITE(6,*) " zImplicit failed in solve ZGTTRS "
               STOP
            ENDIF

            DO iz = 1, Nz
               Psi(ix,iy,iz) = NewPsi(iz)
            ENDDO
         enddo
      enddo

      return
      end

!     This is the CN for y
!234567890
      subroutine yImplicit(Nx,Ny,Nz,tDeltPar,DL, D, DU, DU2, IPIV,&
                           HyDiag,HyOffDiag,Psi)

      IMPLICIT NONE
      integer Nx,Ny,Nz,INFO, IPIV(Ny)
      double precision HyDiag(Ny),HyOffDiag(Ny-1)
      double complex Psi(Nx,Ny,Nz),tDeltPar, tmp(Ny), NewPsi(Ny)

      integer ix,iy,iz
      double complex DL(Ny-1),D(Ny),DU(Ny-1), DU2(Ny-2)

      do ix = 1,Nx
         do iz = 1,Nz
            do iy = 1, Ny
               tmp(iy) = Psi(ix,iy,iz)
            enddo

            call RHS(Ny,tDeltPar,HyDiag,HyOffDiag,tmp, NewPsi)

            CALL ZGTTRS('N',Ny,1,DL,D,DU,DU2,IPIV,NewPsi,Ny,INFO)

            IF ( INFO .NE. 0 ) THEN
               WRITE(6,*) " yImplicit failed in solve ZGTTRS "
               STOP
            ENDIF

            DO iy = 1, Ny
               Psi(ix,iy,iz) = NewPsi(iy)
            ENDDO
         enddo
      enddo

      return
      end

!     This is the CN for x
!234567890
      subroutine xImplicit(Nx,Ny,Nz,tDeltPar,DL, D, DU, DU2, IPIV,&
                           HxDiag,HxOffDiag,Psi)
      IMPLICIT NONE
      integer Nx,Ny,Nz,INFO, IPIV(Nx)
      double precision HxDiag(Nx),HxOffDiag(Nx-1)
      double complex Psi(Nx,Ny, Nz),tDeltPar, tmp(Nx), NewPsi(Nx)

      integer ix,iy, iz
      double complex DL(Nx-1),D(Nx),DU(Nx-1),DU2(Nx-2)

      do iy = 1,Ny
         do iz = 1,Nz
            do ix = 1, Nx
               tmp(ix) = Psi(ix,iy,iz)
            enddo
            call RHS(Nx,tDeltPar,HxDiag,HxOffDiag,tmp,NewPsi)

            CALL ZGTTRS('N',Nx,1,DL,D,DU,DU2,IPIV,NewPsi,Nx,INFO)

            IF ( INFO .NE. 0 ) THEN
               WRITE(6,*) " xImplicit failed in solve ZGTTRS "
               STOP
            ENDIF

            DO ix = 1, Nx
               Psi(ix,iy,iz) = NewPsi(ix)
            ENDDO
         enddo
      enddo

      return
      end
!234567890
      subroutine initPsi(xN, yN, zN, x, y, z, xD, yD, zD, Psi, InitFlag, InitFile)

      IMPLICIT NONE
      integer xN, yN, zN, InitFlag
      integer ix, iy, iz, xNt, yNt, zNt
      double precision x(0:xN+1), y(0:yN+1), z(0:zN+1), Pi
      double precision dr, mp, norm, xD, yD, zD
      double complex Psi(xN, yN, zN)
      character*64 InitFile

      Pi = 3.14159265358979323846d0

      if (InitFlag .eq. 1) then
        write(6,*) "# Reading init file = ",InitFile
        open(unit=7,file=InitFile, form='UNFORMATTED',STATUS='OLD')
        read(7) xNt, yNt, zNt
        Psi= 0.0d0
        DO ix=1, xNt
           DO iy=1, yNt
              DO iz=1, zNt
                 read(7) dr
                 Psi(ix,iy,iz) = dr*(1.0d0,0.0d0)
              ENDDO
           ENDDO
        ENDDO
        close(unit=7)
      else
!       Psi = 1.0d0
        do ix = 1,xN
           do iy = 1,yN
              do iz = 1,zN
                 dr=dsqrt(x(ix)**2+y(iy)**2+z(iz)**2)
                 Psi(ix,iy,iz) = dexp(-dr)*(1.0d0,0.0d0)
              enddo
           enddo
        enddo
      endif

      dr=NORM(xN*yN*zN,Psi)

      PSI=PSI/dsqrt(dr*xD*yD*zD)

      return
      end
!234567890
      subroutine RHS(xN,tDeltPar,HDiag,HOffDiag,Psi,NewPsi)

      IMPLICIT NONE
      integer xN
      double precision HDiag(xN),HOffDiag(xN-1)
      double complex tDeltPar,Psi(xN),NewPsi(xN)

      integer ix

      NewPsi(1) = HDiag(1)*Psi(1)+HOffDiag(1)*Psi(2)
      do ix = 2,xN-1
         NewPsi(ix) = HOffDiag(ix-1)*Psi(ix-1)+HDiag(ix)*Psi(ix)+ &
                      HOffDiag(ix)*Psi(ix+1)
      enddo
      NewPsi(xN) = HOffDiag(xN-1)*Psi(xN-1)+ HDiag(xN)*Psi(xN)

      NewPsi = Psi-tDeltPar*NewPsi

      return
      end
!234567890
      double precision function TAvg1D(N,HDiag,HOffDiag,Psi)

      IMPLICIT NONE
      integer N, i
      double precision HDiag(N),HOffDiag(N-1)
      double complex Psi(N), Energy

      Energy = CONJG(Psi(1))*(HDiag(1)*Psi(1)+HOffDiag(1)*Psi(2))
      do i = 2,N-1
         Energy = Energy + CONJG(Psi(i))*(HOffDiag(i-1)*Psi(i-1)+ &
                           HDiag(i)*Psi(i)+HOffDiag(i)*Psi(i+1))
      enddo
      Energy = Energy + CONJG(Psi(N))*(HOffDiag(N-1)*Psi(N-1)+ &
                                       HDiag(N)*Psi(N))

      TAvg1D = REAL(Energy)
      
      return
      end
