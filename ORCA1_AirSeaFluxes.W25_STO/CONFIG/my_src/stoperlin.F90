MODULE STOPERLIN
   !!=====================================================================
   !!                       ***  MODULE  stoperlin  ***
   !! Purpose : Apply Perlin noise
   !!=====================================================================
   !! sto_perlin        : calculate noise
   !! sto_perlin_init   : initialize the module
   !!=====================================================================
   !!
   !! PERLIN NOISE FORTRAN IMPLEMENTATION
   !! 
   !! This module provides function and suborutines
   !! to generate multiple Perlin noise fields,
   !! for use inside and outside NEMO.
   !! It is specifically conceived for stochastic physics
   !! packages, with an independent perturbation steup/field
   !! for each parameter/tendency.
   !!
   !! Andrea Storto @ CNR.it  -  Aug 2023
   !!
   !! ! Noise generation at each timestep / update
   !! CALL PNOISE_ADVANCE()
   !!
   !! PARAMETERS OF THE PERLIN NOISE TO PASS TO PNOISE_INIT
   !!
   !! N,M  : Size of the Perlin grid
   !! LP   : Logical(2) for cyclic conditions in x,y
   !! NOC  : Number of octaves
   !! NLA  : Lacunarity
   !! RPE  : Persistency
   !! RDC  : Distance to coast decay (m) if greater than 1km
   !!
   !! PERLIN_NOISE : low-level noise generation routine
   !!
   !!----------------------------------------------------------------------
   USE par_kind
   USE dom_oce
   USE lib_mpp
   USE in_out_manager
   USE iom
   USE ioipsl
   USE stoexternal
   USE stoarray
   USE stowhite        ! uncorrelatedi normal  random number generator

IMPLICIT NONE

PRIVATE

PUBLIC sto_perlin, sto_perlin_init

! Internal Compound/Derived type for Perlin Noise
TYPE pnoise_type
        INTEGER  :: sz_x, sz_y, nx, ny
        INTEGER  :: nn_octaves,nn_lacunarity
        INTEGER  :: nn_initmode
        REAL(wp) :: rn_persistency
        LOGICAL  :: ln_per(2)
        REAL(wp) :: rn_tau
        REAL(wp) :: rn_wa
        REAL(wp) :: rn_wb
        REAL(wp) :: rn_dcoast
        REAL(wp), ALLOCATABLE  :: pnout(:,:)
        REAL(wp), ALLOCATABLE  :: mask (:,:)
END TYPE
PUBLIC pnoise_type

! Compound/Derived type for Perlin Noise setup
TYPE pnoise_setup
        INTEGER  :: sz_x
        INTEGER  :: sz_y
        LOGICAL  :: ln_perx
        LOGICAL  :: ln_pery
        INTEGER  :: nn_octaves
        INTEGER  :: nn_lacunarity
        REAL(wp) :: rn_persistency
END TYPE
PUBLIC pnoise_setup

! Compound/Derived variable
TYPE(pnoise_type), ALLOCATABLE :: pnoise(:)

INTEGER, SAVE :: pnsize= 0
INTEGER, ALLOCATABLE :: pnindex(:)

CONTAINS

   SUBROUTINE sto_perlin_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_perlin_init  ***
      !!
      !! ** Purpose :   initialize Perlin noise
      !!----------------------------------------------------------------------
      INTEGER :: jsto, inum, jktot, ji
      REAL(wp), ALLOCATABLE :: zdc(:,:,:) 

      ! Compute maximum number of kernels to superpose
      ALLOCATE ( pnindex(jpsto) )
      jktot = 0
      DO jsto=1, jpsto                                   ! loop on all stochastic fields
         IF (stofields(jsto)%type_xy == 'perlin' ) THEN  ! with kernel option
            jktot = jktot + 1 
            pnindex(jsto) = jktot
         ENDIF
      ENDDO
      ALLOCATE( pnoise(jktot) )
      pnsize = jktot

      ji=0
      DO jsto=1, jpsto                                   ! loop on all stochastic fields
         IF (stofields(jsto)%type_xy == 'perlin' ) THEN  ! with kernel option
            ji=ji+1
            pnoise(ji)%nx=jpi
            pnoise(ji)%ny=jpj
            ALLOCATE ( pnoise(ji)%pnout( jpi, jpj ) )
            ALLOCATE ( pnoise(ji)%mask ( jpi, jpj ) )
            pnoise(ji)%pnout(:,:) = 0._wp
            pnoise(ji)%mask (:,:) = 0._wp
            pnoise(ji)%sz_x = stofields(jsto)%perlin_n
            pnoise(ji)%sz_y = stofields(jsto)%perlin_m
            pnoise(ji)%nn_octaves = stofields(jsto)%perlin_octaves
            pnoise(ji)%nn_lacunarity = stofields(jsto)%perlin_lacunarity
            pnoise(ji)%ln_per = stofields(jsto)%perlin_periodic
            pnoise(ji)%rn_persistency = stofields(jsto)%perlin_persistency
            pnoise(ji)%rn_dcoast = stofields(jsto)%perlin_dcoast
         ENDIF
      ENDDO

      IF( ANY(pnoise(1:pnsize)%rn_dcoast .GT. 1000._wp) ) THEN
              ALLOCATE( zdc(jpi,jpj,jpk) ); zdc=1._wp
              CALL iom_open('dist.coast', inum )
              CALL iom_get(inum,jpdom_auto,'Tcoast',zdc)
              CALL iom_close( inum )
              DO ji=1,pnsize
                   IF( pnoise(ji)%rn_dcoast .gt. 1000._wp ) THEN
                       pnoise(ji)%mask(:,:) = 1._wp - exp(-zdc(:,:,1)/pnoise(ji)%rn_dcoast )
                   ELSE
                       pnoise(ji)%mask(:,:) = tmask(:,:,1)
                   ENDIF
              ENDDO
              DEALLOCATE ( zdc )
      ELSE
              DO ji=1,pnsize
                   pnoise(ji)%mask(:,:) = tmask(:,:,1)
              ENDDO
      ENDIF

      IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sto_perlin init : initialization of Perlin noise'
            WRITE(numout,*) '~~~~~~~~~~~~~~~'
            WRITE(numout,*)
            WRITE(numout,*) ' Number of Perlin noise fields :',pnsize
            WRITE(numout,*)
      ENDIF

   END SUBROUTINE sto_perlin_init

   SUBROUTINE sto_perlin(psto, jsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_perlin  ***
      !!
      !! ** Purpose :   Calculate Perlin noise
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(:,:), INTENT(out) :: psto   ! output stochastic field
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray
      INTEGER :: kp
      !
      kp = pnindex(jsto)
      CALL pnoise_advance(kp,0)
      psto(:,:) = pnoise(kp)%pnout(:,:)
      !
   END SUBROUTINE sto_perlin

   SUBROUTINE pnoise_advance(kper, nn_write)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE pnoise_advance ***
      !!
      !! ** Purpose :   Wrapper for internal Perlin noise calculations
      !!----------------------------------------------------------------------
      !
      IMPLICIT NONE
      INTEGER, INTENT(in) :: kper
      INTEGER, OPTIONAL :: nn_write
      INTEGER :: ji, inum, jsto
      REAL(wp):: zsum, zcnt
      CHARACTER(len=99) :: cfout
      LOGICAL :: ln_iomput=.false.
   
      ji= kper
    
      CALL perlin_noise(pnoise(ji)%sz_x,pnoise(ji)%sz_y,pnoise(ji)%nx,pnoise(ji)%ny,&
         & pnoise(ji)%ln_per,pnoise(ji)%pnout,pnoise(ji)%nn_octaves,&
         & pnoise(ji)%nn_lacunarity,pnoise(ji)%rn_persistency )
    
      zsum =  SUM( pnoise(ji)%pnout, mask=( pnoise(ji)%mask(:,:) .gt. 0._wp ))
      zcnt =  REAL( COUNT( (pnoise(ji)%mask(:,:) .gt. 0._wp )), wp)

      CALL mpp_sum ( 'pnoise_advance', zsum )
      CALL mpp_sum ( 'pnoise_advance', zcnt )

      zsum = zsum / zcnt
    
      zsum =  sum( (pnoise(ji)%pnout-zsum)*(pnoise(ji)%pnout-zsum), mask=( pnoise(ji)%mask(:,:) .gt. 0._wp ))

      CALL mpp_sum ( 'pnoise_advance', zsum )

      zsum = zsum / zcnt
    
      pnoise(ji)%pnout = pnoise(ji)%mask(:,:)*pnoise(ji)%pnout
      pnoise(ji)%pnout = pnoise(ji)%pnout / SQRT( zsum )
    
      IF( PRESENT(nn_write) ) THEN
         IF(nn_write>0) ln_iomput = .TRUE.
      ENDIF
    
      IF( ln_iomput ) THEN
          WRITE(CFOUT,'(A,I2.2,A)') 'pnoise_',JI
          IF (iom_use(TRIM(cfout))) CALL iom_put(TRIM(CFOUT),PNOISE(JI)%pnout(:,:))
      ENDIF

   END SUBROUTINE pnoise_advance

!
! Internal routines for Perlin noise generation
!

SUBROUTINE PERLIN_NOISE_INT(N, M, NN, MM, LPER, PN)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, M, NN, MM
LOGICAL, INTENT(IN) :: LPER(2)
REAL(WP),   INTENT(OUT) :: PN(NN,MM)

INTEGER :: I, J, K, KK, IERR
INTEGER, PARAMETER :: KROOT=0
REAL(WP) :: AR(2,N,M), ZZ(2*N*M), DX, DY
REAL(WP) ::  XS(JPIGLO), YS(JPJGLO)

IF ( narea == 1 ) THEN
    CALL sto_white( psto1d = zz(:) )
ENDIF

CALL broadcast_array( zz(:) )

KK=0
DO K=1,M
  DO J=1,N
    DO I=1,2
       KK=KK+1
       AR(I,J,K) = ZZ(KK)
    ENDDO
  ENDDO
ENDDO

AR = AR / SQRT( SUM( AR*AR ) )

IF(LPER(1)) AR(:,N,:) = AR(:,1,:)
IF(LPER(2)) AR(:,:,M) = AR(:,:,1)

DX = REAL(N-1,WP) / REAL(JPIGLO,WP)
DY = REAL(M-1,WP) / REAL(JPJGLO,WP)

DO K=1,JPIGLO
   XS(K)=1._WP+(K-1)*DX
ENDDO
DO K=1,JPJGLO
   YS(K)=1._WP+(K-1)*DY
ENDDO

DO J=1,JPJ
   DO I=1,JPI
       PN(I,J) = FF(N,M,XS(MIG(I,NN_HLS)),YS(MJG(J,NN_HLS)),AR)
   ENDDO
ENDDO

END SUBROUTINE PERLIN_NOISE_INT

REAL FUNCTION FF(K,L,X,Y,VF)
IMPLICIT NONE

INTEGER,INTENT(IN) :: K,L
REAL(WP), INTENT(IN)   :: X,Y,VF(2,K,L)
INTEGER            :: I,J
REAL(WP), DIMENSION(2) :: UU, V1, V2, V3, V4
REAL(WP), DIMENSION(2) :: U1, U2, U3, U4
REAL(WP)           :: A1, A2, A3, A4, P, Q, B1, B2

I = FLOOR(X)
J = FLOOR(Y)

IF(I<1 .OR. J<1 .OR. I>K .OR. J>L) THEN
   IF(lwp) WRITE(numout,*) 'WRONG INDEX',I,J,K,L
   CALL ctl_stop('STOP','Wrong indexing in Perlin Noise module FF function')
ENDIF

UU = (/X, Y/)

V1(:) = VF(:,I  ,J  )
V2(:) = VF(:,I+1,J  )
V3(:) = VF(:,I  ,J+1)
V4(:) = VF(:,I+1,J+1)

U1(:) = UU(:) - (/I  , J  /)
U2(:) = UU(:) - (/I+1, J  /)
U3(:) = UU(:) - (/I  , J+1/)
U4(:) = UU(:) - (/I+1, J+1/)

A1    = SUM( V1 * U1 )
A2    = SUM( V2 * U2 )
A3    = SUM( V3 * U3 )
A4    = SUM( V4 * U4 )

P     = SS( X - I)
Q     = SS( Y - J)

B1    = (1._WP-P)*A1+P*A2
B2    = (1._WP-P)*A3+P*A4

FF    = (1._WP-Q)*B1 + Q*B2

END FUNCTION FF

REAL FUNCTION SS ( PP )
  IMPLICIT NONE
  REAL(WP) :: PP
  SS = 3*PP*PP - 2*PP*PP*PP
END FUNCTION SS

SUBROUTINE PERLIN_NOISE(KX, KY, NX, NY, LL_PER, PNOISE, NOCT, LACUNARITY, PERSISTENCE)
IMPLICIT NONE

INTEGER, INTENT(IN) :: KX, KY, NX, NY, NOCT
LOGICAL, INTENT(IN) :: LL_PER(2)
INTEGER, INTENT(IN) :: LACUNARITY
REAL(WP),INTENT(IN) :: PERSISTENCE
REAL(WP),INTENT(OUT):: PNOISE(NX, NY)
!
REAL(WP)            :: PNOIST(NX, NY)
REAL(WP)            :: AM
INTEGER             :: JK, FR

! Default to persistence=0.5, lacunarity=2

IF( NOCT .LT.1 .OR. NOCT .GT.8 )  CALL ctl_stop('STOP','Perlin noise: number of octaves not supported')
IF( LACUNARITY.LT. 1. )           CALL ctl_stop('STOP','Perlin noise: lacunarity not supported')
IF( PERSISTENCE .GT.1.)           CALL ctl_stop('STOP','Perlin noise: persistence not supported')

PNOISE = 0.

FR=1
AM=1.
DO JK=1, NOCT
       CALL PERLIN_NOISE_INT(KX*FR, KY*FR, NX, NY, LL_PER, PNOIST)
       PNOISE = PNOISE + AM*PNOIST
       FR = FR * LACUNARITY
       AM = AM * PERSISTENCE
ENDDO

END SUBROUTINE PERLIN_NOISE

END MODULE stoperlin
