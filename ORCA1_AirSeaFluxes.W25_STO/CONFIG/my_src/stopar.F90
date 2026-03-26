MODULE stopar
   !!======================================================================
   !!                       ***  MODULE  stopar  ***
   !! Stochastic parameters : definition and time stepping
   !!=====================================================================
   !! History :  3.3  ! 2011-10 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sto_par       : update the stochastic parameters
   !!   sto_par_init  : define the stochastic parameterization
   !!----------------------------------------------------------------------
   USE stoarray        ! module with stochastic arrays to update
   USE stowhite        ! uncorrelatedi normal  random number generator
   USE stodiff         ! diffusion method to generate new stochastic field
   USE stokernel       ! kernel method to generate new stochastic field
   USE stoperlin       ! Perlin noise to generate new stochastic field
   USE stomarginal     ! transformation to requested marginal distribution
   ! user supplied external resources
   USE stoexternal, only : wp, jpk, nmember, narea, mppsize, lbc_lnk

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sto_par_init    ! called by nemogcm.F90
   PUBLIC   sto_par         ! called by step.F90

   ! Parameters of the autoregressive processes (defining time correlation)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sto0d_abc ! a, b, c parameters (for 0D arrays)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sto2d_abc ! a, b, c parameters (for 2D arrays)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sto3d_abc ! a, b, c parameters (for 3D arrays)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stopar.F90 13255 2020-07-06 15:41:29Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sto_par( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par  ***
      !!
      !! ** Purpose :   update the stochastic parameters
      !!
      !! ** Method  :   model basic stochastic parameters
      !!                as a first order autoregressive process AR(1),
      !!                governed by the equation:
      !!                   X(t) = a * X(t-1) + b * w + c
      !!                where the parameters a, b and c are related
      !!                to expected value, standard deviation
      !!                and time correlation (all stationary in time) by:
      !!                   E   [X(t)]        = c / ( 1 - a )
      !!                   STD [X(t)]        = b / SQRT( 1 - a * a )
      !!                   COR [X(t),X(t-k)] = a ** k
      !!                and w is a Gaussian white noise.
      !!
      !!                We generate processes with mean=0 and std=1
      !!                       => c=0 and b = SQRT( 1 - a * a )
      !!
      !!                Higher order autoregressive proces can be optionally generated
      !!                by replacing the white noise by a lower order process.
      !!
      !!                1) The statistics of the stochastic parameters (X) are assumed
      !!                constant in space (homogeneous) and time (stationary).
      !!                This could be generalized by replacing the constant
      !!                a, b, c parameters by functions of space and time.
      !!
      !!                2) The computation is performed independently for every model
      !!                grid point, which corresponds to assume that the stochastic
      !!                parameters are uncorrelated in space.
      !!                This could be generalized by including a spatial filter: Y = Filt[ X ]
      !!                (possibly non-homgeneous and non-stationary) in the computation,
      !!                or by solving an elliptic equation: L[ Y ] = X.
      !!
      !!                3) The stochastic model for the parameters could also
      !!                be generalized to depend on the current state of the ocean (not done here).
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  :: jk, jsto, jstoidx, jidx, jidx0, jidx1, nupdate, ktmod
      REAL(wp) :: rr

      !!----------------------------------------------------------------------
      !
      ! Update 2D stochastic arrays
      !
      DO jsto = 1, jpsto2d
        IF (stofields(sto2d_idx(jsto))%type_t /= 'constant') THEN
          ! Number of time steps between update of AR processes
          nupdate = stofields(sto2d_idx(jsto))%nar_update
          ! Compute position in update interval
          ktmod = MOD(kt,nupdate)
          IF (nupdate>1) THEN ! AR process not forwarded at each time step
            IF (ktmod==0) THEN
               ! Shift array in time
               DO jidx=1,jpidx2d-1
                 sto2d(:,:,jidx,jsto) = sto2d(:,:,jidx+1,jsto)
               ENDDO
               ! Forward stochastic array in time
               jidx0 = jpidx2d-1 ; jidx1 = jpidx2d
               CALL time_forward_2d( jsto, jidx0, jidx1)
            ENDIF
            ! Interpolate AR process in time
            rr = REAL(ktmod,wp) / REAL(nupdate,wp)
            jidx0 = jpidx2d-1 ; jidx1 = jpidx2d
            CALL time_interp_2d( jsto, jidx0, jidx1, rr)
          ELSE
            ! Forward stochastic array in time
            jidx0 = 1 ; jidx1 = 1
            CALL time_forward_2d( jsto, jidx0, jidx1)
            IF (jpidxlast2d>1) sto2d(:,:,jpidxlast2d,jsto) = sto2d(:,:,jidx1,jsto)
          ENDIF
          ! Transform to requested marginal distributions
          IF (jpidxsup2d>0) THEN
            jstoidx = sto2d_idx(jsto)
            CALL sto_marginal(jstoidx)
          ENDIF
        ENDIF
      END DO
      !
      ! Update 3D stochastic arrays
      !
      DO jsto = 1, jpsto3d
        IF (stofields(sto3d_idx(jsto))%type_t(1:8) /= 'constant') THEN
          ! Number of time steps between update of AR processes
          nupdate = stofields(sto3d_idx(jsto))%nar_update
          ! Compute position in update interval
          ktmod = MOD(kt,nupdate)
          IF (nupdate>1) THEN ! AR process not forwarded at each time step
            IF (ktmod==0) THEN
              ! Shift array in time
              DO jidx=1,jpidx3d-1
                sto3d(:,:,:,jidx,jsto) = sto3d(:,:,:,jidx+1,jsto)
              ENDDO
              ! Forward stochastic array in time
              jidx0 = jpidx3d-1 ; jidx1 = jpidx3d
              CALL time_forward_3d( jsto, jidx0, jidx1)
            ENDIF
            ! Interpolate AR process in time
            rr = REAL(ktmod,wp) / REAL(nupdate,wp)
            jidx0 = jpidx3d-1 ; jidx1 = jpidx3d
            CALL time_interp_3d( jsto, jidx0, jidx1, rr)
          ELSE
            ! Forward stochastic array in time
            jidx0 = 1 ; jidx1 = 1
            CALL time_forward_3d( jsto, jidx0, jidx1)
            IF (jpidxlast3d>1) sto3d(:,:,:,jpidxlast3d,jsto) = sto3d(:,:,:,jidx1,jsto)
          ENDIF
          ! Transform to requested marginal distributions
          IF (jpidxsup3d>0) THEN
            jstoidx = sto3d_idx(jsto)
            CALL sto_marginal(jstoidx)
          ENDIF
        ENDIF
      END DO
      !
      ! Update 0D stochastic numbers
      !
      DO jsto = 1, jpsto0d
        IF (stofields(sto0d_idx(jsto))%type_t(1:8) /= 'constant') THEN
          ! Number of time steps between update of AR processes
          nupdate = stofields(sto0d_idx(jsto))%nar_update
          ! Compute position in update interval
          ktmod = MOD(kt,nupdate)
          IF (nupdate>1) THEN ! AR process not forwarded at each time step
            IF (ktmod==0) THEN
              ! Shift array in time
              DO jidx=1,jpidx0d-1
                sto0d(jidx,jsto) = sto0d(jidx+1,jsto)
              ENDDO
              ! Forward stochastic array in time
              jidx0 = jpidx0d-1 ; jidx1 = jpidx0d
              CALL time_forward_0d( jsto, jidx0, jidx1)
            ENDIF
            ! Interpolate AR process in time
            rr = REAL(ktmod,wp) / REAL(nupdate,wp)
            jidx0 = jpidx0d-1 ; jidx1 = jpidx0d
            CALL time_interp_0d( jsto, jidx0, jidx1, rr)
          ELSE
            ! Forward stochastic array in time
            jidx0 = 1 ; jidx1 = 1
            CALL time_forward_0d( jsto, jidx0, jidx1)
            IF (jpidxlast0d>1) sto0d(jpidxlast0d,jsto) = sto0d(jidx1,jsto)
          ENDIF
          ! Transform to requested marginal distributions
          IF (jpidxsup0d>0) THEN
            jstoidx = sto0d_idx(jsto)
            CALL sto_marginal(jstoidx)
          ENDIF
        ENDIF
      END DO

   END SUBROUTINE sto_par


   SUBROUTINE sto_par_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par_init  ***
      !!
      !! ** Purpose :   define the stochastic parameterization
      !!----------------------------------------------------------------------
      USE storst
      INTEGER :: seedindex

      ! Initialize random number generator
      seedindex = nmember * mppsize + narea
      CALL sto_white_init(seedindex)

      ! Initialize methods to generate spatially correlated random fields
      IF (ANY(stofields(:)%type_xy=='diffusive')) CALL sto_diff_init
      IF (ANY(stofields(:)%type_xy=='kernel'))    CALL sto_kernel_init
      IF (ANY(stofields(:)%type_xy=='perlin'))    CALL sto_perlin_init

      ! Initialize transformations to requested marginal distributions
      CALL sto_marginal_init

      ! Initialize parameters of autoregressive processes
      CALL initialize_abc

      ! Initialize random fields
      IF ( ln_rststo ) THEN
         ! Restart stochastic fields (and seeds) from file
         CALL sto_rst_read
      ELSE
         ! Generate random initial condition
         CALL sto_fields_init
      ENDIF

   END SUBROUTINE sto_par_init


   SUBROUTINE sto_new_2d( psto, jsto, jk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_new_2d  ***
      !!
      !! ** Purpose :   select method to generate new 2D field
      !!
      !! The returned field must have normal marginal distribution N(0,1)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out)           :: psto
      INTEGER,                  INTENT(in)            :: jsto
      INTEGER,                  INTENT(in), OPTIONAL  :: jk
      !!
      LOGICAL :: l2d, l3d
      INTEGER :: type_xy, jstoidx, kdim, jkused
      CHARACTER(len=1) :: sto_typ
      REAL(wp) :: sto_sgn

      l2d = .NOT.PRESENT(jk)
      l3d = PRESENT(jk)

      IF (l2d) THEN
         kdim = 2 ; jkused = 1
         jstoidx = sto2d_idx(jsto)
         type_xy = sto2d_xy(jsto)
         sto_typ = sto2d_typ(jsto)
         sto_sgn = sto2d_sgn(jsto)
      ELSEIF (l3d) THEN
         kdim = 3 ; jkused = jk
         jstoidx = sto3d_idx(jsto)
         type_xy = sto3d_xy(jsto)
         sto_typ = sto3d_typ(jsto)
         sto_sgn = sto3d_sgn(jsto)
      ENDIF

      IF (type_xy==0) THEN
         ! white noise
         CALL sto_white( psto2d = psto )
         CALL lbc_lnk( 'sto_new_2d', psto, sto_typ, sto_sgn, ldfull = .TRUE. )
      ELSEIF (type_xy==1) THEN
         ! diffusion method
         CALL sto_diff( psto,  jstoidx, kdim, jkused )
      ELSEIF (type_xy==2) THEN
         ! kernel method
         CALL sto_kernel( psto, jstoidx )
      ELSEIF (type_xy==3) THEN
         ! Perlin noise
         CALL sto_perlin( psto, jstoidx )
      ELSE
         STOP 'Bad type of method to generate 2D field in sto_par_init'
      ENDIF

   END SUBROUTINE sto_new_2d


   SUBROUTINE sto_fields_init()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_fields_init  ***
      !!
      !! ** Purpose :   initialize random fields
      !!----------------------------------------------------------------------
      INTEGER :: jsto, jk, jidx1, jstoidx

      DO jsto = 1, jpsto2d
         ! Draw random numbers from N(0,1) --> w
         CALL sto_new_2d( sto2d(:,:,1,jsto), jsto )
         ! Forward process in time to fill all time slices
         DO jidx1 = 2, jpidx2d
            CALL time_forward_2d( jsto, jidx1-1, jidx1 )
         ENDDO
         ! Transform to requested marginal distributions
         IF (jpidxsup2d>0) THEN
            sto2d(:,:,jpidxlast2d,jsto) = sto2d(:,:,jpidx2d,jsto)
            jstoidx = sto2d_idx(jsto)
            CALL sto_marginal(jstoidx)
         ENDIF
      END DO
      !
      DO jsto = 1, jpsto3d
         DO jk = 1, jpk
            ! Draw random numbers from N(0,1) --> w
            CALL sto_new_2d( sto3d(:,:,jk,1,jsto), jsto, jk )
         END DO
         ! Forward process in time to fill all time slices
         DO jidx1 = 2, jpidx3d
            CALL time_forward_3d( jsto, jidx1-1, jidx1 )
         ENDDO
         ! Transform to requested marginal distributions
         IF (jpidxsup3d>0) THEN
            sto3d(:,:,:,jpidxlast3d,jsto) = sto3d(:,:,:,jpidx3d,jsto)
            jstoidx = sto3d_idx(jsto)
            CALL sto_marginal(jstoidx)
         ENDIF
      END DO
      !
      DO jsto = 1, jpsto0d
         ! Draw random numbers from N(0,1) --> w
         CALL sto_white( psto0d = sto0d(1,jsto) )
         ! Forward processes in time to fill all time slices
         DO jidx1 = 2, jpidx0d
            CALL time_forward_0d( jsto, jidx1-1, jidx1 )
         ENDDO
         ! Transform to requested marginal distributions
         IF (jpidxsup0d>0) THEN
            sto0d(jpidxlast0d,jsto) = sto0d(jpidx0d,jsto)
            jstoidx = sto0d_idx(jsto)
            CALL sto_marginal(jstoidx)
         ENDIF
      ENDDO

   END SUBROUTINE sto_fields_init


   SUBROUTINE initialize_abc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE initialize_abc  ***
      !!
      !! ** Purpose :   initialize parameters of autoregressive processes
      !!----------------------------------------------------------------------
      INTEGER  :: jsto, jord, jordmax
      REAL(wp) :: a, b, tcor
      
      ! Allocate arrays with parameters
      IF ( jpsto0d > 0 ) ALLOCATE ( sto0d_abc(jpsto0d,2) )
      IF ( jpsto2d > 0 ) ALLOCATE ( sto2d_abc(jpsto2d,2) )
      IF ( jpsto3d > 0 ) ALLOCATE ( sto3d_abc(jpsto3d,2) )

      ! For every stochastic parameter:
      ! -------------------------------
      ! - compute parameters (a, b) of the autoregressive processes
      !   from time correlation (tcor):
      !     a = EXP ( - 1 / tcor )           --> sto2d_abc(:,1)
      !     b = SQRT( 1 - a * a )            --> sto2d_abc(:,2)
      ! - for higher order processes (ARn, n>1), use approximate formula
      !   for the b parameter (valid for tcor>>1 time step)
      ! - we simplified to mean = 0 and std = 1
      DO jsto = 1, jpsto2d
         tcor = sto2d_tcor(jsto)  ! time correlation
         jord = sto2d_ord(jsto)   ! order of process

         CALL eval_ab(a,b,tcor,jord)

         jordmax = stofields(sto2d_idx(jsto))%nar_order
         IF ((jordmax>1).AND.(jord==1)) CALL tune_ab(a,b,jordmax)

         sto2d_abc(jsto,1) = a
         sto2d_abc(jsto,2) = b
      END DO
      !
      DO jsto = 1, jpsto3d
         tcor = sto3d_tcor(jsto)  ! time correlation
         jord = sto3d_ord(jsto)   ! order of process

         CALL eval_ab(a,b,tcor,jord)

         jordmax = stofields(sto3d_idx(jsto))%nar_order
         IF ((jordmax>1).AND.(jord==1)) CALL tune_ab(a,b,jordmax)

         sto3d_abc(jsto,1) = a
         sto3d_abc(jsto,2) = b
      END DO
      !
      DO jsto = 1, jpsto0d
         tcor = sto0d_tcor(jsto)  ! time correlation
         jord = sto0d_ord(jsto)   ! order of process

         CALL eval_ab(a,b,tcor,jord)

         jordmax = stofields(sto0d_idx(jsto))%nar_order
         IF ((jordmax>1).AND.(jord==1)) CALL tune_ab(a,b,jordmax)

         sto0d_abc(jsto,1) = a
         sto0d_abc(jsto,2) = b
      END DO

   END SUBROUTINE initialize_abc


   SUBROUTINE eval_ab(a,b,tcor,jord)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eval_ab  ***
      !!
      !! ** Purpose :   initialize parameters of autoregressive processes
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: a, b
      REAL(wp), INTENT(in)  :: tcor
      INTEGER, INTENT(in)   :: jord

      IF ( tcor == 0._wp ) THEN
         a = 0._wp
      ELSE
         a = EXP ( - 1._wp / tcor )
      ENDIF
      IF ( jord == 1 ) THEN
         ! Exact formula for 1st order process
         b = SQRT ( 1._wp - a * a )
      ELSE
         ! Approximate formula, valid for tcor >> 1
         b = SQRT ( REAL( jord-1 , wp ) / REAL( 2*(2*jord-3) , wp ) )
         b = b * ( 1._wp - a * a )
      ENDIF

   END SUBROUTINE eval_ab


   SUBROUTINE tune_ab(a,b,kjord)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tune_ab  ***
      !!
      !! ** Purpose :   tune b parameter of AR1 process
      !!                to have the right standard deviation of ARn process
      !!
      !! Arguments: a, b  : parameters of the AR1 process
      !!            kjord : order n or the ARn process
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) :: a, b
      INTEGER, INTENT(in) :: kjord

      INTEGER, parameter :: jptest = 100000
      INTEGER :: jtest, jord
      REAL(wp) :: zran, zmean, zstd, zmisfit1, zmisfit2
      REAL(wp), DIMENSION(:), ALLOCATABLE :: zarn, za, zb

      IF (kjord<2) STOP 'Bad call of tune_ab in module stopar'

      ALLOCATE(zarn(kjord),za(kjord),zb(kjord))
      zarn(:) = 0._wp

      ! compute a and b parameters used for AR1 to ARn processes
      za(:) = a ; zb(1) = SQRT ( 1._wp - a * a )
      DO jord=2,kjord
         zb(jord) = SQRT ( REAL( jord-1 , wp ) / REAL( 2*(2*jord-3) , wp ) )
         zb(jord) = zb(jord) * ( 1._wp - a * a )
      ENDDO

      ! initialize mean and standard deviation
      zmean = 0._wp ;  zstd = 0._wp

      DO jtest=1, jptest
         CALL sto_white( psto0d = zran )
         ! update autoregressive processes
         zarn(1) = za(1) * zarn(1) + zb(1) * zran
         DO jord=2,kjord
            zarn(jord) = za(jord) * zarn(jord) + zb(jord) * zarn(jord-1)
         ENDDO
         ! compute mean and summmed square anomaly
         zmisfit1 = zarn(kjord) - zmean
         zmean = zmean + zmisfit1 / jtest
         zmisfit2 = zarn(kjord) - zmean
         zstd = zstd + zmisfit1 * zmisfit2
      ENDDO
      zstd = SQRT(zstd/(jptest-1))

      ! adjust b value of AR1 process to obtain the right std in ARn process
      b = b / zstd

      DEALLOCATE(zarn)

   END SUBROUTINE tune_ab


   SUBROUTINE time_forward_2d( jsto, jidx0, jidx1 )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1

      ! Store array from previous time step
      sto_tmp_2d(:,:) = sto2d(:,:,jidx0,jsto)

      IF ( sto2d_ord(jsto) == 1 ) THEN
        ! Draw new random numbers from N(0,1) --> w
        CALL sto_new_2d( sto2d(:,:,jidx1,jsto), jsto )
      ELSE
        ! Use previous process (one order lower) instead of white noise
        sto2d(:,:,jidx1,jsto) = sto2d(:,:,jidx1,jsto-1)
      ENDIF

      ! Multiply white noise (or lower order process) by b --> b * w
      sto2d(:,:,jidx1,jsto) = sto2d(:,:,jidx1,jsto) * sto2d_abc(jsto,2)
      ! Update autoregressive processes --> a * X(t-1) + b * w
      sto2d(:,:,jidx1,jsto) = sto2d(:,:,jidx1,jsto) + sto_tmp_2d(:,:) * sto2d_abc(jsto,1)

   END SUBROUTINE time_forward_2d


   SUBROUTINE time_forward_3d( jsto, jidx0, jidx1 )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1

      INTEGER :: jk

      DO jk = 1, jpk
        ! Store array from previous time step
        sto_tmp_2d(:,:) = sto3d(:,:,jk,jidx0,jsto)

        IF ( sto3d_ord(jsto) == 1 ) THEN
          ! Draw new random numbers from N(0,1) --> w
          CALL sto_new_2d( sto3d(:,:,jk,jidx1,jsto), jsto, jk )
        ELSE
          ! Use previous process (one order lower) instead of white noise
          sto3d(:,:,jk,jidx1,jsto) = sto3d(:,:,jk,jidx1,jsto-1)
        ENDIF

        ! Multiply white noise by b --> b * w
        sto3d(:,:,jk,jidx1,jsto) = sto3d(:,:,jk,jidx1,jsto) * sto3d_abc(jsto,2)
        ! Update autoregressive processes --> a * X(t-1) + b * w
        sto3d(:,:,jk,jidx1,jsto) = sto3d(:,:,jk,jidx1,jsto) + sto_tmp_2d(:,:) * sto3d_abc(jsto,1)
      END DO

   END SUBROUTINE time_forward_3d


   SUBROUTINE time_forward_0d( jsto, jidx0, jidx1 )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1

      ! Store all numbers from previous time step
      sto_tmp_0d = sto0d(jidx0,jsto)
      ! Draw new random numbers from N(0,1) --> w
      CALL sto_white( psto0d = sto0d(jidx1,jsto) )

      IF ( sto0d_ord(jsto) /= 1 ) THEN
        ! Use previous process (one order lower) instead of white noise
        sto0d(jidx1,jsto) = sto0d(jidx1,jsto-1)
      ENDIF
      ! Multiply white noise (or lower order process) by b --> b * w
      sto0d(jidx1,jsto) = sto0d(jidx1,jsto) * sto0d_abc(jsto,2)
      ! Update autoregressive processes --> a * X(t-1) + b * w
      sto0d(jidx1,jsto) = sto0d(jidx1,jsto) + sto_tmp_0d * sto0d_abc(jsto,1)

   END SUBROUTINE time_forward_0d


   SUBROUTINE time_interp_2d( jsto, jidx0, jidx1, rr )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1
      REAL(wp) :: rr

      INTEGER :: jidxset
      REAL(wp) :: a, b

      jidxset = jpidxlast2d
      a = rr ; b = 1._wp - rr

      sto2d(:,:,jidxset,jsto) = b * sto2d(:,:,jidx0,jsto) + a * sto2d(:,:,jidx1,jsto)

   END SUBROUTINE time_interp_2d


   SUBROUTINE time_interp_3d( jsto, jidx0, jidx1, rr )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1
      REAL(wp) :: rr

      INTEGER :: jidxset
      REAL(wp) :: a, b

      jidxset = jpidxlast3d
      a = rr ; b = 1._wp - rr

      sto3d(:,:,:,jidxset,jsto) = b * sto3d(:,:,:,jidx0,jsto) + a * sto3d(:,:,:,jidx1,jsto)

   END SUBROUTINE time_interp_3d


   SUBROUTINE time_interp_0d( jsto, jidx0, jidx1, rr )
      INTEGER, INTENT(in) :: jsto, jidx0, jidx1
      REAL(wp) :: rr

      INTEGER :: jidxset
      REAL(wp) :: a, b

      jidxset = jpidxlast0d
      a = rr ; b = 1._wp - rr

      sto0d(jidxset,jsto) = b * sto0d(jidx0,jsto) + a * sto0d(jidx1,jsto)

   END SUBROUTINE time_interp_0d

END MODULE stopar

