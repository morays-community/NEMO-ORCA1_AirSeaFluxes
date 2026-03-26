MODULE stodiff
   !!======================================================================
   !!                       ***  MODULE  stodiff  ***
   !! Purpose : apply diffusion method to generate new random field
   !!           with specific horizontal correlation structure
   !!=====================================================================
   !!   sto_diff      : generate new random field with diffusion method
   !!   sto_diff_init : initialize diffusion method
   !!----------------------------------------------------------------------
   USE stoarray        ! module with stochastic arrays to update
   USE stowhite        ! uncorrelatedi normal  random number generator
   ! user supplied external resources
   USE stoexternal, only : wp, jpi, jpj, lbc_lnk, rmask2d, umask2d, vmask2d, fmask2d, &
                         & rmask3d, umask3d, vmask3d, fmask3d, grid_type, use_mask3d, narea

   IMPLICIT NONE
   PRIVATE

   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztu, ztv ! temporary workspace

   PUBLIC sto_diff, sto_diff_init

CONTAINS

   SUBROUTINE sto_diff(psto,jsto,kdim,jk)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_diff  ***
      !!
      !! ** Purpose :   generate new random field with diffusion method
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) :: psto   ! output stochastic field
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray
      INTEGER, INTENT(in) :: kdim   ! dimension of stochastic field
      INTEGER, INTENT(in) :: jk     ! index of vertical level

      INTEGER :: jpasses, jstoidx
      CHARACTER(len=1) :: sto_typ
      REAL(wp) :: sto_sgn

      ! Get index of field in 2d or 3d arrays
      jstoidx = stofields(jsto)%index

      ! Get type of grid point
      IF (kdim==2) THEN
        sto_typ = sto2d_typ( jstoidx) ; sto_sgn = sto2d_sgn( jstoidx)
      ELSEIF (kdim==3) THEN
        sto_typ = sto3d_typ( jstoidx) ; sto_sgn = sto3d_sgn( jstoidx)
      ELSE
        STOP 'Bad dimensions in sto_diff'
      ENDIF

      ! Check size of input array
      IF (SIZE(psto,1).NE.jpi) STOP 'Bad dimensions in sto_diff'
      IF (SIZE(psto,2).NE.jpj) STOP 'Bad dimensions in sto_diff'

      ! fill array with white noise
      CALL sto_white( psto2d = psto )
      ! apply mask and connections to white noise
      CALL apply_mask( psto, jsto, jk, sto_typ, sto_sgn )

      ! apply passes of the diffusion operator
      DO jpasses = 1, stofields(jsto)%diff_passes
         IF (jpasses==1) CALL noise_boundary_condition( psto )
         CALL diff_operator( psto, jsto, jk, sto_typ, sto_sgn )
         ! apply mask and connections to stochastic field
         CALL apply_mask( psto, jsto, jk, sto_typ, sto_sgn )
      END DO

      ! apply factor to restore unit standard deviation
      psto(:,:) = psto(:,:) * sto_fac(jsto)

   END SUBROUTINE sto_diff


   SUBROUTINE sto_diff_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_diff_init  ***
      !!
      !! ** Purpose :   initialize diffusion method
      !!----------------------------------------------------------------------
      INTEGER :: jsto

      ! compute factor to restore unit standard deviation
      DO jsto = 1, jpsto
         sto_fac(jsto) = diff_operator_factor( jsto ) 
      ENDDO

      ! allocate working arrays if needed
      IF (ANY(stofields(:)%diff_type==1)) THEN
         ALLOCATE(ztu(jpi,jpj))
         ALLOCATE(ztv(jpi,jpj))
      ENDIF

   END SUBROUTINE sto_diff_init


   SUBROUTINE diff_operator( psto, jsto, jk, sto_typ, sto_sgn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diff_operator  ***
      !!
      !! ** Purpose :   apply diffusion operator
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout)           ::   psto
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray
      INTEGER, INTENT(in) :: jk     ! index of vertical level
      CHARACTER(len=1), INTENT(in) :: sto_typ
      REAL(wp), INTENT(in) :: sto_sgn
      !!
      INTEGER :: diff_type
      INTEGER  :: ji, jj

      diff_type = stofields(jsto)%diff_type

      IF (diff_type==0) THEN
         ! Laplacian diffusion
         DO jj = 2, jpj-1
         DO ji = 2, jpi-1
            psto(ji,jj) = 0.5_wp * psto(ji,jj) + 0.125_wp * &
                              &  ( psto(ji-1,jj) + psto(ji+1,jj) +  &
                              &    psto(ji,jj-1) + psto(ji,jj+1) )
         END DO
         END DO
      ELSEIF (diff_type==1) THEN
         IF ( sto_typ /= 'T' ) THEN
            STOP 'Bad grid type for this diffusion option in stodiff'
         ENDIF
         IF (use_mask3d) THEN
            IF (grid_type==0) THEN  ! u points on the left, v points below r/t points (as in NEMO)
               ! Laplacian diffusion, with mask taken into account
               psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask3d(1:jpi,1:jpj,jk)
               ! 1. Gradient computation
               DO jj = 1, jpj-1
               DO ji = 1, jpi-1
                  ztu(ji,jj) = ( psto(ji+1,jj  ) - psto(ji,jj) ) * umask3d(ji,jj,jk)
                  ztv(ji,jj) = ( psto(ji  ,jj+1) - psto(ji,jj) ) * vmask3d(ji,jj,jk)
               END DO
               END DO
               ! 2. Divergence computation
               DO jj = 2, jpj-1
               DO ji = 2, jpi-1
               psto(ji,jj) = psto(ji,jj) + 0.125_wp * (  ztu(ji,jj) - ztu(ji-1,jj)   &
                  &                                    + ztv(ji,jj) - ztv(ji,jj-1)  )
               END DO
               END DO
            ELSE                    ! u points on the right, v points above r/t points (as in CROCO)
               ! Laplacian diffusion, with mask taken into account
               psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask3d(1:jpi,1:jpj,jk)
               ! 1. Gradient computation
               DO jj = 2, jpj
               DO ji = 2, jpi
                  ztu(ji,jj) = ( psto(ji,jj) - psto(ji-1,jj) ) * umask3d(ji,jj,jk)
                  ztv(ji,jj) = ( psto(ji,jj) - psto(ji,jj-1) ) * vmask3d(ji,jj,jk)
               END DO
               END DO
               ! 2. Divergence computation
               DO jj = 2, jpj-1
               DO ji = 2, jpi-1
                  psto(ji,jj) = psto(ji,jj) + 0.125_wp * (  ztu(ji+1,jj) - ztu(ji,jj)   &
                     &                                    + ztv(ji,jj+1) - ztv(ji,jj)  )
               END DO
               END DO
            ENDIF
         ELSE
            IF (grid_type==0) THEN  ! u points on the left, v points below r/t points (as in NEMO)
               ! Laplacian diffusion, with mask taken into account
               psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask2d(1:jpi,1:jpj)
               ! 1. Gradient computation
               DO jj = 1, jpj-1
               DO ji = 1, jpi-1
                  ztu(ji,jj) = ( psto(ji+1,jj  ) - psto(ji,jj) ) * umask2d(ji,jj)
                  ztv(ji,jj) = ( psto(ji  ,jj+1) - psto(ji,jj) ) * vmask2d(ji,jj)
               END DO
               END DO
               ! 2. Divergence computation
               DO jj = 2, jpj-1
               DO ji = 2, jpi-1
               psto(ji,jj) = psto(ji,jj) + 0.125_wp * (  ztu(ji,jj) - ztu(ji-1,jj)   &
                  &                                    + ztv(ji,jj) - ztv(ji,jj-1)  )
               END DO
               END DO
            ELSE                    ! u points on the right, v points above r/t points (as in CROCO)
               ! Laplacian diffusion, with mask taken into account
               psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask2d(1:jpi,1:jpj)
               ! 1. Gradient computation
               DO jj = 2, jpj
               DO ji = 2, jpi
                  ztu(ji,jj) = ( psto(ji,jj) - psto(ji-1,jj) ) * umask2d(ji,jj)
                  ztv(ji,jj) = ( psto(ji,jj) - psto(ji,jj-1) ) * vmask2d(ji,jj)
               END DO
               END DO
               ! 2. Divergence computation
               DO jj = 2, jpj-1
               DO ji = 2, jpi-1
                  psto(ji,jj) = psto(ji,jj) + 0.125_wp * (  ztu(ji+1,jj) - ztu(ji,jj)   &
                     &                                    + ztv(ji,jj+1) - ztv(ji,jj)  )
               END DO
               END DO
            ENDIF
         ENDIF
      ELSE
         STOP 'Bad diffusion operator in stodiff'
      ENDIF

   END SUBROUTINE diff_operator


   SUBROUTINE apply_mask( psto, jsto, jk, sto_typ, sto_sgn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE apply_mask  ***
      !!
      !! ** Purpose :   apply mask to input stochastic field
      !!----------------------------------------------------------------------

      REAL(wp), DIMENSION(:,:), INTENT(inout)           ::   psto
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray
      INTEGER, INTENT(in) :: jk     ! index of vertical level
      CHARACTER(len=1), INTENT(in) :: sto_typ
      REAL(wp), INTENT(in) :: sto_sgn

      ! Apply mask to stochastic field
      IF (use_mask3d) THEN
         IF (sto_typ == 'T') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask3d(1:jpi,1:jpj,jk)
         ELSEIF (sto_typ == 'U') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * umask3d(1:jpi,1:jpj,jk)
         ELSEIF (sto_typ == 'V') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * vmask3d(1:jpi,1:jpj,jk)
         ELSEIF (sto_typ == 'F') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * fmask3d(1:jpi,1:jpj,jk)
         ENDIF
      ELSE
         IF (sto_typ == 'T') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * rmask2d(1:jpi,1:jpj)
         ELSEIF (sto_typ == 'U') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * umask2d(1:jpi,1:jpj)
         ELSEIF (sto_typ == 'V') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * vmask2d(1:jpi,1:jpj)
         ELSEIF (sto_typ == 'F') THEN
            psto(1:jpi,1:jpj) = psto(1:jpi,1:jpj) * fmask2d(1:jpi,1:jpj)
         ENDIF
      ENDIF

      ! Apply connections to stochastic field
      CALL lbc_lnk( 'stodiff', psto, sto_typ, sto_sgn, ldfull = .TRUE. )

   END SUBROUTINE apply_mask


   SUBROUTINE noise_boundary_condition( psto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE noise_boundary_condition  ***
      !!
      !! ** Purpose :   apply lateral boundary conditions to input (white) noise
      !!----------------------------------------------------------------------

      REAL(wp), DIMENSION(:,:), INTENT(inout)           ::   psto

      psto(1  ,:) = 0._wp
      psto(jpi,:) = 0._wp

      psto(:,  1) = 0._wp
      psto(:,jpj) = 0._wp

   END SUBROUTINE noise_boundary_condition


   FUNCTION diff_operator_factor( jsto )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION diff_operator_factor  ***
      !!
      !! ** Purpose :   compute factor to restore standard deviation
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: jsto
      REAL(wp) :: diff_operator_factor
      !!
      INTEGER :: diff_type, npasses, jpasses, ji, jj, jflti, jfltj
      INTEGER, DIMENSION(-1:1,-1:1) :: pflt0
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pfltb
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pflta
      REAL(wp) :: ratio

      diff_type = stofields(jsto)%diff_type
      npasses   = stofields(jsto)%diff_passes

      IF ((diff_type==0).OR.(diff_type==1)) THEN
         pflt0(-1,-1) = 0 ; pflt0(-1,0) = 1 ; pflt0(-1,1) = 0
         pflt0( 0,-1) = 1 ; pflt0( 0,0) = 4 ; pflt0( 0,1) = 1
         pflt0( 1,-1) = 0 ; pflt0( 1,0) = 1 ; pflt0( 1,1) = 0
      ELSE
         STOP 'Bad diffusion operator in stodiff'
      ENDIF

      ALLOCATE(pfltb(-npasses-1:npasses+1,-npasses-1:npasses+1))
      ALLOCATE(pflta(-npasses-1:npasses+1,-npasses-1:npasses+1))

      pfltb(:,:) = 0
      pfltb(0,0) = 1
      DO jpasses = 1, npasses
        pflta(:,:) = 0
        DO jflti= -1, 1
        DO jfltj= -1, 1
          DO ji= -npasses, npasses
          DO jj= -npasses, npasses
            pflta(ji,jj) = pflta(ji,jj) + pfltb(ji+jflti,jj+jfltj) * pflt0(jflti,jfltj)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        pfltb(:,:) = pflta(:,:)
      ENDDO

      ratio = SUM(pfltb(:,:))
      ratio = ratio * ratio / SUM(pfltb(:,:)*pfltb(:,:))
      ratio = SQRT(ratio)

      DEALLOCATE(pfltb,pflta)

      diff_operator_factor = ratio

   END FUNCTION diff_operator_factor

END MODULE stodiff
