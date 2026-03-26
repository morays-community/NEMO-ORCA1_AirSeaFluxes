MODULE stoeos
   !!======================================================================
   !!                       ***  MODULE stoeos  ***
   !!
   !! Purpose : Stochastic simulation of the effect of unresolved scales (for T and S)
   !!           in the computation of density by the sea water equation of state
   !!======================================================================
   USE par_oce
   USE dom_oce         ! ocean space and time domain
   USE oce             ! ocean dynamics and active tracers
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE stoarray
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! to use sign with key_nosignedzero
   USE lib_mpp
   USE iom

   IMPLICIT NONE
   PRIVATE

   ! Publicize variables from stoarray
   PUBLIC :: ln_sto_eos

   ! Parameters of EOS stochastic scheme
   INTEGER, PUBLIC :: nn_sto_eos = 0          ! number of degrees of freedom in stochastic equation of state
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosi ! index of stochastic eos parameter (i direction)
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosj ! index of stochastic eos parameter (j direction)
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosk ! index of stochastic eos parameter (k direction)
   REAL(wp)        :: rn_eos_deter = 0.0_wp   ! amplitude of deterministic component
   REAL(wp)        :: rn_eos_stdxy = 0.0_wp   ! random walk horz. standard deviation (in grid points)
   REAL(wp)        :: rn_eos_stdz = 0.0_wp    ! random walk vert. standard deviation (in grid points)
   REAL(wp)        :: rn_eos_tcor = 0.0_wp    ! random walk correlation timescale (in days)
   REAL(wp)        :: rn_eos_lim = 3.0_wp     ! limitation factor
   INTEGER         :: nn_eos_flt = 0          ! number of passes of Laplacian filter
   INTEGER         :: nn_eos_ord = 1          ! order of autoregressive processes

   ! Additional bounding parameters (not in namelist)
   REAL(wp)        :: rn_eos_tmin = -1.8_wp   ! minimum temperature with added fluctuations
   REAL(wp)        :: rn_eos_tmax = 32.0_wp   ! maximum temperature with added fluctuations
   REAL(wp)        :: rn_eos_smin = 28.0_wp   ! minimum salinity with added fluctuations
   REAL(wp)        :: rn_eos_smax = 38.0_wp   ! maximum salinity with added fluctuations

   ! Array with random T and S fluctuations
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:,:), ALLOCATABLE :: pts_ran

   ! Temporary workspace
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztu, ztv

   ! Public routines
   PUBLIC sto_eos, sto_eos_init

   !! * Substitutions
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"

CONTAINS

   SUBROUTINE sto_eos(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_eos  ***
      !!
      !! Output stochastic field
      !! Compute stochastic fluctuations of T and S
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      ! Output one of the stochastic field
      CALL iom_put( 'eos_sto', stofields(jsto_eosi(1))%sto2d(:,:) )

      ! Compute T and S random fluctuations for current time step
      IF ( rn_eos_deter == 0.0_wp ) THEN
         CALL sto_pts( kt, nn_sto_eos, jsto_eosi, jsto_eosj, jsto_eosk,  ts(:,:,:,:,Nnn) )
      ELSE
         CALL sto_pts( kt, nn_sto_eos, jsto_eosi, jsto_eosj, jsto_eosk,  ts(:,:,:,:,Nnn) , rn_eos_deter )
      ENDIF

   END SUBROUTINE sto_eos


   SUBROUTINE sto_eos_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_eos_init  ***
      !!
      !! 
      !! Read namelist namsto_eos
      !! Request and parameterize stochastic arrays
      !! Initialize arrays fo TS fluctuations
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios, jdof
      NAMELIST/namsto_eos/ nn_sto_eos, rn_eos_deter, rn_eos_stdxy, rn_eos_stdz, &
        &                  rn_eos_tcor, nn_eos_ord, nn_eos_flt, rn_eos_lim

      ! Read namelist block corresponding to this stochastic scheme
      READ_NML_REF(numnam_sto,namsto_eos)
      READ_NML_CFG(numnam_sto,namsto_eos)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sto_eos_init : stochastic equation of state'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto_eos :'
         WRITE(numout,*) '      number of degrees of freedom             nn_sto_eos    = ', nn_sto_eos
         WRITE(numout,*) '      deterministic component (in grid points) rn_eos_deter  = ', rn_eos_deter
         WRITE(numout,*) '      random walk horz. std (in grid points)   rn_eos_stdxy  = ', rn_eos_stdxy
         WRITE(numout,*) '      random walk vert. std (in grid points)   rn_eos_stdz   = ', rn_eos_stdz
         WRITE(numout,*) '      random walk tcor (in days)               rn_eos_tcor   = ', rn_eos_tcor
         WRITE(numout,*) '      order of autoregressive  processes       nn_eos_ord    = ', nn_eos_ord
         WRITE(numout,*) '      passes of Laplacian filter               nn_eos_flt    = ', nn_eos_flt
         WRITE(numout,*) '      limitation factor                        rn_eos_lim    = ', rn_eos_lim
      ENDIF

      ! Define characteristics of each stochastic array
      IF( nn_sto_eos > 0 ) THEN
         ALLOCATE(jsto_eosi(nn_sto_eos))
         ALLOCATE(jsto_eosj(nn_sto_eos))
         ALLOCATE(jsto_eosk(nn_sto_eos))
         DO jdof = 1, nn_sto_eos
            ! Request index for a each stochastic array
            CALL sto_array_request_new(jsto_eosi(jdof))
            CALL sto_array_request_new(jsto_eosj(jdof))
            CALL sto_array_request_new(jsto_eosk(jdof))

            ! Set features of the corresponding stochastic field
            stofields(jsto_eosi(jdof))%dim=2
            stofields(jsto_eosj(jdof))%dim=2
            stofields(jsto_eosk(jdof))%dim=2

            stofields(jsto_eosi(jdof))%type_t='arn'
            stofields(jsto_eosj(jdof))%type_t='arn'
            stofields(jsto_eosk(jdof))%type_t='arn'
            stofields(jsto_eosi(jdof))%corr_t=rn_eos_tcor * 86400. / rdt
            stofields(jsto_eosj(jdof))%corr_t=rn_eos_tcor * 86400. / rdt
            stofields(jsto_eosk(jdof))%corr_t=rn_eos_tcor * 86400. / rdt
            stofields(jsto_eosi(jdof))%nar_order=nn_eos_ord
            stofields(jsto_eosj(jdof))%nar_order=nn_eos_ord
            stofields(jsto_eosk(jdof))%nar_order=nn_eos_ord
            stofields(jsto_eosi(jdof))%nar_update=1
            stofields(jsto_eosj(jdof))%nar_update=1
            stofields(jsto_eosk(jdof))%nar_update=1

            stofields(jsto_eosi(jdof))%type_xy='diffusive'
            stofields(jsto_eosj(jdof))%type_xy='diffusive'
            stofields(jsto_eosk(jdof))%type_xy='diffusive'
            stofields(jsto_eosi(jdof))%diff_passes=nn_eos_flt
            stofields(jsto_eosj(jdof))%diff_passes=nn_eos_flt
            stofields(jsto_eosk(jdof))%diff_passes=nn_eos_flt
            stofields(jsto_eosi(jdof))%diff_type=1
            stofields(jsto_eosj(jdof))%diff_type=1
            stofields(jsto_eosK(jdof))%diff_type=1

            stofields(jsto_eosi(jdof))%min=-rn_eos_lim
            stofields(jsto_eosj(jdof))%min=-rn_eos_lim
            stofields(jsto_eosk(jdof))%min=-rn_eos_lim
            stofields(jsto_eosi(jdof))%max=rn_eos_lim
            stofields(jsto_eosj(jdof))%max=rn_eos_lim
            stofields(jsto_eosk(jdof))%max=rn_eos_lim
            stofields(jsto_eosi(jdof))%std=rn_eos_stdxy
            stofields(jsto_eosj(jdof))%std=rn_eos_stdxy
            stofields(jsto_eosk(jdof))%std=rn_eos_stdz
         END DO
      ELSE
         STOP 'Bad nn_sto_eos parameter in namsto_eos'
      ENDIF

      ! Allocate array with random T and S fluctuations
      ALLOCATE(pts_ran(jpi,jpj,jpk,jpts,nn_sto_eos))

      ! Allocate temporary workspace
      IF (nn_eos_flt>0) THEN
         ALLOCATE(ztu(jpi,jpj))
         ALLOCATE(ztv(jpi,jpj))
      ENDIF

   END SUBROUTINE sto_eos_init


   SUBROUTINE sto_pts( kt, kn_sto_eos, ksto_eosi, ksto_eosj, ksto_eosk, pts, keos_deter )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_pts  ***
      !!
      !! ** Purpose :   Compute current stochastic tracer fluctuations
      !!
      !! ** Method  :   Compute tracer differences from a random walk
      !!                around every model grid point
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, kn_sto_eos
      INTEGER, INTENT(in), DIMENSION(kn_sto_eos) :: ksto_eosi, ksto_eosj, ksto_eosk
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pts   ! 1 : potential temperature  [Celsius]
      !                                                               ! 2 : salinity               [psu]
      INTEGER  ::   ji, jj, jk, jts, jdof, jpasses ! dummy loop indices
      REAL(wp), INTENT(in), OPTIONAL :: keos_deter

      INTEGER  ::   jim1, jjm1, jkm1  ! incremented indices
      INTEGER  ::   jip1, jjp1, jkp1  !     -          -
      REAL(wp), DIMENSION(2) ::   zdtsim, zdtsjm, zdtskm         ! temporary arrays
      REAL(wp), DIMENSION(2) ::   zdtsip, zdtsjp, zdtskp, zdts   !     -        -
      REAL(wp) ::   zdrhoi, zdrhoj, zdrhonorm, znbri, znbrj, znbrk  ! temporary scalars
      LOGICAL  ::   l_deter
      !!----------------------------------------------------------------------

      l_deter = PRESENT(keos_deter)

      DO jts = 1, jpts
        CALL lbc_lnk( 'stopts', pts(:,:,:,jts), 'T' , 1._wp )
      ENDDO

      DO_3D( 1, 1, 1, 1, 1, jpkm1 )
         jkm1 = MAX(jk-1,1) ; jkp1 = MIN(jk+1,jpkm1)
         jjm1 = MAX(jj-1,1) ; jjp1 = MIN(jj+1,jpj)
         jim1 = MAX(ji-1,1) ; jip1 = MIN(ji+1,jpi)
         !
         ! compute tracer gradient
         zdtsip(:) = ( pts(jip1,jj,jk,1:2) - pts(ji,jj,jk,1:2) ) * tmask(jip1,jj,jk) * tmask(ji,jj,jk)
         zdtsim(:) = ( pts(ji,jj,jk,1:2) - pts(jim1,jj,jk,1:2) ) * tmask(jim1,jj,jk) * tmask(ji,jj,jk)
         zdtsjp(:) = ( pts(ji,jjp1,jk,1:2) - pts(ji,jj,jk,1:2) ) * tmask(ji,jjp1,jk) * tmask(ji,jj,jk)
         zdtsjm(:) = ( pts(ji,jj,jk,1:2) - pts(ji,jjm1,jk,1:2) ) * tmask(ji,jjm1,jk) * tmask(ji,jj,jk)
         zdtskp(:) = ( pts(ji,jj,jkp1,1:2) - pts(ji,jj,jk,1:2) ) * tmask(ji,jj,jkp1) * tmask(ji,jj,jk)
         zdtskm(:) = ( pts(ji,jj,jk,1:2) - pts(ji,jj,jkm1,1:2) ) * tmask(ji,jj,jkm1) * tmask(ji,jj,jk)
         !
         ! number of values in each direction
         znbri = 1.0_wp / MAX( tmask(jim1,jj,jk) + tmask(jip1,jj,jk) , 1.0_wp )
         znbrj = 1.0_wp / MAX( tmask(ji,jjm1,jk) + tmask(ji,jjp1,jk) , 1.0_wp )
         znbrk = 1.0_wp / MAX( tmask(ji,jj,jkm1) + tmask(ji,jj,jkp1) , 1.0_wp )
         !
         ! compute random tracer fluctuation (zdts)
         DO jdof = 1, kn_sto_eos
            zdts(:)   = ( zdtsip(:) + zdtsim(:) ) * znbri * stofields(ksto_eosi(jdof))%sto2d(ji,jj) + &
                      & ( zdtsjp(:) + zdtsjm(:) ) * znbrj * stofields(ksto_eosj(jdof))%sto2d(ji,jj) + &
                      & ( zdtskp(:) + zdtskm(:) ) * znbrk * stofields(ksto_eosk(jdof))%sto2d(ji,jj)
            pts_ran(ji,jj,jk,:,jdof) = zdts(:) * tmask(ji,jj,jk)
         ENDDO
         !
         ! Add deterministic component
         IF (l_deter) THEN
            ! compute direction of density gradient
            zdrhoi = rab_b(ji,jj,jk,jp_tem) * ( zdtsip(jp_tem) + zdtsim(jp_tem) ) + &
                   & rab_b(ji,jj,jk,jp_sal) * ( zdtsip(jp_sal) + zdtsim(jp_sal) )
            zdrhoj = rab_b(ji,jj,jk,jp_tem) * ( zdtsjp(jp_tem) + zdtsjm(jp_tem) ) + &
                   & rab_b(ji,jj,jk,jp_sal) * ( zdtsjp(jp_sal) + zdtsjm(jp_sal) )
            zdrhonorm = SQRT(zdrhoi**2 + zdrhoj**2)
            IF (zdrhonorm.NE.0.0) THEN
               ! compute T and S variation in direction of density gradient
               zdrhoi  = zdrhoi/zdrhonorm
               zdrhoj  = zdrhoj/zdrhonorm
               zdts(:) = ( zdtsip(:) + zdtsim(:) ) * znbri * zdrhoi + &
                       & ( zdtsjp(:) + zdtsjm(:) ) * znbrj * zdrhoj
               ! compute T and S perturbation in that direction
               zdts(:) = zdts(:) * keos_deter
               DO jdof = 1, kn_sto_eos
                  pts_ran(ji,jj,jk,:,jdof) = pts_ran(ji,jj,jk,:,jdof) + zdts(:) * tmask(ji,jj,jk)
               ENDDO
            ENDIF
          ENDIF
          !
      END_3D
      
      ! Moderate perturbation in low and high latitudes
      !DO_3D( 1, 1, 1, 1, 1, jpkm1 )
      !   pts_ran(ji,jj,jk,:,:) = pts_ran(ji,jj,jk,:,:) * tmask(ji,jj,jk) &
      !                         & * ( SIN( 2.0_wp * gphit(ji,jj) * rad ) )**2
      !END_3D

      ! Bound fluctuations to avoid producing irrealistic values of T and S
      DO jdof = 1, kn_sto_eos
         DO_3D( 1, 1, 1, 1, 1, jpkm1 )
            pts_ran(ji,jj,jk,jp_sal,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_sal,jdof)) ,  &
                                          &      MAX(pts(ji,jj,jk,jp_sal)-rn_eos_smin,0._wp) )     &
                                          &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_sal,jdof))
            pts_ran(ji,jj,jk,jp_tem,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_tem,jdof)) ,  &
                                          &      MAX(pts(ji,jj,jk,jp_tem)-rn_eos_tmin,0._wp) )     &
                                          &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_tem,jdof))
            pts_ran(ji,jj,jk,jp_sal,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_sal,jdof)) ,  &
                                          &      MAX(rn_eos_smax-pts(ji,jj,jk,jp_sal),0._wp) )     &
                                          &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_sal,jdof))
            pts_ran(ji,jj,jk,jp_tem,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_tem,jdof)) ,  &
                                          &      MAX(rn_eos_tmax-pts(ji,jj,jk,jp_tem),0._wp) )     &
                                          &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_tem,jdof))
         END_3D
      END DO

      ! Smooth fluctations with passes of a Laplacian filter
      DO jdof = 1, kn_sto_eos
         DO jts = 1, jpts
            DO jpasses = 1, nn_eos_flt
               DO jk = 1, jpkm1
                  ! Gradient computation
                  DO jj = 1, jpj-1
                  DO ji = 1, jpi-1
                     ztu(ji,jj) = ( pts_ran(ji+1,jj  ,jk,jts,jdof) - pts_ran(ji,jj,jk,jts,jdof) ) * umask(ji,jj,jk)
                     ztv(ji,jj) = ( pts_ran(ji  ,jj+1,jk,jts,jdof) - pts_ran(ji,jj,jk,jts,jdof) ) * vmask(ji,jj,jk)
                  END DO
                  END DO
                  ! Divergence computation
                  DO jj = 2, jpj-1
                  DO ji = 2, jpi-1
                     pts_ran(ji,jj,jk,jts,jdof) = pts_ran(ji,jj,jk,jts,jdof) + 0.125_wp * &
                     &                            (  ztu(ji,jj) - ztu(ji-1,jj)   &
                     &                             + ztv(ji,jj) - ztv(ji,jj-1)  )
                  END DO
                  END DO
               ENDDO
               ! Mask and lateral boundary conditions
               pts_ran(1:jpi,1:jpj,1:jpk,jts,jdof) = pts_ran(1:jpi,1:jpj,1:jpk,jts,jdof) * tmask(1:jpi,1:jpj,1:jpk)
               CALL lbc_lnk( 'stoeos', pts_ran(:,:,:,jts,jdof), 'T' , 1._wp, ldfull = .TRUE. )
            ENDDO
         ENDDO
      ENDDO


      ! Output one of the stochastic T and S flucutation
      CALL iom_put( 'eos_sto_T', pts_ran(:,:,:,jp_tem,1) )
      CALL iom_put( 'eos_sto_S', pts_ran(:,:,:,jp_sal,1) )

      ! Lateral boundary conditions on pts_ran
      DO jdof = 1, kn_sto_eos
         DO jts = 1, jpts
            CALL lbc_lnk( 'stopts', pts_ran(:,:,:,jts,jdof), 'T' , 1._wp , ldfull = .TRUE. )
         END DO
      END DO

   END SUBROUTINE sto_pts

   !!======================================================================
END MODULE stoeos
