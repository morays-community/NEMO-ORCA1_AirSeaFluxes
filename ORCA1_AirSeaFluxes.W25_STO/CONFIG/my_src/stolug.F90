MODULE stolug
   !!======================================================================
   !!                       ***  MODULE stolug  ***
   !!
   !! Purpose : Location Uncertainties by Grid perturbations (lug)
   !!======================================================================
   USE par_oce
   USE dom_oce
   USE domhgr
   USE oce
   USE stoarray
   USE in_out_manager ! I/O manager
   USE lib_mpp
   USE iom

   IMPLICIT NONE
   PRIVATE

   ! Publicize variables from stoarray
   PUBLIC :: ln_sto_lug

   ! Parameters of LUG stochastic scheme
   LOGICAL          :: ln_lug_brtrp = .false.  ! apply location uncertainty in barotropic advection
   LOGICAL, PUBLIC  :: ln_lug_vrtcl = .false.  ! apply location uncertainty in vertical physics
   INTEGER          :: jsto_lugi               ! index of stochastic lug parameter (i direction)
   INTEGER          :: jsto_lugj               ! index of stochastic lug parameter (j direction)
   INTEGER, PUBLIC  :: jsto_lugk               ! index of stochastic lug parameter (k direction)
   REAL(wp)         :: rn_lug_stdxy            ! standard deviation (horizontal)
   REAL(wp)         :: rn_lug_stdz             ! standard deviation (vertical)
   REAL(wp)         :: rn_lug_tcor             ! correlation timescale (in days)
   REAL(wp)         :: rn_lug_lim   = 0.05_wp  ! maximum grid deformation
   INTEGER          :: nn_lug_flt   = 0        ! number of passes of Laplacian filter
   INTEGER          :: nn_lug_ord   = 1        ! order of autoregressive processes
   LOGICAL          :: ln_conserve  = .true.   ! locally conserve volume
   LOGICAL          :: ln_bathy     = .true.   ! moderate perturbation according to bathymetry

   ! Reference metric (without the stochastic perturbations)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ref_e1t, ref_e2t
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ref_e1u, ref_e2u
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ref_e1v, ref_e2v
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ref_e1f, ref_e2f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  e12_rescale

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  sto_fac_loc

   ! Public routines
   PUBLIC sto_lug, sto_lug_init, sto_lug_bn2, sto_lug_sh2

   !! * Substitutions
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE sto_lug(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_lug  ***
      !!
      !! Output stochastic field
      !! Update grid perturbations
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER :: ji, jj, jk

      ! Location uncertainty in barotropic advection:
      ! update perturbation of horizontal metrics
      IF (ln_sto_lug .AND. ln_lug_brtrp ) THEN

         ! Interpolate barotropic velocity at F-points (provisionally in e1u, e1v)
         DO jj = 1, jpj-1
            e1u(:,jj) = 0.5_wp * ( uu_b(:,jj,Nbb) + uu_b(:,jj+1,Nbb) )
         END DO
         DO ji = 1, jpi-1
            e1v(ji,:) = 0.5_wp * ( vv_b(ji,:,Nbb) + vv_b(ji+1,:,Nbb) )
         END DO
         CALL lbc_lnk( 'stolug', e1u(:,:), 'F' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e1v(:,:), 'F' , 1._wp , ldfull = .TRUE. )

         ! Compute norm of barotropic velocity (provisionally in e2u)
         e2u(:,:) = SQRT( e1u(:,:) * e1u(:,:) + e1v(:,:) * e1v(:,:) )

         ! Compute relative perturbation to reference metrics
         ! (proportional to barotropic displacement in correlation time scale)
         e1t(:,:) = stofields(jsto_lugi)%sto2d(:,:) * e2u(:,:) * rdt * stofields(jsto_lugi)%corr_t
         e2t(:,:) = stofields(jsto_lugj)%sto2d(:,:) * e2u(:,:) * rdt * stofields(jsto_lugj)%corr_t
         e1t(:,:) = e1t(:,:) / ref_e1t(:,:)
         e2t(:,:) = e2t(:,:) / ref_e2t(:,:)

         ! Moderate perturbation according to bathymetry
         IF (ln_bathy) THEN
            e1t(:,:) = e1t(:,:) * hf_0(:,:)/gdept_1d(jpk)
            e2t(:,:) = e2t(:,:) * hf_0(:,:)/gdept_1d(jpk)
         ENDIF

         ! Set minimum an maximum grid deformation
         WHERE ( e1t(:,:) < -rn_lug_lim ) e1t(:,:) = -rn_lug_lim
         WHERE ( e1t(:,:) >  rn_lug_lim ) e1t(:,:) =  rn_lug_lim
         WHERE ( e2t(:,:) < -rn_lug_lim ) e2t(:,:) = -rn_lug_lim
         WHERE ( e2t(:,:) >  rn_lug_lim ) e2t(:,:) =  rn_lug_lim

         ! Conserve volume of mesh cells
         IF (ln_conserve) THEN
            e1t(:,:) = 1._wp + e1t(:,:)
            e2t(:,:) = 1._wp + e2t(:,:)
            e12_rescale(:,:) = SQRT( e1t(:,:) * e2t(:,:) )
            e1t(:,:) = e1t(:,:) / e12_rescale(:,:) - 1._wp
            e2t(:,:) = e2t(:,:) / e12_rescale(:,:) - 1._wp
         ENDIF

         ! Output stochastic perturbations
         CALL iom_put( 'lug_sto_x', e1t(:,:) )
         CALL iom_put( 'lug_sto_y', e2t(:,:) )

         ! Interpolate perturbations to other grids
         DO ji = 1, jpi-1
            e1u(ji,:) = 0.5_wp * ( e1t(ji,:) + e1t(ji+1,:) )
            e2u(ji,:) = 0.5_wp * ( e2t(ji,:) + e2t(ji+1,:) )
         END DO
         e1u(jpi,:) = 0.5_wp * e1t(jpi,:)
         e2u(jpi,:) = 0.5_wp * e2t(jpi,:)

         DO jj = 1, jpj-1
            e1v(:,jj) = 0.5_wp * ( e1t(:,jj) + e1t(:,jj+1) )
            e1f(:,jj) = 0.5_wp * ( e1u(:,jj) + e1u(:,jj+1) )
            e2v(:,jj) = 0.5_wp * ( e2t(:,jj) + e2t(:,jj+1) )
            e2f(:,jj) = 0.5_wp * ( e2u(:,jj) + e2u(:,jj+1) )
         END DO
         e1v(:,jpj) = 0.5_wp * e1t(:,jpj)
         e1f(:,jpj) = 0.5_wp * e1u(:,jpj)
         e2v(:,jpj) = 0.5_wp * e2t(:,jpj)
         e2f(:,jpj) = 0.5_wp * e2u(:,jpj)

         ! Lateral boundary conditions on perturbations
         CALL lbc_lnk( 'stolug', e1t(:,:), 'T' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e2t(:,:), 'T' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e1u(:,:), 'U' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e2u(:,:), 'U' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e1v(:,:), 'V' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e2v(:,:), 'V' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e1f(:,:), 'F' , 1._wp , ldfull = .TRUE. )
         CALL lbc_lnk( 'stolug', e2f(:,:), 'F' , 1._wp , ldfull = .TRUE. )

         ! Apply perturbations to reference metrics
         e1t(:,:) = ref_e1t(:,:) * ( 1._wp + e1t(:,:) )
         e2t(:,:) = ref_e2t(:,:) * ( 1._wp + e2t(:,:) )
         e1u(:,:) = ref_e1u(:,:) * ( 1._wp + e1u(:,:) )
         e2u(:,:) = ref_e2u(:,:) * ( 1._wp + e2u(:,:) )
         e1v(:,:) = ref_e1v(:,:) * ( 1._wp + e1v(:,:) )
         e2v(:,:) = ref_e2v(:,:) * ( 1._wp + e2v(:,:) )
         e1f(:,:) = ref_e1f(:,:) * ( 1._wp + e1f(:,:) )
         e2f(:,:) = ref_e2f(:,:) * ( 1._wp + e2f(:,:) )

         ! Re-compute associated horizontal metrics arrays accordingly
         r1_e1t(:,:) = 1._wp / e1t(:,:)     ;   r1_e2t (:,:) = 1._wp / e2t(:,:)
         r1_e1u(:,:) = 1._wp / e1u(:,:)     ;   r1_e2u (:,:) = 1._wp / e2u(:,:)
         r1_e1v(:,:) = 1._wp / e1v(:,:)     ;   r1_e2v (:,:) = 1._wp / e2v(:,:)
         r1_e1f(:,:) = 1._wp / e1f(:,:)     ;   r1_e2f (:,:) = 1._wp / e2f(:,:)

         e1e2t (:,:) = e1t(:,:) * e2t(:,:)  ;   r1_e1e2t(:,:) = 1._wp / e1e2t(:,:)
         e1e2f (:,:) = e1f(:,:) * e2f(:,:)  ;   r1_e1e2f(:,:) = 1._wp / e1e2f(:,:)
         e1e2u (:,:) = e1u(:,:) * e2u(:,:)  ;   r1_e1e2u(:,:) = 1._wp / e1e2u(:,:)
         e1e2v (:,:) = e1v(:,:) * e2v(:,:)  ;   r1_e1e2v(:,:) = 1._wp / e1e2v(:,:)

         e2_e1u(:,:) = e2u(:,:) / e1u(:,:)
         e1_e2v(:,:) = e1v(:,:) / e2v(:,:)

         ! Output perturbed metrics
         CALL iom_put( 'lug_sto_e1t', e1t(:,:) )
         CALL iom_put( 'lug_sto_e2t', e2t(:,:) )

      ENDIF

      ! Location uncertainty in vertical physics
      ! update perturbation factor
      IF (ln_sto_lug .AND. ln_lug_vrtcl ) THEN

         ! location uncertainty applied to metrics dz
         ! proportional to vertical velocity (ww)
         sto_fac_loc(:,:,:) = stofields(jsto_lugk)%sto3d(:,:,:) * ww(:,:,:) * rdt
         sto_fac_loc(:,:,:) = sto_fac_loc(:,:,:) * stofields(jsto_lugk)%corr_t

         ! divide by vertical mesh size to obtain relative pertrubation
         sto_fac_loc(:,:,1) = 0.
         DO_3D( 0, 0, 0, 0, 2, jpkm1 )
            sto_fac_loc(ji,jj,jk) = sto_fac_loc(ji,jj,jk) / E3w_0(ji,jj,jk)
         END_3D
         sto_fac_loc(:,:,jpk) = 0.

         ! Set minimum an maximum relative grid perturbation
         WHERE ( sto_fac_loc(:,:,:) < -rn_lug_lim ) sto_fac_loc(:,:,:) = -rn_lug_lim
         WHERE ( sto_fac_loc(:,:,:) >  rn_lug_lim ) sto_fac_loc(:,:,:) =  rn_lug_lim

         ! Compute perturbation factor
         sto_fac_loc(:,:,:) = 1._wp + sto_fac_loc(:,:,:)

         ! Output stochastic field
         CALL iom_put( 'lug_sto_z', stofields(jsto_lugk)%sto3d(:,:,:) )

         ! Output stochastic perturbation to vertical metrics
         CALL iom_put( 'lug_sto_e3t', sto_fac_loc(:,:,:) )

      ENDIF

   END SUBROUTINE sto_lug


   SUBROUTINE sto_lug_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_lug_init  ***
      !!
      !! Read namelist namsto_lug
      !! Request and parameterize stochastic arrays
      !! Allocate and initialize reference metrics arrays
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios, ierr
      NAMELIST/namsto_lug/ ln_lug_brtrp, ln_lug_vrtcl, rn_lug_stdxy, rn_lug_stdz, &
        &                  rn_lug_tcor, nn_lug_ord, nn_lug_flt, rn_lug_lim, ln_conserve

      ! Read namelist block corresponding to this stochastic scheme
      READ_NML_REF(numnam_sto,namsto_lug)
      READ_NML_CFG(numnam_sto,namsto_lug)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sto_lug_init : location uncertainties by grid perturbations'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto_lug :'
         WRITE(numout,*) '      location uncertainty in barotropic advection ln_lug_brtrp = ', ln_lug_brtrp
         WRITE(numout,*) '      location uncertainty in vertical physics     ln_lug_vrtcl = ', ln_lug_vrtcl
         WRITE(numout,*) '      std of relative uncertainty (horizontal)     rn_lug_stdxy = ', rn_lug_stdxy
         WRITE(numout,*) '      std of relative uncertainty (vertical)       rn_lug_stdz  = ', rn_lug_stdz
         WRITE(numout,*) '      correlation timescale (in days)              rn_lug_tcor  = ', rn_lug_tcor
         WRITE(numout,*) '      order of autoregressive  processes           nn_lug_ord   = ', nn_lug_ord
         WRITE(numout,*) '      passes of Laplacian filter                   nn_lug_flt   = ', nn_lug_flt
         WRITE(numout,*) '      maximum grid deformation                     rn_lug_lim   = ', rn_lug_lim
         WRITE(numout,*) '      locally conserve volume                      ln_conserve  = ', ln_conserve
      ENDIF

      ! Define characteristics of stochastic arrays
      ! needed for location uncertainty in barotropic advection
      IF( ln_sto_lug .AND. ln_lug_brtrp ) THEN
         ! Request index for each stochastic array
         CALL sto_array_request_new(jsto_lugi)
         CALL sto_array_request_new(jsto_lugj)

         ! Set features of the corresponding stochastic field
         stofields(jsto_lugi)%dim=2
         stofields(jsto_lugi)%type_grid='T'
         stofields(jsto_lugj)%dim=2
         stofields(jsto_lugj)%type_grid='T'

         stofields(jsto_lugi)%type_t='arn'
         stofields(jsto_lugj)%type_t='arn'
         stofields(jsto_lugi)%corr_t=rn_lug_tcor * 86400. / rdt
         stofields(jsto_lugj)%corr_t=rn_lug_tcor * 86400. / rdt
         stofields(jsto_lugi)%nar_order=nn_lug_ord
         stofields(jsto_lugj)%nar_order=nn_lug_ord
         stofields(jsto_lugi)%nar_update=1
         stofields(jsto_lugj)%nar_update=1

         stofields(jsto_lugi)%type_xy='diffusive'
         stofields(jsto_lugj)%type_xy='diffusive'
         stofields(jsto_lugi)%diff_passes=nn_lug_flt
         stofields(jsto_lugj)%diff_passes=nn_lug_flt
         stofields(jsto_lugi)%diff_type=1
         stofields(jsto_lugj)%diff_type=1

         stofields(jsto_lugi)%std=rn_lug_stdxy
         stofields(jsto_lugj)%std=rn_lug_stdxy
      ENDIF

      ! Define characteristics of stochastic array
      ! needed for location uncertainty in vertical mixing
      IF( ln_sto_lug .AND. ln_lug_vrtcl ) THEN
         ! Request index for stochastic array
         CALL sto_array_request_new(jsto_lugk)

         ! Set features of the corresponding stochastic field
         stofields(jsto_lugk)%dim=3

         stofields(jsto_lugk)%type_t='arn'
         stofields(jsto_lugk)%corr_t=rn_lug_tcor * 86400. / rdt
         stofields(jsto_lugk)%nar_order=nn_lug_ord
         stofields(jsto_lugk)%nar_update=1

         stofields(jsto_lugk)%type_xy='diffusive'
         stofields(jsto_lugk)%diff_passes=nn_lug_flt
         stofields(jsto_lugk)%diff_type=1

         stofields(jsto_lugk)%std=rn_lug_stdz
      ENDIF

      IF( ln_sto_lug .AND. ln_lug_brtrp ) THEN
         ! Allocate reference metrics arrays
         ALLOCATE( ref_e1t(jpi,jpj) , ref_e2t(jpi,jpj) ,  &
             &     ref_e1u(jpi,jpj) , ref_e2u(jpi,jpj) ,  &
             &     ref_e1v(jpi,jpj) , ref_e2v(jpi,jpj) ,  &
             &     ref_e1f(jpi,jpj) , ref_e2f(jpi,jpj) ,  &
             &     STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'sto_lug : unable to allocate arrays' )

         ! Initialize reference metrics arrays
         ref_e1t(:,:) = e1t(:,:)
         ref_e2t(:,:) = e2t(:,:)
         ref_e1u(:,:) = e1u(:,:)
         ref_e2u(:,:) = e2u(:,:)
         ref_e1v(:,:) = e1v(:,:)
         ref_e2v(:,:) = e2v(:,:)
         ref_e1f(:,:) = e1f(:,:)
         ref_e2f(:,:) = e2f(:,:)

         IF( ln_conserve ) THEN
            ! Allocate rescale factor
            ALLOCATE( e12_rescale(jpi,jpj), STAT=ierr )
            IF( ierr /= 0 ) CALL ctl_stop( 'STOP', 'sto_lug : unable to allocate e12_rescale' )
         ENDIF
      ENDIF

      IF( ln_sto_lug .AND. ln_lug_vrtcl ) THEN
         ! Allocate array for stochastic factor to vertical metrics
         ALLOCATE( sto_fac_loc(jpi,jpj,jpk), STAT=ierr )
         sto_fac_loc(:,:,:) = 1._wp
         IF( ierr /= 0 ) CALL ctl_stop( 'STOP', 'sto_lug : unable to allocate sto_fac_loc' )
      ENDIF

   END SUBROUTINE sto_lug_init


   SUBROUTINE sto_lug_bn2( bn2 )
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_lug_bn2  ***
      !!
      !! Location uncertainty in Brunt-Vaisala frequency
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: bn2   ! Brunt-Vaisala frequency squared [1/s^2]

      IF (ln_sto_lug .AND. ln_lug_vrtcl ) THEN
         ! location uncertainty applied to metrics dz (since bn2 in 1/dz)
         bn2(:,:,:) = bn2(:,:,:) / sto_fac_loc(:,:,:)
      ENDIF

   END SUBROUTINE sto_lug_bn2


   SUBROUTINE sto_lug_sh2( sh2 )
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_lug_sh2  ***
      !!
      !! Location uncertainty in velocity shear
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: sh2   ! Brunt-Vaisala frequency squared [1/s^2]

      IF (ln_sto_lug .AND. ln_lug_vrtcl ) THEN
         ! location uncertainty applied to metrics dz (since sh2 in 1/dz**2)
         sh2(:,:,:) = sh2(:,:,:) / ( sto_fac_loc(:,:,:) * sto_fac_loc(:,:,:) )
      ENDIF

   END SUBROUTINE sto_lug_sh2

   !!======================================================================
END MODULE stolug
