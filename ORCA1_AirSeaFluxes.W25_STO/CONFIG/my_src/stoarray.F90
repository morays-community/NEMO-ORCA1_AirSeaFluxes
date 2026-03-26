MODULE stoarray
   !!======================================================================
   !!                       ***  MODULE  stoarray  ***
   !! Stochastic parameters : definition and storage of stochastic fields
   !!=====================================================================
   !! History :  4.0  ! 2024-01 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sto_param_init         : read main block of parameter from namelist
   !!   sto_array_init         : initialize stochastic arrays
   !!   sto_array_request_size : request maximum number of stochastic fields
   !!   sto_array_request_new  : request new stochastic field
   !!----------------------------------------------------------------------
   USE stoexternal, only : wp, lc, jpi, jpj, jpk, jpiglo, jpjglo, numout, load_namelist_sto, &
   &                       lwm, lwp, numnam_sto_ref, numnam_sto_cfg, numond, ctl_nam, &
   &                       ln_ens_rst_in, cn_member

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sto_param_init
   PUBLIC   sto_array_init
   PUBLIC   sto_array_request_size
   PUBLIC   sto_array_request_new

   ! General parameters of stochastic modules (read from namelist)
   LOGICAL, PUBLIC, SAVE :: ln_sto_eos  = .FALSE.    ! stochastic equation of state
   LOGICAL, PUBLIC, SAVE :: ln_sto_sppt = .FALSE.    ! SPPT scheme
   LOGICAL, PUBLIC, SAVE :: ln_sto_spp  = .FALSE.    ! SPP scheme
   LOGICAL, PUBLIC, SAVE :: ln_sto_skeb = .FALSE.    ! SKEB scheme
   LOGICAL, PUBLIC, SAVE :: ln_sto_lug  = .FALSE.    ! LUG scheme

   LOGICAL, PUBLIC           :: ln_rststo = .FALSE.  ! restart stochastic parameters from restart file
   LOGICAL, PUBLIC           :: ln_rstseed = .FALSE. ! read seed of RNG from restart file
   CHARACTER(len=lc), PUBLIC :: cn_storst_in = "restart_sto"     ! suffix of sto restart name (input)
   CHARACTER(len=lc), PUBLIC :: cn_storst_out = "restart_sto"    ! suffix of sto restart name (output)
   LOGICAL, PUBLIC           :: ln_rng_test = .FALSE.  ! perform test of random number generator
   INTEGER, PUBLIC           :: numstor, numstow     ! logical unit for restart (read and write)

   ! General type with features of stochastic fields
   TYPE sto_fields
     ! primary parameters set be user (from namelist in sto<physics>_init)
     CHARACTER(len=lc) :: stoname       ! Identifier, for instance: name=<physics>_<index
     CHARACTER(len=lc) :: type_t        ! type of time structure (constant, white, arn, spectrum,...)
     CHARACTER(len=lc) :: type_xy       ! type of xy structure (constant, white, diffusive, modes, perlin,...)
     CHARACTER(len=lc) :: type_z        ! type of z structure (constant, white, modes,...)
     CHARACTER(len=lc) :: type_variate  ! type of marginal pdf(normal,lognormal,gamma,beta,...)
     CHARACTER(len=1)  :: type_grid     ! type of grid (T, U, V, F)
     REAL(wp) :: sign_grid              ! sign to use to connect field across domain
     REAL(wp) :: ave                    ! average value, if homogeneous (default=0)
     REAL(wp) :: std                    ! standard deviation, if homogeneous (default=0)
     REAL(wp) :: min                    ! lower bound for bounded or cyclic distributions (defaut=none)
     REAL(wp) :: max                    ! upper bound for bounded or cyclic distributions (defaut=none)
     REAL(wp) :: corr_t                 ! correlation timescale, if homogeneous (default=0)
     REAL(wp) :: corr_xy                ! correlation length, if homogeneous/isotropic (default=0)
     CHARACTER(len=lc) :: corr_t_model  ! type dependent, e.g. type of spectrum to use
     CHARACTER(len=lc) :: corr_xy_model ! type dependent, e.g. which diffusive model
     CHARACTER(len=lc) :: corr_xy_file  ! file to read from
     CHARACTER(len=lc) :: corr_z_file   ! file to read from
     CHARACTER(len=lc) :: corr_t_file   ! file to read from
     CHARACTER(len=lc) :: ave_file      ! file to read from
     CHARACTER(len=lc) :: std_file      ! file to read from
     INTEGER :: nar_order               ! order of autoregressive processes (default=1)
     INTEGER :: nar_update              ! frequency of update of autoregressive processes (default=1)
     INTEGER :: diff_type               ! type of diffusion operator to use in 'diff' method
     INTEGER :: diff_passes             ! number of passes of the diffusive model (default=0)
     INTEGER :: ker_type                ! type of kernel to use in 'kernel' method
     INTEGER :: ker_coord               ! type of coordinates to use in 'kernel' method
     REAL(wp) :: ker_density            ! density of kernels in 'kernel' method
     INTEGER :: perlin_n                ! size of subgrid on x for perlin noise
     INTEGER :: perlin_m                ! size of subgrid on y for perlin noise
     INTEGER :: perlin_octaves          ! Number of octaves for perlin noise
     INTEGER :: perlin_lacunarity       ! Lacunarity for perlin noise
     LOGICAL :: perlin_periodic(2)      ! Periodicity flag for perlin noise
     REAL(wp):: perlin_persistency      ! Persistency for perlin noise
     REAL(wp):: perlin_dcoast           ! Distance to coast decay for perlin noise
     ! derived parameters useful to user (from primary parameters in sto_arrays_init)
     INTEGER :: dim   ! dimension of requested array (derived from types) 
     INTEGER :: index ! index in internal arrays
     REAL(wp), DIMENSION(:,:),   POINTER :: sto2d ! 2D stochastic array (points to slice of sto2d)
     REAL(wp), DIMENSION(:,:,:), POINTER :: sto3d ! 3D stochastic array (points to slice of sto3d)
     REAL(wp),                   POINTER :: sto0d ! 0D stochastic array (points to value in sto0d)
   END TYPE

   ! Arrays with all features of stochastic fields
   TYPE(sto_fields), DIMENSION(:), PUBLIC, allocatable :: stofields

   ! Size of stochastic arrays
   INTEGER, PUBLIC :: jpsto = 0            ! total number of requested stochastic fields
   INTEGER, PUBLIC :: jpstomax = 0         ! maximum number of rrquested stochastic fields
   INTEGER, PUBLIC :: jpsto0d = 0          ! number of 0D stochastic parameters
   INTEGER, PUBLIC :: jpsto2d = 0          ! number of 2D stochastic parameters
   INTEGER, PUBLIC :: jpsto3d = 0          ! number of 3D stochastic parameters
   INTEGER, PUBLIC :: jpidx0d = 1          ! number of time slices in 0D stochastic parameters
   INTEGER, PUBLIC :: jpidx2d = 1          ! number of time slices in 2D stochastic parameters
   INTEGER, PUBLIC :: jpidx3d = 1          ! number of time slices in 3D stochastic parameters
   INTEGER, PUBLIC :: jpidxsup2d = 0       ! supplementary slice for transformed field
   INTEGER, PUBLIC :: jpidxsup3d = 0       ! supplementary slice for transformed field
   INTEGER, PUBLIC :: jpidxsup0d = 0       ! supplementary slice for transformed field
   INTEGER, PUBLIC :: jpidxlast2d = 1      ! last slice for transformed field
   INTEGER, PUBLIC :: jpidxlast3d = 1      ! last slice for transformed field
   INTEGER, PUBLIC :: jpidxlast0d = 1      ! last slice for transformed field

   ! Arrays with stochastic fields and features
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:),   ALLOCATABLE, TARGET :: sto2d ! 2D stochastic parameters
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET :: sto3d ! 3D stochastic parameters
   REAL(wp), PUBLIC, DIMENSION(:,:),       ALLOCATABLE, TARGET :: sto0d ! 0D stochastic parameters

   REAL(wp), PUBLIC, DIMENSION(:,:),     ALLOCATABLE :: sto_tmp_2d ! temporary workspace
   REAL(wp), PUBLIC                                  :: sto_tmp_0d ! temporary workspace

   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_idx  ! index of field in stofield array
   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_idx  ! index of field in stofield array
   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto0d_idx  ! index of field in stofield array

   REAL(wp), PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_tcor ! time correlation (for 2D arrays)
   REAL(wp), PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_tcor ! time correlation (for 3D arrays)
   REAL(wp), PUBLIC, DIMENSION(:), ALLOCATABLE :: sto0d_tcor ! time correlation (for 0D arrays)

   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_ord  ! order of autoregressive process
   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_ord  ! order of autoregressive process
   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto0d_ord  ! order of autoregressive process

   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_xy   ! type of horizontal structure
   INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_xy   ! type of horizontal structure


   CHARACTER(len=1), PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_typ  ! nature of grid point (T, U, V, W, F, I)
   CHARACTER(len=1), PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_typ  ! nature of grid point (T, U, V, W, F, I)

   REAL(wp),         PUBLIC, DIMENSION(:), ALLOCATABLE :: sto2d_sgn  ! control of the sign accross the north fold
   REAL(wp),         PUBLIC, DIMENSION(:), ALLOCATABLE :: sto3d_sgn  ! control of the sign accross the north fold

   REAL(wp),         PUBLIC, DIMENSION(:), ALLOCATABLE :: sto_fac    ! factor to restore std

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stopar.F90 13255 2020-07-06 15:41:29Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sto_param_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_param_init  ***
      !!
      !! ** Purpose :   read main block of parameter from namelist
      !!----------------------------------------------------------------------
      INTEGER :: ios
      ! Namelist with general parameters for stochastic modules
      NAMELIST/namsto/ ln_rststo, ln_rstseed, cn_storst_in, cn_storst_out, ln_rng_test, &
                       ln_sto_eos, ln_sto_sppt, ln_sto_spp, ln_sto_skeb, ln_sto_lug

      ! Load namelist for stochastic parameterization
#if ! defined key_ensemble
      CALL load_namelist_sto
#endif

      ! Read namsto namelist
      READ_NML_REF(numnam_sto,namsto)
      READ_NML_CFG(numnam_sto,namsto)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sto_param_init : stochastic parameterization'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto : stochastic parameterization'
         WRITE(numout,*) '      restart stochastic parameters              ln_rststo     = ', ln_rststo
         WRITE(numout,*) '      read seed of RNG from restart file         ln_rstseed    = ', ln_rstseed
         WRITE(numout,*) '      suffix of sto restart name (input)         cn_storst_in  = ', cn_storst_in
         WRITE(numout,*) '      suffix of sto restart name (output)        cn_storst_out = ', cn_storst_out
         WRITE(numout,*) '      stochastic equation of state               ln_sto_eos    = ', ln_sto_eos
         WRITE(numout,*) '      stochastic SPPT scheme                     ln_sto_sppt   = ', ln_sto_sppt
         WRITE(numout,*) '      stochastic SPP scheme                      ln_sto_spp    = ', ln_sto_spp
         WRITE(numout,*) '      stochastic SKEB scheme                     ln_sto_skeb   = ', ln_sto_skeb
         WRITE(numout,*) '      location uncertainty by grid perturbations ln_sto_lug    = ', ln_sto_lug
         WRITE(numout,*) '      perform test of random number generator    ln_rng_test   = ', ln_rng_test
      ENDIF

   END SUBROUTINE sto_param_init


   SUBROUTINE sto_array_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_array_init  ***
      !!
      !! ** Purpose :   allocate stochastic arrays
      !!                initialize arrays with features of stochastic processes
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                            ! Local integer output status for namelist read
      INTEGER  ::   jsto, jsto2d, jsto3d, jsto0d   ! Index of stochastic field
      INTEGER  ::   jord                           ! Order of stochastic field
      LOGICAL  ::   marginal                       ! do we make call to marginal transformation

      ! Check that requests for stochatic fields have been made before
      IF (jpsto==0) THEN
         STOP 'Error in sto_array_init: sto_array_request_new has not been called'
      ENDIF

      ! define name of input restart file (if ensemble restart requred)
      IF(ln_ens_rst_in) cn_storst_in = cn_member//cn_storst_in

      ! Parameter print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sto_array_init : stochastic parameterization'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto : stochastic parameterization'
         WRITE(numout,*) '      restart stochastic parameters           ln_rststo     = ', ln_rststo
         WRITE(numout,*) '      read seed of RNG from restart file      ln_rstseed    = ', ln_rstseed
         WRITE(numout,*) '      suffix of sto restart name (input)      cn_storst_in  = ', trim(cn_storst_in)
         WRITE(numout,*) '      suffix of sto restart name (output)     cn_storst_out = ', trim(cn_storst_out)
         WRITE(numout,*) ' '
      ENDIF

      ! Loop on stochastic fields and set derived parameters
      ! (dimension of array and index in the memory arrays)
      jpsto2d=0 ; jpsto3d=0 ; jpsto0d=0

      DO jsto=1, jpsto

         ! Set index of stochastic field
         IF (stofields(jsto)%dim==2) THEN
            jpsto2d = jpsto2d + stofields(jsto)%nar_order
            stofields(jsto)%index = jpsto2d
         ELSEIF (stofields(jsto)%dim==3) THEN
            jpsto3d = jpsto3d + stofields(jsto)%nar_order
            stofields(jsto)%index = jpsto3d
         ELSEIF (stofields(jsto)%dim==0) THEN
            jpsto0d = jpsto0d + stofields(jsto)%nar_order
            stofields(jsto)%index = jpsto0d
         ENDIF

         ! Set number of time slices (assuming linear interpolation for now)
         IF (stofields(jsto)%dim==2) THEN
            IF (stofields(jsto)%nar_update>1) THEN
               jpidx2d=2 ; jpidxsup2d=1
            ENDIF
         ELSEIF (stofields(jsto)%dim==3) THEN
            IF (stofields(jsto)%nar_update>1) THEN
               jpidx3d=2 ; jpidxsup3d=1
            ENDIF
         ELSEIF (stofields(jsto)%dim==0) THEN
            IF (stofields(jsto)%nar_update>1) THEN
               jpidx0d=2 ; jpidxsup0d=1
            ENDIF
         ENDIF

         ! Do we make calls to the marginal transformation routines
         marginal = stofields(jsto)%type_variate /= 'normal'
         marginal = marginal .OR. stofields(jsto)%ave /= 0._wp
         marginal = marginal .OR. stofields(jsto)%std /= 1._wp
         marginal = marginal .OR. stofields(jsto)%min /= stofields(jsto)%max
         IF (marginal) THEN
            IF (stofields(jsto)%dim==2) THEN
               jpidxsup2d=1
            ELSEIF (stofields(jsto)%dim==3) THEN
               jpidxsup3d=1
            ELSEIF (stofields(jsto)%dim==0) THEN
               jpidxsup0d=1
            ENDIF
         ENDIF

         ! Check user defined options
         CALL check_options(jsto)

      ENDDO

      ! Allocate global stochastic arrays
      IF ( jpsto > 0 ) THEN
         ALLOCATE ( sto_fac(jpsto) )
      ENDIF

      ! Allocate 2D stochastic arrays
      IF ( jpsto2d > 0 ) THEN
         jpidxlast2d = jpidx2d+jpidxsup2d
         ALLOCATE ( sto2d(jpi,jpj,jpidxlast2d,jpsto2d) )
         ALLOCATE ( sto2d_tcor(jpsto2d) )
         ALLOCATE ( sto2d_ord(jpsto2d) )
         ALLOCATE ( sto2d_typ(jpsto2d) )
         ALLOCATE ( sto2d_sgn(jpsto2d) )
         ALLOCATE ( sto2d_idx(jpsto2d) )
         ALLOCATE ( sto2d_xy(jpsto2d) )
      ENDIF

      ! Allocate 3D stochastic arrays
      IF ( jpsto3d > 0 ) THEN
         jpidxlast3d = jpidx3d+jpidxsup3d
         ALLOCATE ( sto3d(jpi,jpj,jpk,jpidxlast3d,jpsto3d) )
         ALLOCATE ( sto3d_tcor(jpsto3d) )
         ALLOCATE ( sto3d_ord(jpsto3d) )
         ALLOCATE ( sto3d_typ(jpsto3d) )
         ALLOCATE ( sto3d_sgn(jpsto3d) )
         ALLOCATE ( sto3d_idx(jpsto3d) )
         ALLOCATE ( sto3d_xy(jpsto3d) )
      ENDIF

      ! Allocate 0D stochastic arrays
      IF ( jpsto0d > 0 ) THEN
         jpidxlast0d = jpidx0d+jpidxsup0d
         ALLOCATE ( sto0d(jpidxlast0d,jpsto0d) )
         ALLOCATE ( sto0d_tcor(jpsto0d) )
         ALLOCATE ( sto0d_ord(jpsto0d) )
         ALLOCATE ( sto0d_idx(jpsto0d) )
      ENDIF

      ! Allocate temporary workspace
      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN
         ALLOCATE ( sto_tmp_2d(jpi,jpj) ) ; sto_tmp_2d(:,:) = 0._wp
      ENDIF

      ! For every stochastic parameter set pointer to memory array
      DO jsto = 1, jpsto
         IF (stofields(jsto)%dim==2) THEN
            stofields(jsto)%sto2d => sto2d( :, :, jpidx2d+jpidxsup2d, stofields(jsto)%index )
            sto2d_idx(stofields(jsto)%index) = jsto
         ELSEIF (stofields(jsto)%dim==3) THEN
            stofields(jsto)%sto3d => sto3d( :, :, :, jpidx3d+jpidxsup3d, stofields(jsto)%index )
            sto3d_idx(stofields(jsto)%index) = jsto
         ELSEIF (stofields(jsto)%dim==0) THEN
            stofields(jsto)%sto0d => sto0d( jpidx0d+jpidxsup0d, stofields(jsto)%index )
            sto0d_idx(stofields(jsto)%index) = jsto
         ENDIF
      ENDDO

      ! For every stochastic parameter:
      ! - set nature of grid point and control of the sign
      !       across the north fold (sto2d_typ, sto2d_sgn)
      ! - set order of every autoregressive process (sto2d_ord)
      ! - set time correlation
      ! - set type of horizontal correlation structure
      ! These arrays are those needed in stopar
      DO jsto = 1, jpsto
         IF (stofields(jsto)%dim==2) THEN
            jsto2d = stofields(jsto)%index
            sto2d_typ(jsto2d)  = stofields(jsto)%type_grid
            sto2d_sgn(jsto2d)  = stofields(jsto)%sign_grid
            sto2d_ord(jsto2d)  = stofields(jsto)%nar_order
            sto2d_tcor(jsto2d) = stofields(jsto)%corr_t

            SELECT CASE(stofields(jsto)%type_xy)
            CASE('white')
              sto2d_xy(jsto2d) = 0
            CASE('diffusive')
              sto2d_xy(jsto2d) = 1
            CASE('kernel')
              sto2d_xy(jsto2d) = 2
            CASE('perlin')
              sto2d_xy(jsto2d) = 3
            CASE DEFAULT
              STOP 'Bad type of horizontal correlation structure in sto_array_init'
            END SELECT

            DO jord = 1, stofields(jsto)%nar_order - 1
               sto2d_ord(jsto2d-jord)  = stofields(jsto)%nar_order - jord
               sto2d_tcor(jsto2d-jord) = stofields(jsto)%corr_t
               sto2d_idx(jsto2d-jord) = sto2d_idx(jsto2d)
               sto2d_xy(jsto2d-jord) = sto2d_xy(jsto2d)
            ENDDO

         ELSEIF (stofields(jsto)%dim==3) THEN
            jsto3d = stofields(jsto)%index
            sto3d_typ(jsto3d)  = stofields(jsto)%type_grid
            sto3d_sgn(jsto3d)  = stofields(jsto)%sign_grid
            sto3d_ord(jsto3d)  = stofields(jsto)%nar_order
            sto3d_tcor(jsto3d) = stofields(jsto)%corr_t

            SELECT CASE(stofields(jsto)%type_xy)
            CASE('white')
              sto3d_xy(jsto3d) = 0
            CASE('diffusive')
              sto3d_xy(jsto3d) = 1
            CASE('kernel')
              sto3d_xy(jsto3d) = 2
            CASE('perlin')
              sto3d_xy(jsto3d) = 3
            CASE DEFAULT
              STOP 'Bad type of horizontal correlation structure in sto_array_init'
            END SELECT

            DO jord = 1, stofields(jsto)%nar_order - 1
               sto3d_ord(jsto3d-jord)  = stofields(jsto)%nar_order - jord
               sto3d_tcor(jsto3d-jord) = stofields(jsto)%corr_t
               sto3d_idx(jsto3d-jord) = sto3d_idx(jsto3d)
               sto3d_xy(jsto3d-jord) = sto3d_xy(jsto3d)
            ENDDO

         ELSEIF (stofields(jsto)%dim==0) THEN
            jsto0d = stofields(jsto)%index
            sto0d_ord(jsto0d)  = stofields(jsto)%nar_order
            sto0d_tcor(jsto0d) = stofields(jsto)%corr_t

            DO jord = 1, stofields(jsto)%nar_order - 1
               sto0d_ord(jsto0d-jord)  = stofields(jsto)%nar_order - jord
               sto0d_tcor(jsto0d-jord) = stofields(jsto)%corr_t
               sto0d_idx(jsto0d-jord) = sto0d_idx(jsto0d)
            ENDDO

         ENDIF
      END DO

   END SUBROUTINE sto_array_init


   SUBROUTINE sto_array_request_size(kjpsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_array_request_size  ***
      !!
      !! ** Purpose :   request maximum number of stochastic fields
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kjpsto

      jpstomax = kjpsto

      ! Allocate array with features of stochastic fields
      IF (jpstomax>0) THEN
         ALLOCATE ( stofields(jpstomax) )
      ELSE
         STOP 'Bad number of requested stochastic fields in sto_array_init'
      ENDIF

    END SUBROUTINE sto_array_request_size


   SUBROUTINE sto_array_request_new(kjsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_array_request_new  ***
      !!
      !! ** Purpose :   request new stochastic field
      !!----------------------------------------------------------------------
      INTEGER, INTENT(out) :: kjsto
      CHARACTER(len=8) :: default_name

      ! Check that request for amximum number of stochatic fields has been made before
      IF (jpstomax==0) THEN
         STOP 'Error in sto_array_init: sto_array_request_new has not been called'
      ENDIF

      ! Attribute new index of stochastic field
      jpsto = jpsto+1
      kjsto = jpsto

      ! Check size of stochastic arrays
      IF (jpsto>jpstomax) THEN
         STOP 'Insufficient size of stochastic arrays in sto_array_init'
      ENDIF

      ! Set default features of stochastic field
      ! They can be modified by the user until sto_array_init is called
      WRITE(default_name,'(a3,i5.5)') 'sto',kjsto
      stofields(kjsto)%stoname=default_name  ! Identifier, for instance: name=<physics>_<index>
      stofields(kjsto)%type_t='white'        ! type of time structure (constant, white, arn, spectrum,...)
      stofields(kjsto)%type_xy='white'       ! type of xy structure (white, diffusive, perlin,...)
      stofields(kjsto)%type_z='white'        ! type of z structure (white, modes,...)
      stofields(kjsto)%type_variate='normal' ! type of marginal pdf(normal,lognormal,gamma,beta,...)
      stofields(kjsto)%type_grid='T'         ! type of grid (T, U, V, F)
      stofields(kjsto)%sign_grid=1._wp       ! sign to use to connect field across domain
      stofields(kjsto)%ave=0._wp             ! average value, if homogeneous (default=0)
      stofields(kjsto)%std=1._wp             ! standard deviation, if homogeneous (default=0)
      stofields(kjsto)%min=0._wp             ! lower bound for bounded or cyclic distribution (default=none)
      stofields(kjsto)%max=0._wp             ! upper bound for bounded or cyclic distribution (default=none)
      stofields(kjsto)%corr_t=0._wp          ! correlation timescale, if homogeneous (default=0)
      stofields(kjsto)%corr_xy=0._wp         ! correlation length, if homogeneous/isotropic (default=0)
      stofields(kjsto)%corr_t_model='none'   ! type dependent, e.g. type of spectrum to use
      stofields(kjsto)%corr_xy_model='none'  ! type dependent, e.g. which diffusive model
      stofields(kjsto)%corr_xy_file='none'   ! file to read from
      stofields(kjsto)%corr_z_file='none'    ! file to read from
      stofields(kjsto)%corr_t_file='none'    ! file to read from
      stofields(kjsto)%ave_file='none'       ! file to read from
      stofields(kjsto)%std_file='none'       ! file to read from
      stofields(kjsto)%nar_order=1           ! order of autoregressive processes (default=1)
      stofields(kjsto)%nar_update=1          ! frequency of update of autoregressive processes (default=1)
      stofields(kjsto)%diff_type=0           ! type of diffusion operator to use in 'diff' method (default=laplacian)
      stofields(kjsto)%diff_passes=0         ! number of passes of the diffusive model (default=0)
      stofields(kjsto)%ker_type=0            ! type of kernel to use in 'kernel' method (default=gaussian)
      stofields(kjsto)%ker_coord=0           ! type of coordinates to use in 'kernel' method (default=grid)
      stofields(kjsto)%ker_density=0._wp     ! density of kernels in 'kernel' method (default=0.0->not used)
      stofields(kjsto)%perlin_n = 0                     ! Perlin noise Default size
      stofields(kjsto)%perlin_m = 0                     !
      stofields(kjsto)%perlin_octaves  = 1              ! Default number of octaves
      stofields(kjsto)%perlin_lacunarity = 0            ! Default parameters for Perlin
      stofields(kjsto)%perlin_periodic(:) = .true.      !
      stofields(kjsto)%perlin_persistency = 0._wp       !
      stofields(kjsto)%perlin_dcoast = 0._wp            ! Default distance to coast

   END SUBROUTINE sto_array_request_new


   SUBROUTINE check_options(kjsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE check_options  ***
      !!
      !! ** Purpose :   check user defined options for stochastic fields
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kjsto

      SELECT CASE(stofields(kjsto)%type_t)
      CASE('white','arn')
      CASE DEFAULT
         STOP 'Bad type of time structure in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%type_xy)
      CASE('constant','white','kernel','diffusive','perlin')
      CASE DEFAULT
         STOP 'Bad type of xy structure in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%type_z)
      CASE('constant','white')
      CASE DEFAULT
         STOP 'Bad type of z structure in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%type_variate)
      CASE('normal','lognormal','wrapped_normal','bounded_atan')
      CASE DEFAULT
         STOP 'Bad type of marginal distribution in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%type_grid)
      CASE('T','U','V','F')
      CASE DEFAULT
         STOP 'Bad type of grid in stoarray'
      END SELECT

      IF (ABS(stofields(kjsto)%sign_grid) .NE. 1._wp ) THEN
         STOP 'Bad type of grid connection in stoarray'
      ENDIF

      SELECT CASE(stofields(kjsto)%corr_t_model)
      CASE('none')
      CASE DEFAULT
         STOP 'Bad type of time correlation model in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%corr_xy_model)
      CASE('none')
      CASE DEFAULT
         STOP 'Bad type of xy correlation model in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%corr_t_file)
      CASE('none')
      CASE DEFAULT
         STOP 'File option is not yet implemented in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%corr_xy_file)
      CASE('none')
      CASE DEFAULT
         STOP 'File option is not yet implemented in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%corr_z_file)
      CASE('none')
      CASE DEFAULT
         STOP 'File option is not yet implemented in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%ave_file)
      CASE('none')
      CASE DEFAULT
         STOP 'File option is not yet implemented in stoarray'
      END SELECT

      SELECT CASE(stofields(kjsto)%std_file)
      CASE('none')
      CASE DEFAULT
         STOP 'File option is not yet implemented in stoarray'
      END SELECT

      IF (stofields(kjsto)%nar_order < 1) THEN
         STOP 'Bad order of AR process in stoarray'
      ENDIF

      IF (stofields(kjsto)%nar_update < 1) THEN
         STOP 'Bad update time of AR process in stoarray'
      ENDIF

      IF ( (stofields(kjsto)%diff_type < 0) .OR. (stofields(kjsto)%diff_type > 1 ) ) THEN
         STOP 'Bad diffusion operator type in stoarray'
      ENDIF

      IF (stofields(kjsto)%diff_passes < 0) THEN
         STOP 'Bad number of passes of the diffusive model in stoarray'
      ENDIF

      IF ( (stofields(kjsto)%ker_type < 0) .OR. (stofields(kjsto)%ker_type > 5 ) ) THEN
         STOP 'Bad kernel type in stoarray'
      ENDIF

      IF ( (stofields(kjsto)%ker_coord < 0) .OR. (stofields(kjsto)%ker_coord > 2 ) ) THEN
         STOP 'Bad coordinate type for kernel method in stoarray'
      ENDIF

      IF (stofields(kjsto)%ker_density < 0._wp ) THEN
         STOP 'Bad coordinate type for kernel method in stoarray'
      ENDIF

      IF (stofields(kjsto)%type_xy == 'perlin' ) THEN
         !
         IF(stofields(kjsto)%perlin_n <= 1 .OR. stofields(kjsto)%perlin_n >= jpiglo .OR. &
            stofields(kjsto)%perlin_m <= 1 .OR. stofields(kjsto)%perlin_m >= jpjglo ) &
         & STOP 'Bad perturbation grid size in Perlin noise'
         !
         IF(stofields(kjsto)%perlin_octaves < 1 .OR. stofields(kjsto)%perlin_octaves > 5 ) &
         & STOP 'Bad number of octaves in Perlin noise'
         !
         IF(stofields(kjsto)%perlin_lacunarity < 1 .OR. stofields(kjsto)%perlin_lacunarity > 5 ) &
         & STOP 'Bad lacunarity in Perlin noise'
         !
         IF(stofields(kjsto)%perlin_persistency <= 0._wp .OR. stofields(kjsto)%perlin_persistency >= 1._wp ) &
         & STOP 'Bad persistency in Perlin noise'
         !
         IF(stofields(kjsto)%perlin_dcoast <= 0._wp .OR. stofields(kjsto)%perlin_dcoast >= 1.e+6_wp ) &
         & STOP 'Bad coastal decay scale in Perlin noise'
         !
      ENDIF

   END SUBROUTINE check_options

END MODULE stoarray

