MODULE stokernel
   !!======================================================================
   !!                       ***  MODULE  stokernel  ***
   !! Purpose : apply kernel method to generate new random field
   !!           with specific horizontal correlation structure
   !!=====================================================================
   !!   sto_kernel      : generate new random field with kernel method
   !!   sto_kernel_init : initialize kernel method
   !!----------------------------------------------------------------------
   USE stoexternal
   USE stoarray
   USE stowhite
   USE stosobolseq

   IMPLICIT NONE
   PRIVATE

   ! User option to randomized length scales of superposed kernels
   LOGICAL, PUBLIC, SAVE  :: kernel_randomized_length_scale = .FALSE.
   REAL(wp), PUBLIC, SAVE :: kernel_randomized_length_scale_std = 0.25_wp

   ! Size of module arrays
   INTEGER, SAVE :: jpker    ! number of superposed kernels
   INTEGER, SAVE :: jpkermax ! maximum number of superposed kernels

   ! Arrays with kernel features
   REAL(wp), DIMENSION(:), ALLOCATABLE :: xkernel, ykernel ! locations of kernels
   REAL(wp), DIMENSION(:), ALLOCATABLE :: ker_lx, ker_ly   ! kernel length scales
   REAL(wp), DIMENSION(:), ALLOCATABLE :: psto_ampl        ! amplitude of perturbations
   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: psto_norm      ! normalizatin factor

   ! Variables for the Sobol sequence generator
   INTEGER, PARAMETER :: sobol_ndims=2
   TYPE(sobol_state), DIMENSION(sobol_ndims) :: rng

   PUBLIC sto_kernel, sto_kernel_init

CONTAINS

   SUBROUTINE sto_kernel(psto,jsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_kernel  ***
      !!
      !! ** Purpose :   generate new random field with kernel method
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) :: psto   ! output stochastic field
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray

      INTEGER :: jker, ji, jj
      REAL(wp) :: kernel_value

      IF (SIZE(psto,1).NE.jpi) STOP 'Bad dimensions in sto_kernel'
      IF (SIZE(psto,2).NE.jpj) STOP 'Bad dimensions in sto_kernel'

      ! Get number of kernels to superpose for this field
      jpker = get_number_of_kernels(jsto)

      ! Sample kernel locations and amplitudes (on first processors)
      IF ( narea == 1 ) THEN ! calculate 1D global values on zero proc and broadcast to the other procs
         DO jker = 1, jpker
            ! get next kernel location
            CALL get_kernel_location(xkernel(jker),ykernel(jker),jsto)
            ! set kernel length scales
            CALL set_length_scales(ker_lx(jker),ker_ly(jker),jsto)
         END DO

         ! sample Gaussian N(0,1) amplitudes
         CALL sto_white( psto1d = psto_ampl )
      ENDIF

      ! Broadcast random locations, length scales and amplitudes to all processors
      CALL broadcast_array( xkernel )
      CALL broadcast_array( ykernel )
      CALL broadcast_array( ker_lx )
      CALL broadcast_array( ker_ly )
      CALL broadcast_array( psto_ampl )

      ! Compute stochastic fields by superposing the kernels at random locations
      ! with random normal amplitude
      psto(:,:) = 0._wp ; psto_norm(:,:) = 0._wp
      DO jker = 1, jpker ! loop on kernels
         DO jj = 1, jpj
         DO ji = 1, jpi
            kernel_value = kernel(ji,jj,jker,jsto)
            psto(ji,jj) = psto(ji,jj) + psto_ampl(jker) * kernel_value
            psto_norm(ji,jj) = psto_norm(ji,jj) + kernel_value * kernel_value
         END DO
         END DO
      END DO

      ! Normalize stochastic field
      WHERE (psto_norm(:,:) /= 0) psto(:,:) = psto(:,:) / SQRT( psto_norm(:,:) )

   END SUBROUTINE sto_kernel


   SUBROUTINE sto_kernel_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_kernel_init  ***
      !!
      !! ** Purpose :   initialize kernel method
      !!----------------------------------------------------------------------
      INTEGER :: jsto

      ! Compute maximum number of kernels to superpose
      jpkermax = 0
      DO jsto=1, jpsto                                   ! loop on all stochastic fields
         IF (stofields(jsto)%type_xy == 'kernel' ) THEN  ! with kernel option
            ! Get number of kernels required for this field
            jpker = get_number_of_kernels(jsto) 
            jpkermax = MAX(jpker,jpkermax)
            ! Check that grid arrays have correct size (if required)
            IF (stofields(jsto)%ker_coord == 0 ) THEN
              IF (SIZE(mig,1).NE.jpi) STOP 'Incorrect grid in stokernel'
              IF (SIZE(mjg,1).NE.jpj) STOP 'Incorrect grid in stokernel'
            ELSE
              IF (SIZE(glamt,1).NE.jpi) STOP 'Incorrect grid in stokernel'
              IF (SIZE(glamt,2).NE.jpj) STOP 'Incorrect grid in stokernel'
              IF (SIZE(gphit,1).NE.jpi) STOP 'Incorrect grid in stokernel'
              IF (SIZE(gphit,2).NE.jpj) STOP 'Incorrect grid in stokernel'
            ENDIF
         ENDIF
      ENDDO

      IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sto_kernel_init : initialization of kernel method'
            WRITE(numout,*) '~~~~~~~~~~~~~~~'
            WRITE(numout,*)
            WRITE(numout,*) ' -- maximum number of superposed kernels -- '
            WRITE(numout,*) '    jpkermax = ', jpkermax
      ENDIF

      ! allocate location of kernels
      ALLOCATE( xkernel(jpkermax) )
      ALLOCATE( ykernel(jpkermax) )
      ! allocate kernel length scales
      ALLOCATE( ker_lx(jpkermax) )
      ALLOCATE( ker_ly(jpkermax) )
      ! allocate amplitudes
      ALLOCATE( psto_ampl(jpkermax) )
      ! allocate normalization factor
      ALLOCATE( psto_norm(jpi,jpj) )

      ! initialize sobol sequence
      ! (used to generate the quasi random locations of the kernels)
      CALL initialize_sobol_sequence()

   END SUBROUTINE sto_kernel_init


   FUNCTION get_number_of_kernels(jsto)
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION get_number_of_kernels  ***
      !!
      !! ** Purpose :   return the number of kernel to use for a given stochastic field
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: jsto
      INTEGER :: get_number_of_kernels

      INTEGER, PARAMETER :: nker_min = 4  ! minimum number of superposed kernels
      INTEGER :: nker_max                 ! maximum number of superposed kernels
      INTEGER :: nker

      REAL(wp) :: lx, ly, sizex, sizey
      INTEGER :: coord_type

      ! Get type of coordinates to use
      coord_type = stofields(jsto)%ker_coord

      ! Set maximum number of kernels to grid size
      nker_max = jpiglo * jpjglo

      ! Transform length scales in grid points
      ! to get the number of kernels per grid point
      lx = stofields(jsto)%corr_xy
      ly = stofields(jsto)%corr_xy
      IF (coord_type>0) THEN        ! grid coordinates
         sizex = glamtglomax - glamtglomin
         sizey = gphitglomax - gphitglomin
         IF ( (sizex==0._wp) .OR. (sizey==0._wp) ) THEN
             STOP 'Bad coordinats in stokernel'
         ENDIF
         lx = lx * REAL(jpiglo,wp) / sizex
         ly = ly * REAL(jpjglo,wp) / sizey
      ENDIF

      ! Option to modify the kernel density
      IF (stofields(jsto)%ker_density > 0._wp) THEN
         lx = lx / stofields(jsto)%ker_density
         ly = ly / stofields(jsto)%ker_density
      ENDIF

      ! Calculate number of kernels to superpose
      lx = MAX( lx, 1._wp ) ; ly = MAX( ly, 1._wp )
      nker = NINT( REAL(nker_max,wp) / ( lx * ly ) )

      ! Bound the number of kernels between min and max
      get_number_of_kernels = MAX( MIN(nker,nker_max), nker_min )

   END FUNCTION get_number_of_kernels


   FUNCTION kernel(ji,jj,jker,jsto)
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION kernel  ***
      !!
      !! ** Purpose :   return the value of the kernel function
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: ji,jj,jker,jsto
      REAL(wp) :: kernel
      
      INTEGER :: ker_type
      REAL(wp) :: ds2

      ! Get type of kernel
      ker_type = stofields(jsto)%ker_type

      ! Get scaled square distance between kernel location and grid point
      ds2 = get_distance(ji,jj,jker,jsto)

      ! Evaluate kernel function
      IF (ker_type==0) THEN
         ! Gaussian kernel
         kernel = EXP( - 0.5_wp * ds2 ) 
      ELSEIF (ker_type==1) THEN
         ! Laplacian kernel
         kernel = 1._wp / ( 1._wp + ds2 )
      ELSEIF (ker_type==2) THEN
         ! Box kernel
         kernel = 0.
         IF (ds2 < 1._wp) kernel = 1._wp
      ELSEIF (ker_type==3) THEN
         ! Triangle kernel
         kernel = 0.
         IF (ds2 < 1._wp) kernel = 1._wp - SQRT(ds2)
      ELSEIF (ker_type==4) THEN
         ! Mexican hat wavelet (Ricker wavelet)
         kernel = ( 1._wp - ds2 ) * EXP( - 0.5_wp * ds2 ) 
      ELSEIF (ker_type==5) THEN
         ! Morlet wavelet (with specific choice of frequency, adjust if needed)
         kernel = EXP( - 0.5_wp * ds2 ) * COS( SQRT(ds2) )
      ELSE
         STOP 'Bad kernel type in stokernel'
      ENDIF

   END FUNCTION kernel


   FUNCTION get_distance(ji,jj,jker,jsto)
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION get_distance  ***
      !!
      !! ** Purpose :   return squared scaled distance
      !!                between grid point and kernel location
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: ji,jj,jker,jsto
      REAL(wp) :: get_distance

      INTEGER :: coord_type
      REAL(wp) :: dx, dy, ds2, ds, az, lx, ly, l2

      ! Get type of coordinates to use
      coord_type = stofields(jsto)%ker_coord

      ! Compute square distance
      IF (coord_type==0) THEN
         ! distance in grid points
         ! -----------------------
         dx = REAL(mig(ji,nn_hls),wp) - xkernel(jker)
         dy = REAL(mjg(jj,nn_hls),wp) - ykernel(jker)
         ! scale with length scales
         dx = dx / ker_lx(jker)
         dy = dy / ker_ly(jker)
         ds2 = dx * dx + dy * dy
      ELSEIF (coord_type==1) THEN
         ! Cartesian coordinates
         ! ---------------------
         dx = glamt(ji,jj) - xkernel(jker)
         dy = gphit(ji,jj) - ykernel(jker)
         ! scale with length scales
         dx = dx / ker_lx(jker)
         dy = dy / ker_ly(jker)
         ds2 = dx * dx + dy * dy
      ELSEIF (coord_type==2) THEN
         ! spherical coordinates
         ! ---------------------
         IF (kernel_randomized_length_scale) THEN
            ! get distance (ds) and azimuth (az) between the two points
            ds = sph_distance(glamt(ji,jj),gphit(ji,jj),xkernel(jker),ykernel(jker),azimuth=az)
            ! compute length scale in direction az
            lx = ker_lx(jker) * SIN(az)
            ly = ker_ly(jker) * COS(az)
            l2 = lx * lx + ly * ly
         ELSE
            ! get distance (ds) between the two points
            ds = sph_distance(glamt(ji,jj),gphit(ji,jj),xkernel(jker),ykernel(jker))
            ! compute length scale
            lx = ker_lx(jker)
            ly = ker_ly(jker)
            l2 = lx * lx + ly * ly
         ENDIF
         ! scale with length scales
         ds2 = ds * ds / l2
      ELSE
         STOP 'Bad coordinate type in stokernel'
      ENDIF

      get_distance = ds2

   END FUNCTION get_distance
   

   SUBROUTINE get_kernel_location(x,y,jsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE get_kernel_location  ***
      !!
      !! ** Purpose :   get next location of superposed kernel
      !!                (uniformly placed in global domain)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: x, y
      INTEGER, INTENT(in) :: jsto

      INTEGER  :: coord_type
      INTEGER  :: ix, iy
      REAL(wp) :: rx, ry
      REAL(wp), DIMENSION(2)  :: u

      ! Get type of coordinates to use
      coord_type = stofields(jsto)%ker_coord

      ! get next kernel location in grid
      ! from quasi random Sobol sequence
      ! (uniform between 0 and 1 in a 2D square)
      u(1) = rng(1)%next()
      u(2) = rng(2)%next()

      ! get location in the required coordinate system
      IF (coord_type == 0) THEN
         ! coordinates in grid points
         x = u(1) * (jpiglo+1)
         y = u(2) * (jpjglo+1)
      ELSEIF (coord_type > 0) THEN
         ! Cartesian and spherical coordinates
         x = glamtglomin + u(1) * (glamtglomax - glamtglomin)
         y = gphitglomin + u(2) * (gphitglomax - gphitglomin)
      ELSE
         STOP 'Bad coordinate type in stokernel'
      ENDIF

   END SUBROUTINE get_kernel_location


   SUBROUTINE set_length_scales(lx,ly,jsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE set_length_scales  ***
      !!
      !! ** Purpose :   set length scales of kernel functions
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: lx, ly    ! length scale
      INTEGER, INTENT(in) :: jsto        ! index of stochastic field

      REAL(wp), DIMENSION(2) :: gran

      ! Get length scales for field with index jsto
      ! (provisionnally isotropic)
      lx = MAX( stofields(jsto)%corr_xy, 1._wp )
      ly = MAX( stofields(jsto)%corr_xy, 1._wp )

      ! Randomize length scale if required
      ! (this makes the individual kernels anisotropic)
      IF (kernel_randomized_length_scale) THEN
         CALL sto_white( psto1d = gran )
         lx = lx * EXP(kernel_randomized_length_scale_std * gran(1))
         ly = ly * EXP(kernel_randomized_length_scale_std * gran(2))
      ENDIF

   END SUBROUTINE set_length_scales


   SUBROUTINE initialize_sobol_sequence()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE initialize_sobol_sequence  ***
      !!
      !! ** Purpose :   initialize sobol sequence algorithm
      !!                see https://github.com/DaanVanVugt/sobseq/
      !!                for more details     
      !!----------------------------------------------------------------------
      ! Parameters needed in 2 dimensions
      INTEGER, parameter, dimension(1:sobol_ndims) :: s = (/1,2/)
      INTEGER, parameter, dimension(1:sobol_ndims) :: a = (/0,1/)
      INTEGER, parameter, dimension(5,sobol_ndims) :: m = reshape((/1,0,0,0,0, &
                                                              1,3,0,0,0/), (/5,2/))
      INTEGER :: i

      ! Initialization in 2 dimensions
      DO i=1, sobol_ndims
         CALL rng(i)%initialize(s(i), a(i), m(:,i))
      ENDDO

   END SUBROUTINE initialize_sobol_sequence


   FUNCTION sph_distance(lon1d, lat1d, lon2d, lat2d, azimuth)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sph_distance  ***
      !!
      !! ** Purpose :   compute distance between two location on the sphere
      !!                (with option to compute azimuth)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(IN) :: lon1d, lat1d, lon2d, lat2d
      REAL(wp), INTENT(OUT), OPTIONAL :: azimuth
      REAL(wp) :: sph_distance

      ! REAL(wp), PARAMETER :: EarthRadius = 6371.0 ! Earth radius in kilometers
      REAL(wp), PARAMETER :: pi = 3.14159265358979323846_wp
      REAL(wp), PARAMETER :: toRadians = pi / 180.0_wp
      REAL(wp) :: lon1, lat1, lon2, lat2

      REAL(wp) :: dlat, dlon, a, c, distance

      ! Convert latitude and longitude from degrees to radians
      lat1 = toRadians*lat1d
      lon1 = toRadians*lon1d
      lat2 = toRadians*lat2d
      lon2 = toRadians*lon2d

      ! Compute differences
      dlat = lat2 - lat1
      dlon = lon2 - lon1

      ! Haversine formula
      a = SIN(dlat / 2)**2 + COS(lat1) * COS(lat2) * SIN(dlon / 2)**2
      c = 2 * ATAN2(SQRT(a), SQRT(1 - a))

      ! Azimuth calculation
      IF (PRESENT(azimuth)) THEN
         azimuth = ATAN2(SIN(dlon) * COS(lat2), COS(lat1) * SIN(lat2) - SIN(lat1) * COS(lat2) * COS(dlon))
      ENDIF

      !sph_distance = EarthRadius * c
      sph_distance = c / toRadians

   END FUNCTION sph_distance

END MODULE stokernel
