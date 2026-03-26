MODULE stoexternal
   !!======================================================================
   !!                       ***  MODULE  stoexternal  ***
   !! Purpose        : external resources provide by the model (user supplied)
   !!=====================================================================
   USE par_oce         ! ocean parameters
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE lib_mpp
   USE lib_fortran     ! to use sign with key_nosignedzero
   USE mppens
#if ! defined key_mpi_off
   USE MPI
#endif


   ! Description of the global grid
   REAL(wp), PUBLIC, SAVE :: glamtglomin, gphitglomin   ! minimum of global latitude and longitude
   REAL(wp), PUBLIC, SAVE :: glamtglomax, gphitglomax   ! maximum of global latitude and longitude

   ! Description of the mask (pointers used in stochastic codes)
   LOGICAL, PUBLIC, SAVE  :: use_mask3d = .TRUE.  ! use 2D or 3D masks in stochastic code
   INTEGER, PUBLIC, SAVE  :: grid_type = 0        ! position of u,v points relative to r,t points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:),   POINTER :: rmask2d ! land/ocean mask at T-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:),   POINTER :: umask2d ! land/ocean mask at U-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:),   POINTER :: vmask2d ! land/ocean mask at V-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:),   POINTER :: fmask2d ! land/ocean mask at F-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: rmask3d ! land/ocean mask at T-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: umask3d ! land/ocean mask at U-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: vmask3d ! land/ocean mask at V-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: fmask3d ! land/ocean mask at F-points

   PUBLIC load_namelist_sto, initialize_grid, initialize_mask, broadcast_array

CONTAINS

   SUBROUTINE load_namelist_sto

      ! Load full namelist in memory
      CALL load_nml( numnam_sto_ref, 'namelist_sto_ref' , numout, lwm )
      CALL load_nml( numnam_sto_cfg, 'namelist_sto_cfg' , numout, lwm )

   END SUBROUTINE load_namelist_sto


   SUBROUTINE initialize_grid
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE initialize_grid  ***
      !!
      !! ** Purpose :   initialization of grid features
      !!----------------------------------------------------------------------
      INTEGER :: iflag, mpi_type

      glamtglomin = MINVAL(glamt)
      gphitglomin = MINVAL(gphit)
      glamtglomax = MAXVAL(glamt)
      gphitglomax = MAXVAL(gphit)

#if ! defined key_mpi_off
      IF (wp == dp) THEN
         mpi_type = mpi_double_precision
      ELSE
         mpi_type = mpi_real
      ENDIF

      CALL MPI_ALLREDUCE( MPI_IN_PLACE, glamtglomin, 1, mpi_type, MPI_MIN, mpi_comm_oce, iflag)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, gphitglomin, 1, mpi_type, MPI_MIN, mpi_comm_oce, iflag)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, glamtglomax, 1, mpi_type, MPI_MAX, mpi_comm_oce, iflag)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, gphitglomax, 1, mpi_type, MPI_MAX, mpi_comm_oce, iflag)
#endif

   END SUBROUTINE initialize_grid


   SUBROUTINE initialize_mask
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE initialize_mask  ***
      !!
      !! ** Purpose :   initialization of mask features
      !!----------------------------------------------------------------------

      rmask3d => tmask
      umask3d => umask
      vmask3d => vmask
      fmask3d => fmask

   END SUBROUTINE initialize_mask


   SUBROUTINE broadcast_array( ptab )
      REAL(wp), DIMENSION(:), INTENT(in) :: ptab   ! array to broadcast

      INTEGER :: iflag, mpi_type

#if ! defined key_mpi_off
      IF (wp == dp) THEN
         mpi_type = mpi_double_precision
      ELSE
         mpi_type = mpi_real
      ENDIF

      ! Insert MPI broadcast code here
      CALL MPI_BCAST(ptab, SIZE(ptab,1), mpi_type, 0, mpi_comm_oce, iflag)
      CALL MPI_BARRIER(mpi_comm_oce, iflag)
#endif

   END SUBROUTINE broadcast_array

END MODULE stoexternal
