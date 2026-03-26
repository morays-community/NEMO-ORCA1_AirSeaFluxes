MODULE mppens
   !!======================================================================
   !!                       ***  MODULE  mppens  ***
   !! Ensemble simulation by allocating a communicator to each ensemble member
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   mpp_start_ensemble : divide input communicator in ensemble members and assign local communicator
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   USE timing
   USE lib_mpp
#if ! defined key_mpi_off
   USE MPI
#endif

   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC   mpp_start_ensemble

   ! Parameters from namelist
   LOGICAL, PUBLIC    :: ln_ensemble = .FALSE.   ! perform ensemble simulation (T/F)
   LOGICAL, PUBLIC    :: ln_ens_diag = .FALSE.   ! ensemble online diagnostic
   LOGICAL, PUBLIC    :: ln_ens_rst_in = .FALSE. ! use ensemble (T) or single (F) input restart file
   INTEGER, PUBLIC    :: nn_ens_size = 1         ! ensemble size
   INTEGER, PUBLIC    :: nn_ens_start = 1        ! index of the first ensemble member

   ! MPI communicators
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE         :: ncomm_member   ! communicator for every ensemble member
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE :: ncomm_area     ! communicator for every subdomain

   ! Member index
   INTEGER, PUBLIC ::   nmember = 1       !: ensemble member index
   CHARACTER(len=3), PUBLIC :: cn_member  !: character string with ensemble member index

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE mpp_start_ensemble( localComm )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_start_ensemble  ***
      !!
      !! ** Purpose : initialize MPI if not yet done, and then
      !!              divide input communicator (localComm) in ensemble members
      !!              and call mpp_start for member communicator (memberComm)
      !!----------------------------------------------------------------------
      INTEGER         , OPTIONAL   , INTENT(in   ) ::   localComm    ! available input communicator
      !
      INTEGER ::   memberComm   ! communicator for current member
      INTEGER ::   icomm_global, impprank, ierr, ios
      LOGICAL ::   llmpi_init
      !
      NAMELIST/namens/ ln_ensemble, ln_ens_rst_in, nn_ens_size, nn_ens_start, ln_ens_diag
      !!----------------------------------------------------------------------
#if ! defined key_mpi_off
      !
      CALL mpi_initialized ( llmpi_init, ierr )
      IF( ierr /= MPI_SUCCESS ) CALL ctl_stop( 'STOP', ' lib_mpp: Error in routine mpi_initialized' )

      ! initialize MPI if not yet done
      IF( .NOT. llmpi_init ) THEN
         IF( PRESENT(localComm) ) THEN
            WRITE(ctmp1,*) ' lib_mpp: You cannot provide a local communicator '
            WRITE(ctmp2,*) '          without calling MPI_Init before ! '
            CALL ctl_stop( 'STOP', ctmp1, ctmp2 )
         ENDIF
         CALL mpi_init( ierr )
         IF( ierr /= MPI_SUCCESS ) CALL ctl_stop( 'STOP', ' lib_mpp: Error in routine mpi_init' )
      ENDIF

      ! set available global communicator for ensemble simulation (icomm_global)
      IF( .NOT. PRESENT(localComm) ) THEN
         CALL mpi_comm_dup( mpi_comm_world, icomm_global, ierr)
         IF( ierr /= MPI_SUCCESS ) CALL ctl_stop( 'STOP', ' lib_mpp: Error in routine mpi_comm_dup' )
      ELSE
         icomm_global = localComm
      ENDIF

      ! read namelist and broadcast to all domains of all members
      mpi_comm_oce = icomm_global
      CALL mpi_comm_rank( localComm, impprank, ierr )
      CALL load_nml( numnam_sto_ref, 'namelist_sto_ref', -1, impprank == 0 )
      CALL load_nml( numnam_sto_cfg, 'namelist_sto_cfg', -1, impprank == 0 )
      READ_NML_REF(numnam_sto,namens)
      READ_NML_CFG(numnam_sto,namens)

      IF(impprank == 0) THEN
         PRINT *, '--- Ensemble configuration ---'
         PRINT *, ' ln_ensemble:   ', ln_ensemble
         PRINT *, ' ln_ens_rst_in: ', ln_ens_rst_in
         PRINT *, ' nn_ens_size:   ', nn_ens_size
         PRINT *, ' nn_ens_start:  ', nn_ens_start
         PRINT *, ' ln_ens_diag:   ', ln_ens_diag
         PRINT *, '------------------------------'
      ENDIF

      ! set ensemble communicators
      ! assign current communicator and current member index
      IF ( ln_ensemble ) THEN
         IF( Agrif_Root() ) THEN
            CALL mpp_set_ensemble( nmember, memberComm, icomm_global )
# if defined key_agrif
         ELSE
            nmember = Agrif_Parent(nmember)
# endif
         ENDIF
      ELSE
         nmember = 1
         memberComm = icomm_global
      ENDIF

      IF ( ln_ensemble ) THEN
         WRITE(cn_member,'(I3.3)') nmember
         cxios_context = TRIM(cn_member)//TRIM(cxios_context)
      ENDIF

      ! Start communicator for current member
      CALL mpp_start( memberComm )

#else
      nmember = 1
      IF( PRESENT( localComm ) ) THEN
         CALL mpp_start( localComm )
      ELSE
         CALL mpp_start( )
      ENDIF
#endif

   END SUBROUTINE mpp_start_ensemble


   SUBROUTINE mpp_set_ensemble( kmember, memberComm, localComm )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_set_ensemble  ***
      !!
      !! ** Purpose : divide input communicator (localComm) in ensemble members
      !!              and assign local communicator (memberComm)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(out  ) ::   kmember      ! current member index
      INTEGER, INTENT(out  ) ::   memberComm   ! communicator for current member
      INTEGER, INTENT(in   ) ::   localComm    ! available input communicator
      !
      INTEGER ::   igrp_global, imember, imppsize, impprank, ierr
      INTEGER ::   jmember, jjproc, jpproc, jarea
      LOGICAL ::   llmpi_init
      !
      INTEGER, DIMENSION(:), ALLOCATABLE ::   irank_member ! list of processors for current member (dim: jpnij)
      INTEGER, DIMENSION(:), ALLOCATABLE ::   igrp_member  ! group ID for every ensemble member
      INTEGER, DIMENSION(:), ALLOCATABLE ::   irank_area   ! list of processors for current subdomain (dim: nn_ens_size)
      INTEGER, DIMENSION(:), ALLOCATABLE ::   igrp_area    ! group ID for every subdomain

      !!----------------------------------------------------------------------

      ! compute number of processors available for each member (jpproc)
      CALL mpi_comm_size( localComm, imppsize, ierr )
      IF( MOD(imppsize,nn_ens_size) /= 0 ) CALL ctl_stop( 'STOP', ' lib_mpp: nbr of proc not multiple of ensemble size' )
      jpproc = imppsize / nn_ens_size

      ALLOCATE( ncomm_member(nn_ens_size)  )   ! communicator for members
      ALLOCATE( irank_member(jpproc), igrp_member(nn_ens_size)  ) ! temporary array of domain rank for a specific member

      ! Create global group
      CALL mpi_comm_group( localComm, igrp_global, ierr )

      ! Create one communicator for each ensemble member
      DO jmember = 1, nn_ens_size
         irank_member(:) = (/ (jjproc, jjproc=0,jpproc-1) /)
         irank_member(:) = irank_member(:) + ( jmember - 1 ) * jpproc
         ! Create the group for current member from the global group
         CALL mpi_group_incl( igrp_global, jpproc, irank_member, igrp_member(jmember), ierr )
         ! Create communicator for jmember
         CALL mpi_comm_create( localComm, igrp_member(jmember), ncomm_member(jmember), ierr )
      ENDDO

      ! Deallocate arrays
      DEALLOCATE( irank_member , igrp_member )

      ! Get rank of processor in global communicator
      ! and identify to which member communicator it belongs
      CALL mpi_comm_rank( localComm, impprank, ierr )
      imember = 1 + impprank / jpproc

      ! Set current member communicator
      memberComm = ncomm_member(imember)

      ! Return index of current ensemble member
      kmember = nn_ens_start + imember - 1

      ! Define communicators for ensemble diagnostics
      IF (ln_ens_diag) THEN
         ALLOCATE( ncomm_area(jpproc)  )     ! communicator for members
         ALLOCATE( irank_area(nn_ens_size),  igrp_area(jpproc) )
         ! Create one communicator for every subdomain
         DO jarea = 1, jpproc
             irank_area(:) = (/ (jjproc, jjproc=jarea-1,nn_ens_size*jpproc-1,jpproc) /)
             ! Create the group for current subdomain from the global group
             CALL mpi_group_incl( igrp_global, nn_ens_size, irank_area, igrp_area(jarea), ierr )
             ! Create the communicator for current subdomain
             CALL mpi_comm_create( localComm, igrp_area(jarea), ncomm_area(jarea), ierr )
         ENDDO
         ! Deallocate arrays
         DEALLOCATE( irank_area , igrp_area )
      ENDIF

   END SUBROUTINE mpp_set_ensemble

END MODULE mppens
