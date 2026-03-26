MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE pyfld  ***
   !! Python module :   variables defined in core memory
   !!======================================================================
   !! History :  4.2  ! 2025-11  (A. Barge)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   pyfld_alloc : allocation of fields arrays for Python coupling module (pycpl)
   !!----------------------------------------------------------------------
   !!=====================================================
   USE oce            ! ocean fields
   USE dom_oce        ! ocean metrics fields
   USE par_oce        ! ocean parameters
   USE lib_mpp        ! MPP library
   USE pycpl          ! Python coupling module
   USE iom

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!                    2D Python coupling Module fields
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: ext_qs, ext_ql, ext_taux, ext_tauy !: fields to store 2D Python fields

   !!----------------------------------------------------------------------
   !!                    3D Python coupling Module fields
   !!----------------------------------------------------------------------

#  include "do_loop_substitute.h90"

CONTAINS

   SUBROUTINE init_python_fields()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE init_python_fields  ***
      !!
      !! ** Purpose :   Initialisation of the Python module
      !!
      !! ** Method  :   * Allocate arrays for Python fields
      !!                * Configure Python coupling
      !!----------------------------------------------------------------------
      !
      ! Allocate fields
      ALLOCATE( ext_taux(A2D(0)), ext_tauy(A2D(0)), ext_qs(A2D(0)), ext_ql(A2D(0)) )
      !
      ! configure coupling
      CALL init_python_coupling()
      !
   END SUBROUTINE init_python_fields


   SUBROUTINE finalize_python_fields()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE finalize_python_fields  ***
      !!
      !! ** Purpose :   Free memory used by Python module
      !!
      !! ** Method  :   * deallocate arrays for Python fields
      !!                * deallocate Python coupling
      !!----------------------------------------------------------------------
      !
      ! Free memory
      DEALLOCATE( ext_taux, ext_tauy, ext_qs, ext_ql )
      !
      ! terminate coupling environment
      CALL finalize_python_coupling()
      !
   END SUBROUTINE finalize_python_fields


   SUBROUTINE inferences( kt, wndx, wndy, tair, sst, hum, slp )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE inferences  ***
      !!
      !! ** Purpose :   update the ocean data with the ML based models
      !!
      !! ** Method  :   *
      !!                *
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt               ! ocean time step
      REAL(wp), DIMENSION(A2D(0)), INTENT(in) :: wndx, wndy, tair, hum, sst, slp ! surface fields
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('Python')
      !
      ! Send inputs
      CALL send_to_python('ux', wndx, kt)
      CALL send_to_python('uy', wndy, kt)
      CALL send_to_python('tair', tair, kt)
      CALL send_to_python('toce', sst, kt)
      CALL send_to_python('q', hum, kt)
      CALL send_to_python('p', slp, kt)
      !
      ! ... external operations ...
      !
      ! Receive results
      CALL receive_from_python('Ql', ext_ql, kt)
      CALL receive_from_python('Qs', ext_qs, kt)
      CALL receive_from_python('taux', ext_taux, kt)
      CALL receive_from_python('tauy', ext_tauy, kt)
      !
      ! Write returned results
      CALL iom_put( "ext_ql" , ext_ql )
      CALL iom_put( "ext_qs" , ext_qs )
      CALL iom_put( "ext_taux" , ext_taux )
      CALL iom_put( "ext_tauy" , ext_tauy )
      !
      IF( ln_timing )   CALL timing_stop('Python')
      !
   END SUBROUTINE inferences

END MODULE pyfld
