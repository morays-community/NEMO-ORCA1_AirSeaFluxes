MODULE stomod
   !!======================================================================
   !!                       ***  MODULE stomod  ***
   !!
   !! Implement stochastic parameterizations in the model
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   sto_mod          : apply parameterizations at each model time step
   !!   sto_mod_init     : initialize stochastic parameterizations
   !!   sto_mod_finalize : finalize stochastic parameterizations
   !!----------------------------------------------------------------------

   USE storng_check
   USE stoarray
   USE stopar
   USE stoeos
   USE stosppt
   USE stospp
   USE stoskeb
   USE stolug
   USE stotemplate
   USE stoexternal
   USE storst
   USE par_oce
   USE dom_oce
   USE in_out_manager
   USE lib_mpp
   USE pyfld

   IMPLICIT NONE
   PRIVATE

   ! maximum total number of stochastic arrays that can be used in NEMO
   ! (to increase if needed, it could be put as parameter in namelist)
   INTEGER :: jpstomax_nemo = 100

   ! Public routines
   PUBLIC sto_mod, sto_mod_init, sto_mod_finalize

CONTAINS

   SUBROUTINE sto_mod(kt)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod  ***
      !!
      !! Purpose : apply parameterizations at each model time step
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      IF (jpsto>0) THEN  ! There are active stochastic parameterizations

         ! Update stochastic fields
         CALL sto_par(kt)

         ! Apply every stochastic parameterization
         IF (ln_sto_eos)  CALL sto_eos(kt)
         IF (ln_sto_sppt) CALL sto_sppt(kt)
         IF (ln_sto_spp)  CALL sto_spp(kt)
         IF (ln_sto_skeb) CALL sto_skeb(kt)
         IF (ln_sto_lug)  CALL sto_lug(kt)

      ENDIF

   END SUBROUTINE sto_mod


   SUBROUTINE sto_mod_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod_init  ***
      !!
      !! Purpose : initialize stochastic parameterizations
      !!
      !!----------------------------------------------------------------------

      ! Stop if AGRIF if used
      IF ( lk_agrif ) CALL ctl_stop('EOS stochastic parametrization is not compatible with AGRIF')

      ! Initialize grid and mask
      CALL initialize_grid
      CALL initialize_mask

      ! Read general parameters from namelist
      CALL sto_param_init

      ! Perform test of random number generators
      IF (ln_rng_test) CALL sto_rng_test

      ! Request maximum number of stochastic arrays
      CALL sto_array_request_size(jpstomax_nemo)

      ! Initialization of the various stochastic schemes
      ! (including requests for stochastic arrays using sto_array_request_new)
      ! There must be one such routine for each stochastic scheme.
      IF (ln_sto_eos)  CALL sto_eos_init
      IF (ln_sto_sppt) CALL sto_sppt_init
      IF (ln_sto_spp)  CALL sto_spp_init
      IF (ln_sto_skeb) CALL sto_skeb_init
      IF (ln_sto_lug)  CALL sto_lug_init

      ! Custom module for W25 ML model
      CALL sto_w25_init()

      IF (jpsto>0) THEN  ! There are active stochastic parameterizations

         ! Initialize stochastic arrays
         CALL sto_array_init

         ! Initialize time iteration of stochastic arrays
         CALL sto_par_init

      ENDIF

   END SUBROUTINE sto_mod_init


   SUBROUTINE sto_mod_finalize( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_mod_finalize  ***
      !!
      !! Purpose : finalize stochastic parameterizations
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kt
      INTEGER :: jsto, jidx

      IF (jpsto>0) THEN  ! There are active stochastic parameterizations

         ! write final restart file
         CALL sto_rst_write( kt )

      ENDIF

   END SUBROUTINE sto_mod_finalize

   !!======================================================================
END MODULE stomod
