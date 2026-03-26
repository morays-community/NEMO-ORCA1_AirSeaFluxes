MODULE stospp
   !!======================================================================
   !!                       ***  MODULE stospp  ***
   !!
   !! Purpose : Implement SPP stochastic parameterization
   !!======================================================================
   USE stoarray
   USE par_oce
   USE dom_oce
   USE sbc_oce
   USE in_out_manager
   USE iom
   USE fldread        ! read input fields
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   ! Switches to activate SPP scheme
   LOGICAL, PUBLIC, SAVE :: ln_spp = .FALSE.

   ! SPP switches (with default value)
   INTEGER, PUBLIC, SAVE :: nn_spp_bfr = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_dqdt = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_dedt = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_avt = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_avm = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_qsi0 = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_relw = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_arnf = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_geot = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_aevd = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahubbl = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahvbbl = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tkelc = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tkedf = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tkeds = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tkebb = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tkefr = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahtu = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahtv = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahtw = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahtt = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahm1 = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_ahm2 = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_tdmp = 0
   ! SPP Sea-ice switches
   INTEGER, PUBLIC, SAVE :: nn_spp_icealb = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_icestr = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_blkd = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_blkh = 0
   INTEGER, PUBLIC, SAVE :: nn_spp_blke = 0

   ! SPP parameters
   REAL(wp), PUBLIC, SAVE :: rn_bfr_sd, rn_dqdt_sd, rn_dedt_sd, rn_avt_sd
   REAL(wp), PUBLIC, SAVE :: rn_avm_sd, rn_qsi0_sd, rn_relw_sd, rn_arnf_sd
   REAL(wp), PUBLIC, SAVE :: rn_geot_sd, rn_aevd_sd, rn_ahubbl_sd, rn_ahvbbl_sd
   REAL(wp), PUBLIC, SAVE :: rn_tkelc_sd,rn_tkedf_sd,rn_tkeds_sd,rn_tkebb_sd,rn_tkefr_sd
   REAL(wp), PUBLIC, SAVE :: rn_ahtu_sd, rn_ahtv_sd, rn_ahtw_sd, rn_ahtt_sd
   REAL(wp), PUBLIC, SAVE :: rn_ahm1_sd, rn_ahm2_sd, rn_tdmp_sd
   REAL(wp), PUBLIC, SAVE :: rn_icestr_sd, rn_icealb_sd
   REAL(wp), PUBLIC, SAVE :: rn_blkd_sd, rn_blkh_sd, rn_blke_sd

   ! SPP Options
   REAL(wp), SAVE :: rn_spp_tau, rn_spp_stdev
   INTEGER :: spp_filter_pass

   ! SPP stochastic array indices
   INTEGER, PUBLIC :: jsto_spp_alb
   INTEGER, PUBLIC :: jsto_spp_rhg
   INTEGER, PUBLIC :: jsto_spp_relw
   INTEGER, PUBLIC :: jsto_spp_dqdt
   INTEGER, PUBLIC :: jsto_spp_deds
   INTEGER, PUBLIC :: jsto_spp_arnf
   INTEGER, PUBLIC :: jsto_spp_geot
   INTEGER, PUBLIC :: jsto_spp_qsi0
   INTEGER, PUBLIC :: jsto_spp_bfr
   INTEGER, PUBLIC :: jsto_spp_aevd
   INTEGER, PUBLIC :: jsto_spp_avt
   INTEGER, PUBLIC :: jsto_spp_avm
   INTEGER, PUBLIC :: jsto_spp_tkelc
   INTEGER, PUBLIC :: jsto_spp_tkedf
   INTEGER, PUBLIC :: jsto_spp_tkeds
   INTEGER, PUBLIC :: jsto_spp_tkebb
   INTEGER, PUBLIC :: jsto_spp_tkefr

   INTEGER, PUBLIC :: jsto_spp_ahtu
   INTEGER, PUBLIC :: jsto_spp_ahtv
   INTEGER, PUBLIC :: jsto_spp_ahtw
   INTEGER, PUBLIC :: jsto_spp_ahtt

   INTEGER, PUBLIC :: jsto_spp_ahubbl
   INTEGER, PUBLIC :: jsto_spp_ahvbbl

   INTEGER, PUBLIC :: jsto_spp_ahm1
   INTEGER, PUBLIC :: jsto_spp_ahm2
   INTEGER, PUBLIC :: jsto_spp_ahm3
   INTEGER, PUBLIC :: jsto_spp_ahm4

   INTEGER, PUBLIC :: jsto_spp_blkd
   INTEGER, PUBLIC :: jsto_spp_blkh
   INTEGER, PUBLIC :: jsto_spp_blke

   INTEGER, PUBLIC :: jsto_spp_tdmp

   ! Public routines
   PUBLIC sto_spp, sto_spp_init

   !! * Substitutions
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"

CONTAINS

   SUBROUTINE sto_spp(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_spp  ***
      !!
      !! This routine is called at every time step
      !! to make appropriate use of the stochastic field
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

   END SUBROUTINE sto_spp


   SUBROUTINE sto_spp_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_spp_init  ***
      !!
      !! This routine is calle at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios
      NAMELIST/namsto_spp/&
      spp_filter_pass,rn_spp_stdev,rn_spp_tau,&
      rn_icestr_sd, rn_icealb_sd, &
      nn_spp_bfr,nn_spp_dqdt,nn_spp_dedt,nn_spp_avt,nn_spp_avm,nn_spp_qsi0,&
      nn_spp_ahtu,nn_spp_ahtv,nn_spp_ahm1,nn_spp_ahm2,&
      nn_spp_blkd,nn_spp_blkh,nn_spp_blke,&
      nn_spp_ahtw,nn_spp_ahtt,nn_spp_relw,nn_spp_arnf,nn_spp_geot,nn_spp_aevd,nn_spp_ahubbl,nn_spp_ahvbbl,&
      nn_spp_tkelc,nn_spp_tkedf,nn_spp_tkeds,nn_spp_tkebb,nn_spp_tkefr,&
      rn_bfr_sd,rn_dqdt_sd,rn_dedt_sd,rn_avt_sd,rn_avm_sd,rn_qsi0_sd,rn_ahtu_sd,rn_ahtv_sd,&
      rn_ahm1_sd,rn_ahm2_sd,&
      rn_ahtw_sd,rn_ahtt_sd, rn_relw_sd, rn_arnf_sd,rn_geot_sd, rn_aevd_sd,rn_ahubbl_sd,rn_ahvbbl_sd,&
      rn_tkelc_sd,rn_tkedf_sd,rn_tkeds_sd,rn_tkebb_sd,rn_tkefr_sd,&
      rn_blkd_sd,rn_blkh_sd,rn_blke_sd,&
      nn_spp_tdmp,rn_tdmp_sd

      ! Read namelist block corresponding to this stochastic scheme
      READ_NML_REF(numnam_sto,namsto_spp)
      READ_NML_CFG(numnam_sto,namsto_spp)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stospp : Stochastic SPP scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto_spp: '
         WRITE(numout,*)
      ENDIF

      ! Define characteristics of each stochastic array
      IF( nn_spp_bfr == 1 ) THEN
         ! Request index for a each stochastic array
         CALL sto_array_request_new(jsto_spp_bfr)

         ! Set features of the corresponding stochastic field
         stofields(jsto_spp_bfr)%type_t='arn'
         stofields(jsto_spp_bfr)%corr_t=rn_spp_tau * 86400._wp /rdt
         stofields(jsto_spp_bfr)%nar_order=1
         stofields(jsto_spp_bfr)%nar_update=1

         stofields(jsto_spp_bfr)%type_xy='diffusive'
         stofields(jsto_spp_bfr)%diff_passes=spp_filter_pass
         stofields(jsto_spp_bfr)%diff_type=1

         stofields(jsto_spp_bfr)%std=rn_spp_stdev
      ENDIF

      ! -> use stofields(jsto_spp_bfr)%sto2d(:,:) to perturb bfr...
      ! ... and repeat the same block to get a stichastic array
      !     for each parameter to perturb....

   END SUBROUTINE sto_spp_init


!   SUBROUTINE spp_aht ( kt, coeff,nn_type, rn_sd, kspp  )
   SUBROUTINE spp_aht
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_aht ***
      !!
      !! ** Purpose :   Perturbing diffusivity parameters
      !!                As spp_gen, but specifically designed for diffusivity
      !!                where coeff can be 2D or 3D
      !!----------------------------------------------------------------------

      ! insert updated routine here
   END SUBROUTINE


!   SUBROUTINE spp_ahm ( kt, coeff,nn_type, rn_sd, kspp  )
   SUBROUTINE spp_ahm
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_ahm ***
      !!
      !! ** Purpose :   Perturbing viscosity parameters
      !!                As spp_gen, but specifically designed for viscosity
      !!                where coeff can be 2D or 3D
      !!----------------------------------------------------------------------

      ! insert updated routine here
   END SUBROUTINE

   !!======================================================================
END MODULE stospp
