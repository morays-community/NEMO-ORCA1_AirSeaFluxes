MODULE stosppt
   !!======================================================================
   !!                       ***  MODULE stosppt  ***
   !!
   !! Purpose : Implement SPPT stochastic parameterization
   !!=====================================================================
   !! History :  3.4  ! 2018-05 (A. Storto)  Original code
   !! History :  5.x  ! 2025-01 (A. Storto)  Updated code for STO
   !!----------------------------------------------------------------------
   !!
   !! The SPPT implementation perturbs tendencies collinearly to the 
   !! unperturbed tendencies, for selected processes that the user can 
   !! switch on.
   !! More details in https://doi.org/10.1002/qj.3990
   !! 
   !!----------------------------------------------------------------------
   !!   
   !! List of public routines:
   !!   sto_sppt, sto_sppt_init, sppt_collect, tra_sppt, 
   !!   dyn_sppt, tau_sppt, sppt_ice
   !!----------------------------------------------------------------------
   !
   USE par_oce
   USE dom_oce
   USE trd_oce
   USE sbc_oce
   USE stoarray
   USE in_out_manager
   USE iom
   USE lib_mpp
#if defined key_si3
   USE par_ice
   USE ice
#endif

   !! * Substitutions
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"

   IMPLICIT NONE
   PRIVATE

   ! Publicize variables from stoarray
   PUBLIC :: ln_sto_sppt
   !
   ! SPPT Logical switches for individual tendencies
   LOGICAL, PUBLIC, SAVE :: ln_sppt_traldf   = .FALSE.  ! Switch for lateral diffusion
   LOGICAL, PUBLIC, SAVE :: ln_sppt_trazdf   = .FALSE.  ! Switch for vertical diffusion
   LOGICAL, PUBLIC, SAVE :: ln_sppt_trazdfp  = .FALSE.  ! Switch for pure vertical diffusion
   LOGICAL, PUBLIC, SAVE :: ln_sppt_traevd   = .FALSE.  ! Switch for enhanced vertical diffusion
   LOGICAL, PUBLIC, SAVE :: ln_sppt_trabbc   = .FALSE.  ! Switch for bottom boundary condition
   LOGICAL, PUBLIC, SAVE :: ln_sppt_trabbl   = .FALSE.  ! Switch for bottom boundary layer
   LOGICAL, PUBLIC, SAVE :: ln_sppt_tranpc   = .FALSE.  ! Switch for non-penetrative convection
   LOGICAL, PUBLIC, SAVE :: ln_sppt_tradmp   = .FALSE.  ! Switch for tracer damping
   LOGICAL, PUBLIC, SAVE :: ln_sppt_traqsr   = .FALSE.  ! Switch for solar radiation
   LOGICAL, PUBLIC, SAVE :: ln_sppt_transr   = .FALSE.  ! Switch for non-solar radiation / freshwater flux
   !
   LOGICAL, PUBLIC, SAVE :: ln_sppt_taumap   = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_distcoast     = .FALSE.
   !
   LOGICAL, PUBLIC, SAVE :: ln_sppt_dynldf   = .FALSE.  ! Switch for lateral viscosity
   LOGICAL, PUBLIC, SAVE :: ln_sppt_dynzdf   = .FALSE.  ! Switch for vertical viscosity
   LOGICAL, PUBLIC, SAVE :: ln_sppt_dynbfr   = .FALSE.  ! Switch for bottom friction
   !
   LOGICAL, PUBLIC, SAVE :: ln_sppt_glocon   = .FALSE.
   !
   LOGICAL, PUBLIC, SAVE :: ln_sppt_tau      = .FALSE.  ! Switch for wind stress
   LOGICAL, PUBLIC, SAVE :: ln_sppt_icehdf   = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_icelat   = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_icezdf   = .FALSE.
   !
   ! Global switches activating SPPT schemes
   LOGICAL, PUBLIC, SAVE :: ln_sppt_tra = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_dyn = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_ice = .FALSE.
   !
   ! SPPT Options (with default values)
   INTEGER, PARAMETER :: nn_sppt_ord = 1
   INTEGER, PARAMETER, PUBLIC :: ktra = 1, kdyn = 2
   INTEGER :: sppt_filter_pass = 30
   INTEGER :: sppt_step
   INTEGER :: nn_vwei = 0
   INTEGER :: nn_sppt_freq = 1
   REAL(wp), SAVE :: rn_sppt_tau = 5._wp
   LOGICAL,  SAVE :: ln_sppt_bound = .false.
   REAL(wp), SAVE :: rn_sppt_bound = 3._wp
   REAL(wp), SAVE :: rn_distcoast = 10000._wp
   REAL(wp), SAVE :: rn_sppt_stdev = 1._wp
   REAL(wp), SAVE :: rn_uv_infl = 1._wp
   REAL(wp), SAVE :: rn_ice_infl = 1._wp
   REAL(wp), ALLOCATABLE :: zicewrk(:,:,:), sppt_mask(:,:,:), sppt_mask2(:,:)
   REAL(wp), ALLOCATABLE :: sppt_tra(:,:,:,:), sppt_dyn(:,:,:,:), zwper(:,:,:)
   !
   ! Indices of SPPT stochastic arrays
   INTEGER, PUBLIC :: jsto_sppt
   INTEGER, PUBLIC :: jsto_sppt_ice
   INTEGER, PUBLIC :: jsto_sppt_tau
   !
   ! Public routines
   PUBLIC tra_sppt, dyn_sppt, sto_sppt, sto_sppt_init, tau_sppt, sppt_collect, sppt_ice
   LOGICAL, PUBLIC, ALLOCATABLE :: ln_sppt_tra_ten(:), ln_sppt_dyn_ten(:)
   LOGICAL, SAVE, ALLOCATABLE :: ln_sppt_tra_ten_save(:), ln_sppt_dyn_ten_save(:)
   !
   INTERFACE sppt_glocon
          MODULE PROCEDURE sppt_glocon_3d, sppt_glocon_2d
   END INTERFACE sppt_glocon
   !
CONTAINS

   SUBROUTINE sto_sppt(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_sppt  ***
      !!
      !! This routine is called at the beginning of every timestep
      !! to setup tendency calculation
      !!
      !!----------------------------------------------------------------------
      INTEGER  :: kt
      !
      ! Output stochastic field
      IF (ln_sppt_tra .OR. ln_sppt_dyn) &
                     & CALL iom_put( 'sppt_sto',     stofields(jsto_sppt)%sto2d(:,:) )
      IF (ln_sppt_ice) CALL iom_put( 'sppt_ice_sto', stofields(jsto_sppt_ice)%sto2d(:,:) )
      IF (ln_sppt_tau) CALL iom_put( 'sppt_tau_sto', stofields(jsto_sppt_tau)%sto2d(:,:) )
      !
      IF( MOD( kt - 1, nn_sppt_freq ) == 0 .OR. kt == nit000 ) THEN
               ln_sppt_tra_ten(:) = ln_sppt_tra_ten_save(:)
               ln_sppt_dyn_ten(:) = ln_sppt_dyn_ten_save(:)
      ELSE
               ln_sppt_tra_ten(:) = .false.
               ln_sppt_dyn_ten(:) = .false.
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE sto_sppt_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_sppt_init  ***
      !!
      !! This routine is called at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------
      INTEGER  :: nn_sppt_tra, nn_sppt_dyn, nn_spp_aht, nn_sppt_ice
      INTEGER  :: ios, inum, ji, jj, nbl
      REAL(wp), ALLOCATABLE :: zdc(:,:,:), wei(:)
      REAL(wp) :: w1, w2
      !!
      NAMELIST/namsto_sppt/ ln_sppt_taumap, rn_sppt_tau, &
      sppt_filter_pass, ln_sppt_bound, rn_sppt_bound, sppt_step, rn_sppt_stdev, &
      ln_distcoast, rn_distcoast, nn_vwei, nn_sppt_freq, &
      ln_sppt_glocon, rn_uv_infl, rn_ice_infl, ln_sppt_traldf, &
      ln_sppt_trazdf, ln_sppt_trazdfp, ln_sppt_traevd, ln_sppt_trabbc, ln_sppt_trabbl, &
      ln_sppt_tranpc, ln_sppt_tradmp, ln_sppt_traqsr, ln_sppt_transr, &
      ln_sppt_dynldf, ln_sppt_dynzdf, ln_sppt_dynbfr, &
      ln_sppt_icehdf, ln_sppt_icelat, ln_sppt_icezdf, ln_sppt_tau, &
      ln_sppt_tra, ln_sppt_dyn, ln_sppt_ice

      ! Read namelist block corresponding to this stochastic scheme
      READ_NML_REF(numnam_sto,namsto_sppt)
      READ_NML_CFG(numnam_sto,namsto_sppt)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stosppt : Stochastic SPPT scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto_sppt: '
         WRITE(numout,*)
         WRITE(numout,*) '       SPPT step (0: beginning; 1: before ZVD; 2: after ZVD)  sppt_step :', sppt_step
         WRITE(numout,*) '       Use map of decorr. time scale                     ln_sppt_taumap :', ln_sppt_taumap
         WRITE(numout,*) '       If ln_sppt_taumap FALSE, use this constant (in days) rn_sppt_tau :', rn_sppt_tau
         WRITE(numout,*) '       Number of filter passes to correlate in space   sppt_filter_pass :', sppt_filter_pass
         WRITE(numout,*) '       Standard deviation of the white noise              rn_sppt_stdev :', rn_sppt_stdev
         WRITE(numout,*) '       Apply global conservation constraints             ln_sppt_glocon :', ln_sppt_glocon
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on tracers:'
         WRITE(numout,*) '       Switch for lateral diffusion                      ln_sppt_traldf :', ln_sppt_traldf
         WRITE(numout,*) '       Switch for vertical diffusion                     ln_sppt_trazdf :', ln_sppt_trazdf
         WRITE(numout,*) '       Switch for pure vertical diffusion               ln_sppt_trazdfp :', ln_sppt_trazdfp
         WRITE(numout,*) '       Switch for enhanced vertical diffusion            ln_sppt_traevd :', ln_sppt_traevd
         WRITE(numout,*) '       Switch for bottom boundary condition              ln_sppt_trabbc :', ln_sppt_trabbc
         WRITE(numout,*) '       Switch for bottom boundary layer                  ln_sppt_trabbl :', ln_sppt_trabbl
         WRITE(numout,*) '       Switch for non-penetrative convection             ln_sppt_tranpc :', ln_sppt_tranpc
         WRITE(numout,*) '       Switch for tracer damping                         ln_sppt_tradmp :', ln_sppt_tradmp
         WRITE(numout,*) '       Switch for solar radiation                        ln_sppt_traqsr :', ln_sppt_traqsr
         WRITE(numout,*) '       Switch for non-solar rad. / freshwater flx flux   ln_sppt_transr :', ln_sppt_transr
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on dynamics:'
         WRITE(numout,*) '       Inflation coefficient for SPPT on DYN             rn_uv_infl     :', rn_uv_infl
         WRITE(numout,*) '       Switch for lateral diffusion                      ln_sppt_dynldf :', ln_sppt_dynldf
         WRITE(numout,*) '       Switch for vertical diffusion                     ln_sppt_dynzdf :', ln_sppt_dynzdf
         WRITE(numout,*) '       Switch for bottom friction                        ln_sppt_dynbfr :', ln_sppt_dynbfr
         WRITE(numout,*)
         WRITE(numout,*) '       Switch for tau perturbation                       ln_sppt_tau    :', ln_sppt_tau
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on sea-ice:'
         WRITE(numout,*) '       Inflation coefficient for SPPT on ICE             rn_ice_infl    :', rn_ice_infl
         WRITE(numout,*) '       Switch for sea-ice diffusivity                    ln_sppt_icehdf :', ln_sppt_icehdf
         WRITE(numout,*) '       Switch for sea-ice lateral accretion              ln_sppt_icelat :', ln_sppt_icelat
         WRITE(numout,*) '       Switch for sea-ice vertical thermodyn.            ln_sppt_icezdf :', ln_sppt_icezdf
         WRITE(numout,*)
         WRITE(numout,*) '       Flag for bounding perturbations                    ln_sppt_bound :', ln_sppt_bound
         WRITE(numout,*) '       Bound for perturbations                            rn_sppt_bound :', rn_sppt_bound
         WRITE(numout,*) '       Smoothing of perturbation close to coast            ln_distcoast :', ln_distcoast
         WRITE(numout,*) '       Spatial scale of the smoothing near coasts (m)      rn_distcoast :', rn_distcoast
         WRITE(numout,*) '       Type of vertical weight:                                 nn_vwei :', nn_vwei
         WRITE(numout,*) '                            0 : No weight '
         WRITE(numout,*) '                            1 : Top/Bottom smoothing'
         WRITE(numout,*) '                            2 : Bottom smoothing'
      ENDIF

      ! Check which SPPT schemes are activated
      nn_sppt_tra = COUNT( (/ln_sppt_traldf, ln_sppt_trazdf, &
                             ln_sppt_trazdfp, ln_sppt_traevd, ln_sppt_trabbc, &
                             ln_sppt_trabbl, ln_sppt_tranpc, ln_sppt_tradmp, &
                             ln_sppt_traqsr, ln_sppt_transr/) )
      nn_sppt_dyn = COUNT( (/ln_sppt_dynldf, ln_sppt_dynzdf, ln_sppt_dynbfr/) )
      nn_sppt_ice = COUNT( (/ln_sppt_icehdf, ln_sppt_icelat, ln_sppt_icezdf/) )

      ln_sppt_tra = ( nn_sppt_tra > 0 )
      ln_sppt_dyn = ( nn_sppt_dyn > 0 )
      ln_sppt_ice = ( nn_sppt_ice > 0 )

      ! Define characteristics of required stochastic arrays 
      ! - we assume the same random fields for dyn and tra
      ! - we use a diffrent one (in 2D) for the ice
      ! - we use a diffrent one (in 2D) for the tau
      IF( ln_sppt_tra .OR. ln_sppt_dyn ) THEN
         ! Request index for a each stochastic array
         CALL sto_array_request_new(jsto_sppt)

         ! Set features of the corresponding stochastic field
         stofields(jsto_sppt)%dim=2

         stofields(jsto_sppt)%type_t='arn'
         stofields(jsto_sppt)%corr_t=rn_sppt_tau * 86400. / rdt
         stofields(jsto_sppt)%nar_order=nn_sppt_ord
         stofields(jsto_sppt)%nar_update=nn_sppt_freq

         stofields(jsto_sppt)%type_xy='diffusive'
         stofields(jsto_sppt)%diff_passes=sppt_filter_pass
         stofields(jsto_sppt)%diff_type=1

         stofields(jsto_sppt)%std=rn_sppt_stdev
      ENDIF

      IF( ln_sppt_ice ) THEN 
         ! Request index for a each stochastic array
         CALL sto_array_request_new(jsto_sppt_ice)

         ! Set features of the corresponding stochastic field
         stofields(jsto_sppt_ice)%dim=2

         stofields(jsto_sppt_ice)%type_t='arn'
         stofields(jsto_sppt_ice)%corr_t=rn_sppt_tau * 86400. / rdt
         stofields(jsto_sppt_ice)%nar_order=nn_sppt_ord
         stofields(jsto_sppt_ice)%nar_update=nn_sppt_freq

         stofields(jsto_sppt_ice)%type_xy='diffusive'
         stofields(jsto_sppt_ice)%diff_passes=sppt_filter_pass
         stofields(jsto_sppt_ice)%diff_type=1

         stofields(jsto_sppt_ice)%std=rn_sppt_stdev * rn_ice_infl
      ENDIF

      IF( ln_sppt_tau ) THEN
         ! Request index for a each stochastic array
         CALL sto_array_request_new(jsto_sppt_tau)

         ! Set features of the corresponding stochastic field
         stofields(jsto_sppt_tau)%dim=2

         stofields(jsto_sppt_tau)%type_t='arn'
         stofields(jsto_sppt_tau)%corr_t=rn_sppt_tau * 86400. / rdt
         stofields(jsto_sppt_tau)%nar_order=nn_sppt_ord
         stofields(jsto_sppt_tau)%nar_update=nn_sppt_freq

         !stofields(jsto_sppt_tau)%type_xy='diffusive'
         !stofields(jsto_sppt_tau)%diff_passes=sppt_filter_pass
         !stofields(jsto_sppt_tau)%diff_type=1
         stofields(jsto_sppt_tau)%type_xy='kernel'
         stofields(jsto_sppt_tau)%corr_xy=20.
         stofields(jsto_sppt_tau)%ker_coord=2
         stofields(jsto_sppt_tau)%ker_density=0.0

         stofields(jsto_sppt_tau)%std=rn_sppt_stdev
      ENDIF

      ALLOCATE ( sppt_mask(jpi,jpj,jpk) )
      ALLOCATE ( sppt_mask2(jpi,jpj) )

      IF( ln_distcoast ) THEN
              ALLOCATE( zdc(jpi,jpj,jpk) ); zdc=1._wp
              CALL iom_open('dist.coast', inum )
              CALL iom_get(inum,jpdom_auto,'Tcoast',zdc)
              CALL iom_close( inum )
              sppt_mask(:,:,:) = 1._wp - exp(-zdc(:,:,:)/rn_distcoast )
              DEALLOCATE ( zdc )
      ELSE
              sppt_mask(:,:,:) = 1._wp
      ENDIF
      sppt_mask2(:,:) = sppt_mask(:,:,1)

      IF( nn_vwei .gt. 0 ) THEN
         IF( nn_vwei .eq. 1 ) THEN
             w1 = 0._wp
             w2 =0.5_wp
         ELSEIF( nn_vwei .eq. 2 ) THEN
             w1 = 1._wp
             w2 = 1._wp
         ELSE
             CALL ctl_stop('STOP','SPPT: nn_vwei > 2')
         ENDIF
         !
         ALLOCATE (wei(jpk))
         !
         DO_2D( 0, 0, 0, 0)
               nbl = ht_0(ji,jj)
               wei(:) = 0._wp
               IF(nbl>1) THEN
                  wei(1:nbl) = 1._wp
                  wei(1) = w1
                  wei(2) = w2
                  wei(nbl-1) = 0.5_wp
                  wei(nbl) = 0._wp
                  sppt_mask(ji,jj,:) = sppt_mask(ji,jj,:) * wei(:)
               ENDIF
         END_2D
         !
         DEALLOCATE (wei)
         !
      ENDIF

      ALLOCATE ( ln_sppt_tra_ten(jptot_tra) )
      ALLOCATE ( ln_sppt_tra_ten_save(jptot_tra) )
      ln_sppt_tra_ten (:) = .false.
      IF( ln_sppt_traldf ) ln_sppt_tra_ten(jptra_ldf) = .true.
      IF( ln_sppt_trazdf ) ln_sppt_tra_ten(jptra_zdf) = .true.
      IF( ln_sppt_trazdfp) ln_sppt_tra_ten(jptra_zdfp)= .true.
      IF( ln_sppt_traevd ) ln_sppt_tra_ten(jptra_evd )= .true.
      IF( ln_sppt_trabbc ) ln_sppt_tra_ten(jptra_bbc )= .true. 
      IF( ln_sppt_trabbl ) ln_sppt_tra_ten(jptra_bbl )= .true.
      IF( ln_sppt_tranpc ) ln_sppt_tra_ten(jptra_npc )= .true.
      IF( ln_sppt_tradmp ) ln_sppt_tra_ten(jptra_dmp )= .true.
      IF( ln_sppt_traqsr ) ln_sppt_tra_ten(jptra_qsr )= .true.
      IF( ln_sppt_transr ) ln_sppt_tra_ten(jptra_nsr )= .true.

      ALLOCATE ( ln_sppt_dyn_ten(jptot_dyn) )
      ALLOCATE ( ln_sppt_dyn_ten_save(jptot_dyn) )
      ln_sppt_dyn_ten (:) = .false.
      IF( ln_sppt_dynldf ) ln_sppt_dyn_ten(jpdyn_ldf) = .true.
      IF( ln_sppt_dynzdf ) ln_sppt_dyn_ten(jpdyn_zdf) = .true.
      IF( ln_sppt_dynbfr ) ln_sppt_dyn_ten(jpdyn_bfr) = .true.

      ln_sppt_tra_ten_save(:) = ln_sppt_tra_ten(:) 
      ln_sppt_dyn_ten_save(:) = ln_sppt_dyn_ten(:) 

      ALLOCATE( sppt_tra(jpi,jpj,jpk,2), sppt_dyn(jpi,jpj,jpk,2), zwper(jpi,jpj,jpk) )
      !
   END SUBROUTINE

   SUBROUTINE sppt_collect ( ktype, ptrd1, ptrd2, ktrd, kt, Kmm)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sppt_apply ***
      !!
      !! ** Purpose :   Apply perturbation to tracers
      !!                Step 0/1 for tendencies, step 2 for variables
      !!                after timestepping
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: ktype, ktrd, kt, Kmm
      REAL(wp), INTENT(IN) :: ptrd1(:,:,:), ptrd2(:,:,:)
      !
      IF(ktype .eq. ktra) THEN
         sppt_tra(:,:,:,jp_tem) = sppt_tra(:,:,:,jp_tem) + ptrd1
         sppt_tra(:,:,:,jp_sal) = sppt_tra(:,:,:,jp_sal) + ptrd2
      ELSEIF(ktype .eq. kdyn) THEN
         sppt_dyn(:,:,:,1) = sppt_dyn(:,:,:,1) + ptrd1
         sppt_dyn(:,:,:,2) = sppt_dyn(:,:,:,2) + ptrd2
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE tra_sppt ( kstp, Nbb, Nnn, Nrhs, pts, Naa  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sppt_apply ***
      !!
      !! ** Purpose :   Apply perturbation to tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN)    :: kstp, Nbb, Nnn, Nrhs, Naa
      REAL(wp),INTENT(INOUT) :: pts(jpi,jpj,jpk,jpts,jpt)
      REAL(wp), ALLOCATABLE  :: zwr(:,:,:)
      INTEGER :: jk
      !
      DO jk = 1, jpk
         zwper(:,:,jk) = sppt_mask(:,:,jk) * stofields(jsto_sppt)%sto2d(:,:) * tmask(:,:,jk)
      ENDDO
      !
      IF( ln_sppt_bound ) THEN
          WHERE( zwper .lt. -ABS(rn_sppt_bound) ) zwper = -ABS(rn_sppt_bound)
          WHERE( zwper .gt.  ABS(rn_sppt_bound) ) zwper =  ABS(rn_sppt_bound)
      ENDIF
      !
      IF( ln_sppt_glocon ) THEN
          ALLOCATE ( zwr(jpi,jpj,jpk) )
          !
          zwr(:,:,:) = zwper(:,:,:) * sppt_tra(:,:,:,jp_tem)
          CALL sppt_glocon( zwr, 'T' )
          pts(:,:,:,jp_tem,Naa) = pts(:,:,:,jp_tem,Naa) + zwr(:,:,:)*rDt
          !
          zwr(:,:,:) = zwper(:,:,:) * sppt_tra(:,:,:,jp_sal)
          CALL sppt_glocon( zwr, 'T' )
          pts(:,:,:,jp_sal,Naa) = pts(:,:,:,jp_sal,Naa) + zwr(:,:,:)*rDt
          !
          DEALLOCATE ( zwr )
      ELSE
          pts(:,:,:,jp_tem,Naa) = pts(:,:,:,jp_tem,Naa) + sppt_tra(:,:,:,jp_tem)*zwper(:,:,:)*rDt
          pts(:,:,:,jp_sal,Naa) = pts(:,:,:,jp_sal,Naa) + sppt_tra(:,:,:,jp_sal)*zwper(:,:,:)*rDt
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE dyn_sppt ( kstp, Nbb, Nnn, Nrhs, pu, pv, Naa  ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_sppt ***
      !!
      !! ** Purpose :   Apply perturbation to momentum
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN)    :: kstp, Nbb, Nnn, Nrhs, Naa
      REAL(wp),INTENT(INOUT) :: pu(jpi,jpj,jpk,jpt), pv(jpi,jpj,jpk,jpt)
      REAL(wp), ALLOCATABLE  :: zwr(:,:,:)
      INTEGER :: jk
      !
      DO jk = 1, jpk
         zwper(:,:,jk) = sppt_mask(:,:,jk) * stofields(jsto_sppt)%sto2d(:,:) * umask(:,:,jk)*vmask(:,:,jk)
      ENDDO
      !
      IF( ln_sppt_bound ) THEN
          WHERE( zwper .gt.  ABS(rn_sppt_bound) ) zwper =  ABS(rn_sppt_bound)
      ENDIF
      !
      IF( ln_sppt_glocon ) THEN
          ALLOCATE ( zwr(jpi,jpj,jpk) )
          !
          zwr(:,:,:) = zwper(:,:,:) * sppt_dyn(:,:,:,1)
          CALL sppt_glocon( zwr, 'U' )
          pu(:,:,:,Naa) = pu(:,:,:,Naa) + zwr(:,:,:)*rDt
          !
          zwr(:,:,:) = zwper(:,:,:) * sppt_dyn(:,:,:,2)
          CALL sppt_glocon( zwr, 'V' )
          pv(:,:,:,Naa) = pv(:,:,:,Naa) + zwr(:,:,:)*rDt
          !
          DEALLOCATE ( zwr )
      ELSE
          pu(:,:,:,Naa) = pu(:,:,:,Naa) + sppt_dyn(:,:,:,1)*zwper(:,:,:)*rDt
          pv(:,:,:,Naa) = pv(:,:,:,Naa) + sppt_dyn(:,:,:,2)*zwper(:,:,:)*rDt
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE tau_sppt( kt, tau, tau_b, tav, tav_b )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tau_sppt ***
      !!
      !! ** Purpose :   Apply perturbation to wind stress
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: kt 
      REAL(wp), INTENT(INOUT), DIMENSION(jpi,jpj) :: tau, tau_b, tav, tav_b
      REAL(wp), ALLOCATABLE  :: zwr(:,:)
      !
      zwper(:,:,1) = sppt_mask2(:,:)*stofields(jsto_sppt_tau)%sto2d(:,:)*umask(:,:,1)*vmask(:,:,1)
      !
      IF( ln_sppt_bound ) THEN
          WHERE( zwper(:,:,1) .lt. -ABS(rn_sppt_bound) ) zwper(:,:,1) = -ABS(rn_sppt_bound)
          WHERE( zwper(:,:,1) .gt.  ABS(rn_sppt_bound) ) zwper(:,:,1) =  ABS(rn_sppt_bound)
      ENDIF
      !
      IF( ln_sppt_glocon ) THEN
          !
          ALLOCATE ( zwr(jpi,jpj) )
          !
          zwr(:,:) = zwper(:,:,1) * (tau(:,:)-tau_b(:,:))
          CALL sppt_glocon( zwr, 'U' )
          tau(:,:) = tau(:,:) + zwr(:,:)
          !
          zwr(:,:) = zwper(:,:,1) * (tav(:,:)-tav_b(:,:))
          CALL sppt_glocon( zwr, 'V' )
          tav(:,:) = tav(:,:) + zwr(:,:)
          !
          DEALLOCATE ( zwr )
      ELSE
          tau(:,:) = tau(:,:) + (tau(:,:)-tau_b(:,:))*zwper(:,:,1)
          tav(:,:) = tav(:,:) + (tav(:,:)-tav_b(:,:))*zwper(:,:,1)
      ENDIF
      !
   END SUBROUTINE

#if defined key_si3
   SUBROUTINE sppt_ice( kstep, kl  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_ice ***
      !!
      !! ** Purpose :   Apply collinear perturbation to ice fields
      !!                For specific processes coded in IC3
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: kstep         ! Start (1) or End (2) of ice routine
      INTEGER, INTENT(IN) :: kl            ! Category
      INTEGER :: jmt                       ! Number of sea-ice variables (depends on LIM version, process)
      REAL(wp), DIMENSION(jpi,jpj) :: zmul
      INTEGER :: jk
      !
      jmt=11+2*nlay_s+3*nlay_i
      !
      IF( kstep == 1 ) THEN        ! Store the before values
        !
        IF ( .NOT. ALLOCATED ( zicewrk) ) ALLOCATE ( zicewrk(jpi,jpj,jmt) )
        !
        zicewrk(:,:,1) = at_i(:,:)
        zicewrk(:,:,2) = a_i (:,:,kl)
        zicewrk(:,:,3) = h_i (:,:,kl)
        zicewrk(:,:,4) = h_s (:,:,kl)
        zicewrk(:,:,5) = t_su(:,:,kl)
        zicewrk(:,:,6) = s_i (:,:,kl)
        DO jk = 1, nlay_s
           zicewrk(:,:,6+jk) = t_s(:,:,jk,kl)
        ENDDO
        DO jk = 1, nlay_s
           zicewrk(:,:,6+jk+nlay_s) = e_s(:,:,jk,kl)
        ENDDO
        DO jk = 1, nlay_i
           zicewrk(:,:,6+jk+2*nlay_s) = t_i  (:,:,jk,kl)
        ENDDO
        DO jk = 1, nlay_i
           zicewrk(:,:,6+jk+2*nlay_s+nlay_i) = e_i  (:,:,jk,kl)
        ENDDO
        DO jk = 1, nlay_i
           zicewrk(:,:,6+jk+2*nlay_s+2*nlay_i) = sz_i  (:,:,jk,kl)
        ENDDO
        zicewrk(:,:,6+2*nlay_s+3*nlay_i+1) = sst_m(:,:)
        zicewrk(:,:,6+2*nlay_s+3*nlay_i+2) = sss_m(:,:)
        zicewrk(:,:,6+2*nlay_s+3*nlay_i+3) = frq_m(:,:)
        zicewrk(:,:,6+2*nlay_s+3*nlay_i+4) = o_i (:,:,kl)
        zicewrk(:,:,6+2*nlay_s+3*nlay_i+5) = oa_i(:,:,kl)
        !
      ELSEIF ( kstep == 2 ) THEN
        !
        zmul = stofields(jsto_sppt_ice)%sto2d(:,:)*sppt_mask2(:,:)
        !
        at_i(:,:)     = at_i(:,:)   + ( at_i(:,:)   - zicewrk(:,:,1) ) * zmul(:,:)
        a_i(:,:,kl)   = a_i(:,:,kl) + ( a_i(:,:,kl) - zicewrk(:,:,2) ) * zmul(:,:)
        h_i(:,:,kl)   = h_i(:,:,kl) + ( h_i(:,:,kl) - zicewrk(:,:,3) ) * zmul(:,:)
        h_s(:,:,kl)   = h_s(:,:,kl) + ( h_s(:,:,kl) - zicewrk(:,:,4) ) * zmul(:,:)
        t_su(:,:,kl)  = t_su(:,:,kl)+ ( t_su(:,:,kl)- zicewrk(:,:,5) ) * zmul(:,:)
        s_i(:,:,kl)   = s_i(:,:,kl) + ( s_i(:,:,kl) - zicewrk(:,:,6) ) * zmul(:,:)
        DO jk = 1, nlay_s
           t_s(:,:,jk,kl) = t_s(:,:,jk,kl) + (t_s(:,:,jk,kl)-zicewrk(:,:,6+jk)) * zmul(:,:)
        ENDDO
        DO jk = 1, nlay_s
           e_s(:,:,jk,kl) = e_s(:,:,jk,kl) + (e_s(:,:,jk,kl)-zicewrk(:,:,6+jk+nlay_s)) * zmul(:,:)
        ENDDO
        DO jk = 1, nlay_i
           t_i(:,:,jk,kl) = t_i(:,:,jk,kl) + (t_i(:,:,jk,kl)-zicewrk(:,:,6+jk+2*nlay_s)) * zmul(:,:)
        ENDDO
        DO jk = 1, nlay_i
           e_i(:,:,jk,kl) = e_i(:,:,jk,kl) + (e_i(:,:,jk,kl)-zicewrk(:,:,6+jk+2*nlay_s+nlay_i)) * zmul(:,:)
        ENDDO
        DO jk = 1, nlay_i
          sz_i(:,:,jk,kl) =sz_i(:,:,jk,kl) + (sz_i(:,:,jk,kl)-zicewrk(:,:,6+jk+2*nlay_s+2*nlay_i)) * zmul(:,:)
        ENDDO
        sst_m(:,:)     = sst_m(:,:)   + (sst_m(:,:)   - zicewrk(:,:,6+2*nlay_s+3*nlay_i+1) ) * zmul(:,:)
        sss_m(:,:)     = sss_m(:,:)   + (sss_m(:,:)   - zicewrk(:,:,6+2*nlay_s+3*nlay_i+2) ) * zmul(:,:)
        frq_m(:,:)     = frq_m(:,:)   + (frq_m(:,:)   - zicewrk(:,:,6+2*nlay_s+3*nlay_i+3) ) * zmul(:,:)
        o_i  (:,:,kl)  = o_i  (:,:,kl)+ (o_i  (:,:,kl)- zicewrk(:,:,6+2*nlay_s+3*nlay_i+4) ) * zmul(:,:)
        oa_i (:,:,kl)  = oa_i (:,:,kl)+ (oa_i (:,:,kl)- zicewrk(:,:,6+2*nlay_s+3*nlay_i+5) ) * zmul(:,:)
        !
        v_i(:,:,kl) = h_i(:,:,kl) * a_i(:,:,kl)
        v_s(:,:,kl) = h_s(:,:,kl) * a_i(:,:,kl)
        sv_i(:,:,kl) = s_i(:,:,kl) * v_i(:,:,kl)
        oa_i(:,:,kl) = o_i (:,:,kl) * a_i(:,:,kl)
        !
      ELSE
        CALL ctl_stop ('STOP', ' Step in sppt_ice is not valid')
      ENDIF
   END SUBROUTINE
#endif
   !
   SUBROUTINE sppt_glocon_3d ( ztu, cd_type, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_glocon ***
      !!
      !! ** Purpose :   Apply global conservation constraint
      !!                Note: not coded for vvl (code will stop at initializ.)
      !!----------------------------------------------------------------------
      !
      REAL(wp), INTENT(INOUT)      :: ztu(jpi,jpj,jpk)
      CHARACTER(len=1), INTENT(IN) :: cd_type
      INTEGER, INTENT(IN), OPTIONAL :: Kmm
      !
      INTEGER :: ji,jj,jk,jm
      REAL(wp) :: zsc, zww
      !
      jm = 1      ;   IF(PRESENT(Kmm)) jm = Kmm
      !
      zsc = 0._wp ;             zww = 0._wp
      !
      IF( cd_type .eq. 'T') THEN
          DO_3D( 0, 0, 0, 0, 1, jpkm1 )
             zsc = zsc + ztu(ji,jj,jk)*e3t(ji,jj,jk,jm)*e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk)*tmask_i(ji,jj)
             zww = zww +               e3t(ji,jj,jk,jm)*e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk)*tmask_i(ji,jj)
          END_3D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:,:) = ztu(:,:,:) - zsc*tmask(:,:,:)
      ELSEIF( cd_type .eq. 'U') THEN
          DO_3D( 0, 0, 0, 0, 1, jpkm1 )
             zsc = zsc + ztu(ji,jj,jk)*e3u(ji,jj,jk,jm)*e1u(ji,jj)*e2u(ji,jj)*umask(ji,jj,jk)*tmask_i(ji,jj)
             zww = zww +               e3u(ji,jj,jk,jm)*e1u(ji,jj)*e2u(ji,jj)*umask(ji,jj,jk)*tmask_i(ji,jj)
          END_3D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:,:) = ztu(:,:,:) - zsc*umask(:,:,:)
      ELSEIF( cd_type .eq. 'V') THEN
          DO_3D( 0, 0, 0, 0, 1, jpkm1 )
             zsc = zsc + ztu(ji,jj,jk)*e3v(ji,jj,jk,jm)*e1v(ji,jj)*e2v(ji,jj)*vmask(ji,jj,jk)*tmask_i(ji,jj)
             zww = zww +               e3v(ji,jj,jk,jm)*e1v(ji,jj)*e2v(ji,jj)*vmask(ji,jj,jk)*tmask_i(ji,jj)
          END_3D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:,:) = ztu(:,:,:) - zsc*vmask(:,:,:)
      ENDIF
      !
   END SUBROUTINE
   !
   SUBROUTINE sppt_glocon_2d ( ztu, cd_type, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_glocon ***
      !!
      !! ** Purpose :   Apply global conservation constraint
      !!                Note: not coded for vvl (code will stop at initializ.)
      !!----------------------------------------------------------------------
      !
      REAL(wp), INTENT(INOUT)      :: ztu(jpi,jpj)
      CHARACTER(len=1), INTENT(IN) :: cd_type
      INTEGER, INTENT(IN), OPTIONAL :: Kmm
      !
      INTEGER :: ji,jj
      REAL(wp) :: zsc, zww
      !
      zsc = 0._wp ;             zww = 0._wp
      !
      IF( cd_type .eq. 'T') THEN
          DO_2D( 0, 0, 0, 0 )
             zsc = zsc + ztu(ji,jj)*e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)*tmask_i(ji,jj)
             zww = zww +            e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)*tmask_i(ji,jj)
          END_2D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:) = ztu(:,:) - zsc*tmask(:,:,1)
      ELSEIF( cd_type .eq. 'U') THEN
          DO_2D( 0, 0, 0, 0 )
             zsc = zsc + ztu(ji,jj)*e1u(ji,jj)*e2u(ji,jj)*umask(ji,jj,1)*tmask_i(ji,jj)
             zww = zww +            e1u(ji,jj)*e2u(ji,jj)*umask(ji,jj,1)*tmask_i(ji,jj)
          END_2D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:) = ztu(:,:) - zsc*umask(:,:,1)
      ELSEIF( cd_type .eq. 'V') THEN
          DO_2D( 0, 0, 0, 0 )
             zsc = zsc + ztu(ji,jj)*e1v(ji,jj)*e2v(ji,jj)*vmask(ji,jj,1)*tmask_i(ji,jj)
             zww = zww +            e1v(ji,jj)*e2v(ji,jj)*vmask(ji,jj,1)*tmask_i(ji,jj)
          END_2D
          CALL mpp_sum('sppt_glocon',zsc)
          CALL mpp_sum('sppt_glocon',zww)
          IF( zww .gt. 0._wp ) zsc = zsc / zww
          ztu(:,:) = ztu(:,:) - zsc*vmask(:,:,1)
      ENDIF
      !
   END SUBROUTINE
   !
   !!======================================================================
END MODULE stosppt
