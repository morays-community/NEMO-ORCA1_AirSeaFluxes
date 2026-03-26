MODULE stoskeb
   !!======================================================================
   !!                       ***  MODULE stoskeb  ***
   !!
   !! Purpose : Implement SKBE stochastic parameterization
   !!======================================================================
   USE stoarray
   USE par_oce
   USE dom_oce
   USE zdf_oce
   USE oce
   USE phycst
   USE in_out_manager
   USE iom
   USE fldread        ! read input fields
   USE lib_mpp
   USE ldfdyn

   IMPLICIT NONE
   PRIVATE

   ! Publicize variables from stoarray
   PUBLIC :: ln_sto_skeb

   ! Switches to activate SKEB scheme
   LOGICAL, SAVE :: ln_active_skeb = .FALSE.

   ! SKEB Options (with default values)
   REAL(wp), SAVE :: rn_skeb = 0.4_wp, rn_kh =1._wp, rn_kc=1._wp
   INTEGER,  SAVE :: nn_skeb_freq
   REAL(wp), SAVE :: rn_skeb_tau = 5._wp
   REAL(wp), SAVE :: rn_skeb_stdev = 1.0_wp
   INTEGER, SAVE  :: skeb_filter_pass = 50
   LOGICAL, SAVE  :: ln_skeb_num, ln_skeb_con, ln_skeb_eke
   REAL(wp), SAVE :: rn_beta_num = 1._wp
   REAL(wp), SAVE :: rn_beta_con = 1._wp
   REAL(wp), SAVE :: rn_beta_eke = 1._wp

   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: psi, ustar, vstar, dnum, dcon, deke
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:  ) :: rn_kc2, rn_kh2
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_uvc

   ! SKEB stochastic array indices
   INTEGER, PUBLIC :: jsto_skeb

   ! Public routines
   PUBLIC sto_skeb, sto_skeb_init, skeb_apply

   !! * Substitutions
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE sto_skeb ( kstp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_skeb  ***
      !!
      !! ** Purpose :   Computation of SKEB perturbations
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)    :: kstp
      !
      ! Output stochastic field
      CALL iom_put( 'skeb_sto',     stofields(jsto_skeb)%sto2d(:,:) )
      !
      CALL skeb_comp( kstp )
      !
   END SUBROUTINE

   SUBROUTINE skeb_apply ( kstp, Nbb, Nnn, Nrhs, pu, pv, Naa  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_apply  ***
      !!
      !! ** Purpose :   Application of SKEB perturbation
      !!                Convective and Numerical energy dissipation are
      !!                multiplied by the beta terms
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)    :: kstp, Nbb, Nnn, Nrhs, Naa
      REAL(wp),INTENT(INOUT) :: pu(jpi,jpj,jpk,jpt), pv(jpi,jpj,jpk,jpt)
      INTEGER :: ji,jj,jk
      !
      psi = 0._wp
      DO jk=1,jpkm1
          !IF( ln_skeb_num ) psi(:,:,jk) = rn_skeb * sqrt(rn_beta_num*dnum(:,:,jk) + rn_beta_con*dcon(:,:,jk)+ &
          !& rn_beta_eke * abs( deke(:,:,jk) )) * stofields(jsto_skeb)%sto2d(:,:)
          IF( ln_skeb_num ) psi(:,:,jk) = rn_beta_num*dnum(:,:,jk) 
          IF( ln_skeb_con ) psi(:,:,jk) = psi(:,:,jk) + rn_beta_con*dcon(:,:,jk) 
          IF( ln_skeb_eke ) psi(:,:,jk) = psi(:,:,jk) + rn_beta_eke*abs(deke(:,:,jk))
          psi(:,:,jk) = rn_skeb * sqrt( psi(:,:,jk) ) * stofields(jsto_skeb)%sto2d(:,:) * tmask(:,:,jk) 
      ENDDO
      !
      ustar = 0._wp
      vstar = 0._wp
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
           ustar(ji,jj,jk) =   ( psi(ji,jj,jk)-psi(ji,jj-1,jk) )*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
           vstar(ji,jj,jk) = - ( psi(ji,jj,jk)-psi(ji-1,jj,jk) )*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj)
      END_3D
      ! 
      CALL lbc_lnk("skeb_apply", ustar,'U',-1._wp)
      CALL lbc_lnk("skeb_apply", vstar,'V',-1._wp)
      ! 
      WHERE( .NOT. abs(ustar) .lt. 1.e+6 ) ustar = 0._wp
      WHERE( .NOT. abs(vstar) .lt. 1.e+6 ) vstar = 0._wp
      !
      ustar(:,:,:) = ustar(:,:,:) * umask(:,:,:)
      vstar(:,:,:) = vstar(:,:,:) * vmask(:,:,:)
      !
      pu(:,:,:,Naa) = pu(:,:,:,Naa) + ustar(:,:,:)*umask(:,:,:)*rDt
      pv(:,:,:,Naa) = pv(:,:,:,Naa) + vstar(:,:,:)*vmask(:,:,:)*rDt
      !
!      ub(:,:,:) = ub(:,:,:) + ustar(:,:,:)*umask(:,:,:)
!      vb(:,:,:) = vb(:,:,:) + vstar(:,:,:)*vmask(:,:,:)
      !
   END SUBROUTINE

   SUBROUTINE sto_skeb_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_skeb_init  ***
      !!
      !! This routine is calle at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------
      INTEGER  :: ios
      CHARACTER(len=lc)  ::   cn_dir
      TYPE(FLD_N)        ::  sn_uc, sn_vc
      !
      NAMELIST/namsto_skeb/&
      rn_skeb, nn_skeb_freq, rn_kh, rn_kc, skeb_filter_pass, &
      rn_skeb_stdev, rn_skeb_tau,&
      rn_beta_num, rn_beta_con, rn_beta_eke, sn_uc, sn_vc, cn_dir, &
      ln_skeb_num, ln_skeb_con, ln_skeb_eke

      ! Read namelist block corresponding to this stochastic scheme
      READ_NML_REF(numnam_sto,namsto_skeb)
      READ_NML_CFG(numnam_sto,namsto_skeb)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stoskeb : Stochastic SKEB scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto_skeb: '
         WRITE(numout,*)
      ENDIF

      ln_active_skeb = .FALSE.
      IF( ln_skeb_num .OR. ln_skeb_con .OR. ln_skeb_eke ) THEN
          ln_active_skeb = .TRUE. 
      ENDIF
      IF( .NOT. ln_active_skeb ) &
      & CALL ctl_stop ( 'STOP', 'ln_sto_skeb TRUE without energy terms' )

      ! Define characteristics of each stochastic array
      IF( ln_active_skeb ) THEN
         ! Request index for a each stochastic array
         CALL sto_array_request_new(jsto_skeb)

         ! Set features of the corresponding stochastic field
         stofields(jsto_skeb)%dim=2

         stofields(jsto_skeb)%type_t='arn'
         stofields(jsto_skeb)%corr_t=rn_skeb_tau * 86400._wp /rdt
         stofields(jsto_skeb)%nar_order=1
         stofields(jsto_skeb)%nar_update=1

         stofields(jsto_skeb)%type_xy='diffusive'
         stofields(jsto_skeb)%diff_passes=skeb_filter_pass
         stofields(jsto_skeb)%diff_type=1

         stofields(jsto_skeb)%std=rn_skeb_stdev

         ALLOCATE( psi(jpi,jpj,jpk), ustar(jpi,jpj,jpk), vstar(jpi,jpj,jpk) )

      ENDIF

      IF( ln_skeb_num ) THEN
              ALLOCATE( dnum(jpi,jpj,jpk) )
              dnum(:,:,:) = 0._wp
              ALLOCATE( rn_kh2(jpi,jpj) )
              IF( rn_kh .LE. 0._wp ) THEN
                  CALL iom_open('skeb', ios )
                  CALL iom_get(ios,jpdom_auto,'rn_kh',rn_kh2)
                  CALL iom_close( ios )
              ELSE
                  rn_kh2(:,:) = rn_kh
              ENDIF
      ENDIF

      IF( ln_skeb_con ) THEN
              ALLOCATE( dcon(jpi,jpj,jpk) )
              dcon(:,:,:) = 0._wp
              ALLOCATE( rn_kc2(jpi,jpj) )
              IF( rn_kc .LE. 0._wp ) THEN
                  CALL iom_open('skeb', ios )
                  CALL iom_get(ios,jpdom_auto,'rn_kc',rn_kc2)
                  CALL iom_close( ios )
              ELSE
                  rn_kc2(:,:) = rn_kc
              ENDIF
      ENDIF

      IF( ln_skeb_eke ) THEN
              ALLOCATE( deke(jpi,jpj,jpk) )
              deke(:,:,:) = 0._wp
              ALLOCATE( sf_uvc(2) ) 
              ALLOCATE( sf_uvc(1)%fnow(jpi,jpj,1) ) 
              ALLOCATE( sf_uvc(2)%fnow(jpi,jpj,1) ) 
              IF( sn_uc%ln_tint ) ALLOCATE( sf_uvc(1)%fdta(jpi,jpj,1,2) )
              IF( sn_vc%ln_tint ) ALLOCATE( sf_uvc(2)%fdta(jpi,jpj,1,2) )
              CALL fld_fill( sf_uvc, (/sn_uc, sn_vc/), cn_dir, 'sto_skeb_init',&
              & 'Climatological current vector', 'namsto_skeb' )
              sf_uvc(1)%zsgn = -1._wp   ;   sf_uvc(2)%zsgn = -1._wp
      ENDIF

   END SUBROUTINE sto_skeb_init

   SUBROUTINE skeb_comp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_comp  ***
      !!
      !! ** Purpose :   Computation of energy dissipation terms
      !!                This is a wrapper to the enrgy terms computation
      !!                routines
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: kt
      !
      IF(ln_skeb_num) CALL skeb_dnum ( kt )
      IF(ln_skeb_con) CALL skeb_dcon ( kt )
      IF(ln_skeb_eke) CALL skeb_deke ( kt )
      !
   END SUBROUTINE skeb_comp

   SUBROUTINE skeb_dnum ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_dnum  ***
      !!
      !! ** Purpose :   Computation of numerical energy dissipation
      !!                For later use in the SKEB scheme
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: kt
      REAL(wp) :: ds,dt,dtot
      INTEGER :: ji,jj,jk
      !
      IF( kt .eq. nit000 .OR. MOD( kt - 1, nn_skeb_freq ) == 0 ) THEN
        !
        DO_3D( 0, 0, 0, 0, 1, jpkm1 )
                 ! Shear
                 ds = (vv(ji,jj,jk,Nbb)-vv(ji-1,jj,jk,Nbb))*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj) + &
                      (uu(ji,jj,jk,Nbb)-uu(ji,jj-1,jk,Nbb))*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
                 ! Tension
                 dt = (vv(ji,jj,jk,Nbb)-vv(ji-1,jj,jk,Nbb))*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e2v(ji,jj) + &
                      (uu(ji,jj,jk,Nbb)-uu(ji,jj-1,jk,Nbb))*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e1u(ji,jj)
                 !
                 dtot = sqrt( ds*ds + dt*dt ) * tmask(ji,jj,jk)
                 dnum(ji,jj,jk) = dtot*dtot*dtot*rn_kh2(ji,jj)*rn_kh2(ji,jj)*e1t(ji,jj)*e2t(ji,jj)
                 !
        END_3D
        !
        CALL lbc_lnk("skeb_dnum",dnum,'T',1._wp)
        !
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE skeb_deke ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_eked  ***
      !!
      !! ** Purpose :   Computation of eddy dissipation
      !!                For later use in the SKEB scheme
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: kt
      REAL(wp) :: dudx,dvdy,Ah
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ut, vt
      INTEGER :: ji,jj,jk
      ! 
      IF( kt .eq. nit000 .OR. MOD( kt - 1, nn_skeb_freq ) == 0 ) THEN
         !
         ALLOCATE ( ut(jpi,jpj,jpk), vt(jpi,jpj,jpk) )
         ut(:,:,: ) = 0._wp
         vt(:,:,: ) = 0._wp
         !
         CALL fld_read( kt, nn_skeb_freq, sf_uvc )
         !
         ut = ( uu(:,:,:,Nbb) - sf_uvc(1)%fnow(:,:,:) ) * umask
         vt = ( vv(:,:,:,Nbb) - sf_uvc(2)%fnow(:,:,:) ) * vmask
         !
         ! Horizontal eddy dissipation
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
                  ! 
                  dudx = (ut(ji,jj,jk)-ut(ji-1,jj,jk))/e1t(ji,jj) * &
                        & tmask(ji,jj,jk) * umask(ji-1,jj,jk)
                  dvdy = (vt(ji,jj,jk)-vt(ji,jj-1,jk))/e2t(ji,jj) * &
                        & tmask(ji,jj,jk) * vmask(ji,jj-1,jk)
                  Ah  = ahmt(ji,jj,jk)
                  ! 
                  deke(ji,jj,jk) = rho0*(Ah)*(dudx*dudx+dvdy*dvdy)
                  ! 
         END_3D
         !
         DEALLOCATE (ut, vt)
         CALL lbc_lnk("skeb_deke", deke,'T',1._wp)
         !
      ENDIF
   
   END SUBROUTINE

   SUBROUTINE skeb_dcon ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_dcon  ***
      !!
      !! ** Purpose :   Computation of convective energy dissipation
      !!                For later use in the SKEB scheme
      !!                The formulation relies on convective mass flux.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: kt
      REAL(wp) :: mf1,mf2
      INTEGER  :: ji,jj,jk
      ! 
      IF( kt .eq. nit000 .OR. MOD( kt - 1, nn_skeb_freq ) == 0 ) THEN
         ! 
         IF( ln_zad_Aimp ) THEN
             DO_3D( 0, 0, 0, 0, 2, jpkm1 )
                 ! 
                 mf1 = (ww(ji,jj,jk+1)+wi(ji,jj,jk+1))*e1t(ji,jj)*e2t(ji,jj)
                 mf2 = (ww(ji,jj,jk-1)+wi(ji,jj,jk-1))*e1t(ji,jj)*e2t(ji,jj)
                 dcon(ji,jj,jk) = rn_kc2(ji,jj)*(mf1-mf2)*(mf1-mf2)*tmask(ji,jj,jk+1)/(E3w_0(ji,jj,jk)*rho0*rho0)
                 ! 
             END_3D
         ELSE
             DO_3D( 0, 0, 0, 0, 2, jpkm1 )
                 ! 
                 mf1 = ww(ji,jj,jk+1)*e1t(ji,jj)*e2t(ji,jj)
                 mf2 = ww(ji,jj,jk-1)*e1t(ji,jj)*e2t(ji,jj)
                 dcon(ji,jj,jk) = rn_kc2(ji,jj)*(mf1-mf2)*(mf1-mf2)*tmask(ji,jj,jk+1)/(E3w_0(ji,jj,jk)*rho0*rho0)
                 ! 
             END_3D
         ENDIF
         ! 
         CALL lbc_lnk("skeb_dcon", dcon,'T',1._wp)
         ! 
      ENDIF
      ! 
   END SUBROUTINE

END MODULE stoskeb
