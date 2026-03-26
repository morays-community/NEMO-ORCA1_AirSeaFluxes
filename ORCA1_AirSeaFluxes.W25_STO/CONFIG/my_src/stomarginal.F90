MODULE stomarginal
   !!======================================================================
   !!                       ***  MODULE  stomarginaldiff  ***
   !! Purpose : transform marginal N(0,1) variates
   !!           to requested marginal distribution
   !!=====================================================================
   !!   sto_marginal      : transform to requested marginal distribution
   !!   sto_marginal_init : intitialize transformations
   !!----------------------------------------------------------------------
   USE stoarray        ! module with stochastic arrays to update
   ! user supplied external resources
   USE stoexternal, only : wp, lc, jpi, jpj, jpk

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: pi = 3.14159265358979323846_wp
   REAL(wp), PARAMETER :: halfpi = pi * 0.5_wp
   REAL(wp), PARAMETER :: invpi = 1.0_wp / pi

   PUBLIC sto_marginal, sto_marginal_init

CONTAINS

   SUBROUTINE sto_marginal(jsto)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_marginal  ***
      !!
      !! ** Purpose :   transform marginal N(0,1) variates
      !!                to requested marginal distribution
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: jsto   ! index of stochastic field in stoarray

      INTEGER :: jstoidx, kdim, jidx, nwrap
      REAL(wp) :: ave, std, var, zmin, zmax
      CHARACTER(len=lc) :: type_variate

      ! Get dimension of field
      kdim    = stofields(jsto)%dim
      ! Get index of field in 2d or 3d arrays
      jstoidx = stofields(jsto)%index

      ! Get parameters of the distributins
      type_variate = stofields(jsto)%type_variate
      ave = stofields(jsto)%ave
      std = stofields(jsto)%std
      zmin = stofields(jsto)%min
      zmax = stofields(jsto)%max

      ! Set index of field to update
      IF (kdim==2) THEN
         jidx = jpidx2d + jpidxsup2d
      ELSEIF (kdim==3) THEN
         jidx = jpidx3d + jpidxsup3d
      ELSEIF (kdim==0) THEN
         jidx = jpidx0d + jpidxsup0d
      ENDIF

      ! apply transformation according to requested distribution
      IF ( type_variate == 'normal' ) THEN
         ! for normal marginal distributions
         IF (kdim==2) THEN
            IF (zmin /= zmax ) THEN
               WHERE ( sto2d(:,:,jidx,jstoidx) < zmin )
                  sto2d(:,:,jidx,jstoidx) = zmin
               ELSEWHERE ( sto2d(:,:,jidx,jstoidx) > zmax )
                  sto2d(:,:,jidx,jstoidx) = zmax
               ENDWHERE
            ENDIF
            sto2d(:,:,jidx,jstoidx) = sto2d(:,:,jidx,jstoidx) * std + ave
         ELSEIF (kdim==3) THEN
            IF (zmin /= zmax ) THEN
               WHERE ( sto3d(:,:,:,jidx,jstoidx) < zmin )
                  sto3d(:,:,:,jidx,jstoidx) = zmin
               ELSEWHERE ( sto3d(:,:,:,jidx,jstoidx) > zmax )
                  sto3d(:,:,:,jidx,jstoidx) = zmax
               ENDWHERE
            ENDIF
            sto3d(:,:,:,jidx,jstoidx) = sto3d(:,:,:,jidx,jstoidx) * std + ave
         ELSEIF (kdim==0) THEN
            IF (zmin /= zmax ) THEN
               IF ( sto0d(jidx,jstoidx) < zmin ) THEN
                  sto0d(jidx,jstoidx) = zmin
               ELSEIF ( sto0d(jidx,jstoidx) > zmax ) THEN
                  sto0d(jidx,jstoidx) = zmax
               ENDIF
            ENDIF
            sto0d(jidx,jstoidx) = sto0d(jidx,jstoidx) * std + ave
         ENDIF
      ELSEIF ( type_variate == 'lognormal' ) THEN
         ! for lognormal distribution
         IF (ave <= 0._wp) STOP 'Bad parameter of lognormal distribution'
         ! compute parameters of reference normal distribution
         ! to have the requested mean and std for the lognormal distribution
         var = LOG( 1._wp + std*std/(ave*ave) )
         ave = LOG(ave) - 0.5_wp * var
         std = SQRT(var)
         ! apply transformation
         IF (kdim==2) THEN
            sto2d(:,:,jidx,jstoidx) = EXP( sto2d(:,:,jidx,jstoidx) * std + ave )
         ELSEIF (kdim==3) THEN
            sto3d(:,:,:,jidx,jstoidx) = EXP( sto3d(:,:,:,jidx,jstoidx) * std + ave )
         ELSEIF (kdim==0) THEN
            sto0d(jidx,jstoidx) = EXP( sto0d(jidx,jstoidx) * std + ave )
         ENDIF
      ELSEIF ( type_variate == 'wrapped_normal' ) THEN
         ! for wrapped normal marginal distributions (e.g. cyclic for angles, with zmin=0 and zmax=2pi)
         IF (zmin == zmax ) STOP 'Bad parameter of wrapped normal distribution'
         ! wrap normal distribution to fit in required interval
         IF (kdim==2) THEN
            sto2d(:,:,jidx,jstoidx) = sto2d(:,:,jidx,jstoidx) * std + ave
            WHERE ( (sto2d(:,:,jidx,jstoidx) < zmin) .OR. (sto2d(:,:,jidx,jstoidx) > zmax) )
               sto2d(:,:,jidx,jstoidx) = sto2d(:,:,jidx,jstoidx) - ( zmax - zmin ) * &
               &                         FLOOR( ( sto2d(:,:,jidx,jstoidx) - zmin ) / ( zmax - zmin ) )
            ENDWHERE
         ELSEIF (kdim==3) THEN
            sto3d(:,:,:,jidx,jstoidx) = sto3d(:,:,:,jidx,jstoidx) * std + ave
            WHERE ( (sto3d(:,:,:,jidx,jstoidx) < zmin) .OR. (sto3d(:,:,:,jidx,jstoidx) > zmax) )
               sto3d(:,:,:,jidx,jstoidx) = sto3d(:,:,:,jidx,jstoidx) - ( zmax - zmin ) * &
               &                           FLOOR( ( sto3d(:,:,:,jidx,jstoidx) - zmin ) / ( zmax - zmin ) )
            ENDWHERE
         ELSEIF (kdim==0) THEN
            sto0d(jidx,jstoidx) = sto0d(jidx,jstoidx) * std + ave
            IF ( (sto0d(jidx,jstoidx) < zmin) .OR. (sto0d(jidx,jstoidx) > zmax) ) THEN
               nwrap = FLOOR( ( sto0d(jidx,jstoidx) - zmin ) / ( zmax - zmin ) )
               sto0d(jidx,jstoidx) = sto0d(jidx,jstoidx) - nwrap * ( zmax - zmin )
            ENDIF
         ENDIF
      ELSEIF ( type_variate == 'bounded_atan' ) THEN
         ! for bounded distributions with atan transformation
         IF (zmin == zmax ) STOP 'Bad parameter of wrapped normal distribution'
         ! compute parameters of reference normal distribution
         ! to have the requested features for the resulting distribution
         ! Warning: input ave is not the mean, only the transformed value of the normal 0
         !          input std is not the std, only a measure of the spread for the user variable
         ave = TAN( pi * ( ave - zmin ) / ( zmax - zmin ) - halfpi )
         std = std * pi / ( zmax - zmin ) * ( 1._wp + ave * ave )
         ! apply transformation
         IF (kdim==2) THEN
            sto2d(:,:,jidx,jstoidx) = ATAN( sto2d(:,:,jidx,jstoidx) * std + ave )
            sto2d(:,:,jidx,jstoidx) = ( sto2d(:,:,jidx,jstoidx) + halfpi ) * invpi
            sto2d(:,:,jidx,jstoidx) = zmin + sto2d(:,:,jidx,jstoidx) * ( zmax - zmin )
         ELSEIF (kdim==3) THEN
            sto3d(:,:,:,jidx,jstoidx) = ATAN( sto3d(:,:,:,jidx,jstoidx) * std + ave )
            sto3d(:,:,:,jidx,jstoidx) = ( sto3d(:,:,:,jidx,jstoidx) + halfpi ) * invpi
            sto3d(:,:,:,jidx,jstoidx) = zmin + sto3d(:,:,:,jidx,jstoidx) * ( zmax - zmin )
         ELSEIF (kdim==0) THEN
            sto0d(jidx,jstoidx) = ATAN( sto0d(jidx,jstoidx) * std + ave )
            sto0d(jidx,jstoidx) = ( sto0d(jidx,jstoidx) + halfpi ) * invpi
            sto0d(jidx,jstoidx) = zmin + sto0d(jidx,jstoidx) * ( zmax - zmin )
         ENDIF
      ELSE
         STOP 'Bad name of marginal distribution in stomarginal'
      ENDIF

   END SUBROUTINE sto_marginal


   SUBROUTINE sto_marginal_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_marginal_init  ***
      !!
      !! ** Purpose :   initialize transformation to requested marginal distribution
      !!----------------------------------------------------------------------
      INTEGER :: jsto

      ! do nothing
   END SUBROUTINE sto_marginal_init

END MODULE stomarginal
