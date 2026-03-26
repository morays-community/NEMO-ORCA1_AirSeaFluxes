MODULE stowhite
   !!======================================================================
   !!                       ***  MODULE  stowhite  ***
   !! Stochastic parameters : get uncorrelated normal random numbers
   !!=====================================================================
   !! History :  4.0  ! 2024-01 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   sto_white      : get array of uncorrelated normal random numbers
   !!   sto_white_init : initialize random number generator
   !!----------------------------------------------------------------------
   USE stoexternal, only : wp, i4, i8
   USE storng_kiss      ! KISS random number generator
   USE storng_ziggurat  ! Ziggurat algorithm to generate normal numbers

   IMPLICIT NONE
   PRIVATE

   ! Public variable defining method to use
   ! CHARACTER(len=8), SAVE, PUBLIC :: c_rngtype='kiss64'
   CHARACTER(len=8), SAVE, PUBLIC :: c_rngtype='kiss32'
   ! CHARACTER(len=8), SAVE, PUBLIC :: c_rngtype='shr3'
   CHARACTER(len=8), SAVE, PUBLIC :: c_normal_algo='ziggurat'
   !CHARACTER(len=8), SAVE, PUBLIC :: c_normal_algo='polar'
   INTEGER :: normal_algo

   ! Public routines
   PUBLIC :: sto_white, sto_white_init

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stopar.F90 13255 2020-07-06 15:41:29Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sto_white( psto0d, psto1d, psto2d, psto3d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_white  ***
      !!
      !! ** Purpose :   fill input array with uncorrelated Gaussian numbers
      !!                (mean = 0, standard deviation = 1)
      !!----------------------------------------------------------------------
      REAL(wp)                  , INTENT(out), OPTIONAL :: psto0d
      REAL(wp), DIMENSION(:    ), INTENT(out), OPTIONAL :: psto1d
      REAL(wp), DIMENSION(:,:  ), INTENT(out), OPTIONAL :: psto2d
      REAL(wp), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: psto3d

      LOGICAL :: l0d, l1d, l2d, l3d   ! which optional arguments are present
      INTEGER :: ji, jj, jk           ! array indices
      INTEGER :: jp1di                ! size of 1d array
      INTEGER :: jp2di, jp2dj         ! size of 2d array
      INTEGER :: jp3di, jp3dj, jp3dk  ! size of 3d array

      ! Check presence of optional arguments
      l0d = PRESENT(psto0d)
      l1d = PRESENT(psto1d)
      l2d = PRESENT(psto2d)
      l3d = PRESENT(psto3d)

      ! Get size of arrays in argument
      IF (l1d) THEN
        jp1di = SIZE(psto1d,1)
      ENDIF
      IF (l2d) THEN
        jp2di = SIZE(psto2d,1)
        jp2dj = SIZE(psto2d,2)
      ENDIF
      IF (l3d) THEN
        jp3di = SIZE(psto3d,1)
        jp3dj = SIZE(psto3d,2)
        jp3dk = SIZE(psto3d,3)
      ENDIF

      ! Select type of random number generator to use
      IF (normal_algo==0) THEN
         IF (l0d) THEN
            psto0d = zig_normal( )
         ENDIF
         IF (l1d) THEN
            DO ji=1,jp1di
               psto1d(ji) = zig_normal( )
            ENDDO
         ENDIF
         IF (l2d) THEN
            DO ji=1,jp2di
            DO jj=1,jp2dj
               psto2d(ji,jj) = zig_normal( )
            ENDDO
            ENDDO
         ENDIF
         IF (l3d) THEN
            DO ji=1,jp3di
            DO jj=1,jp3dj
            DO jk=1,jp3dk
               psto3d(ji,jj,jk) = zig_normal( )
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ELSEIF (normal_algo==1) THEN
         IF (l0d) THEN
            psto0d = kiss_normal()
         ENDIF
         IF (l1d) THEN
            DO ji=1,jp1di
               psto1d(ji) = kiss_normal()
            ENDDO
         ENDIF
         IF (l2d) THEN
            DO ji=1,jp2di
            DO jj=1,jp2dj
               psto2d(ji,jj) = kiss_normal()
            ENDDO
            ENDDO
         ENDIF
         IF (l3d) THEN
            DO ji=1,jp3di
            DO jj=1,jp3dj
            DO jk=1,jp3dk
               psto3d(ji,jj,jk) = kiss_normal()
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ELSEIF (normal_algo==2) THEN
         IF (l0d) THEN
            psto0d = shr3_normal( )
         ENDIF
         IF (l1d) THEN
            DO ji=1,jp1di
               psto1d(ji) = shr3_normal( )
            ENDDO
         ENDIF
         IF (l2d) THEN
            DO ji=1,jp2di
            DO jj=1,jp2dj
               psto2d(ji,jj) = shr3_normal( )
            ENDDO
            ENDDO
         ENDIF
         IF (l3d) THEN
            DO ji=1,jp3di
            DO jj=1,jp3dj
            DO jk=1,jp3dk
               psto3d(ji,jj,jk) = shr3_normal( )
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ELSE
         STOP 'Bad type of algorithm in stowhite'
      ENDIF

   END SUBROUTINE sto_white


   SUBROUTINE sto_white_init( kseedindex )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_white_init  ***
      !!
      !! ** Purpose :   initialize and seed random number generator
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kseedindex
      INTEGER(KIND=i8) :: zseed1_8, zseed2_8, zseed3_8, zseed4_8
      INTEGER(KIND=i4) :: zseed1_4, zseed2_4, zseed3_4, zseed4_4
      INTEGER :: jseed

      ! Set type of algorithm to generate normal random numbers
      SELECT CASE (c_normal_algo)
      CASE('ziggurat')
         normal_algo = 0
      CASE('polar')
         SELECT CASE (c_rngtype)
         CASE('kiss64')
            normal_algo = 1
         CASE('kiss32')
            normal_algo = 1
         CASE('shr3')
            normal_algo = 2
         CASE DEFAULT
            STOP 'Bad type of random number generator in stowhite_init'
         END SELECT
      CASE DEFAULT
         STOP 'Bad type of algorithm in stowhite_init'
      END SELECT

      ! Seed random number generator
      SELECT CASE (c_rngtype)
      CASE('kiss64')
         ! Compute seed for the 64-bit KISS random number generator
         kiss_32bits=.FALSE.
         CALL kiss_reset( )
         DO jseed = 0, kseedindex
            zseed1_8 = kiss64() ; zseed2_8 = kiss64()
            zseed3_8 = kiss64() ; zseed4_8 = kiss64()
         END DO
         CALL kiss_seed( zseed1_8, zseed2_8, zseed3_8, zseed4_8 )
      CASE('kiss32')
         ! Compute seed for the 32-bit KISS random number generator
         kiss_32bits=.TRUE.
         CALL kiss_reset( )
         DO jseed = 0, kseedindex
            zseed1_4 = kiss32() ; zseed2_4 = kiss32()
            zseed3_4 = kiss32() ; zseed4_4 = kiss32()
         END DO
         CALL kiss_seed( zseed1_4, zseed2_4, zseed3_4, zseed4_4 )
      CASE('shr3')
         ! Compute seed for the 32-bit shr3 random number generator
         CALL shr3_reset( )
         DO jseed = 0, kseedindex
            zseed1_4 = kiss32()
         END DO
         CALL shr3_seed(zseed1_4)
      CASE DEFAULT
         STOP 'Bad type of random number generator in stowhite_init'
      END SELECT

      ! Initialize Ziggurat algorithm
      IF (normal_algo==0) THEN
         zig_rngtype=c_rngtype
         call zig_set()
      ENDIF

   END SUBROUTINE sto_white_init

END MODULE stowhite

