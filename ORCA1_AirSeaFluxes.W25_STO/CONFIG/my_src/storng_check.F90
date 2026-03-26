MODULE storng_check
   !!======================================================================
   !!                       ***  MODULE storng_check  ***
   !!
   !!  Purpose: check relative performance of random number generators
   !!
   !!======================================================================
   USE storng_kiss
   USE storng_ziggurat

   IMPLICIT NONE
   PRIVATE

   PUBLIC sto_rng_test

CONTAINS

   SUBROUTINE sto_rng_test
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE sto_rng_test  ***
      !!
      !! ** `Purpose : Test efficiency of the various options
      !!               for the random number generator
      !!
      !! We test 3 random number generators:
      !!    -kiss64 (period ~ 2^250, with four 64-bit seeds), 
      !!    -kiss32 (period ~ 2^123, with four 32-bit seeds),
      !!    -shr3   (period ~ 2^64,  with one  32-bit seed)
      !! together with 2 methods to generate normal numbers:
      !!    -the polar method (Marsaglia, 1964),
      !!    -the ziggurat method (Marsaglia, 2000).
      !!
      !! Typical results obtained (in seconds), with jpran=10^8
      !! polar    (Marsaglia, 1964) + kiss64 (Marsaglia, 2009)   1.89253000000000
      !! polar    (Marsaglia, 1964) + kiss32 (Marsaglia, 1992)   1.50936100000000
      !! polar    (Marsaglia, 1964) + shr3 (simple & fast rng)   1.14710500000000
      !! ziggurat (Marsaglia, 2000) + kiss64 (Marsaglia, 2009)   1.26209300000000     
      !! ziggurat (Marsaglia, 2000) + kiss32 (Marsaglia, 1992)  0.693772999999999
      !! ziggurat (Marsaglia, 2000) + shr3 (simple & fast rng)  0.523565000000000
      !! (the last 4 are obtained with real compuations in single precision)
      !!
      !! In view of these results, the default choice used in the code is:
      !! ziggurat (Marsaglia, 2000) + kiss32 (Marsaglia, 1992)
      !! which combines a long enough period with good efficiency.
      !! Options are included in the namelist to modify the default.
      !!----------------------------------------------------------------------
      USE stoexternal, only : sp, dp, wp, i4, i8

      INTEGER, PARAMETER :: jpran_seed=100  ! Loop iterations to compute the seed
      INTEGER, PARAMETER :: jpran=100000000 ! Number of call to rng in the test

      ! Kiss seeds
      INTEGER(KIND=i8) :: zseed1_8, zseed2_8, zseed3_8, zseed4_8
      INTEGER(KIND=i4) :: zseed1_4, zseed2_4, zseed3_4, zseed4_4

      ! Timing variables
      REAL(wp) :: start_time, end_time, elapsed_time

      ! Dummy variables 
      INTEGER :: jran
      REAL(dp) :: zran8
      REAL(sp) :: zran4

      print *, '==================================================='
      print *, 'Test of performance of the random number generators'
      print *, '(comment it out in sto_mod_init to skip the test)'
      print *, '==================================================='

      ! Compute seed for the 64-bit KISS random number generator
      kiss_32bits=.FALSE.
      CALL kiss_reset( )
      DO jran = 1, jpran_seed
         zseed1_8 = kiss64() ; zseed2_8 = kiss64() ; zseed3_8 = kiss64() ; zseed4_8 = kiss64()
      END DO

      ! Compute seed for the 32-bit KISS random number generator
      kiss_32bits=.TRUE.
      CALL kiss_reset( )
      DO jran = 1, jpran_seed
         zseed1_4 = kiss32() ; zseed2_4 = kiss32() ; zseed3_4 = kiss32() ; zseed4_4 = kiss32()
      END DO


      ! Test 64-bit KISS, with classic polar method to obtain normal numbers
      kiss_32bits = .FALSE.
      CALL kiss_seed( zseed1_8, zseed2_8, zseed3_8, zseed4_8 )

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran8 = kiss_normal()
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'polar mt (Marsaglia, 1964) + kiss64 (Marsaglia, 2009)', elapsed_time

      ! Test 32-bit KISS, with classic polar method to obtain normal numbers
      kiss_32bits = .TRUE.
      CALL kiss_seed( zseed1_4, zseed2_4, zseed3_4, zseed4_4 )

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran8 = kiss_normal()
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'polar mt (Marsaglia, 1964) + kiss32 (Marsaglia, 1992)', elapsed_time

      ! Test the shr3 simple and fast rng, with classic polar method to obtain normal numbers
      CALL shr3_seed(zseed1_4)

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran4 = shr3_normal( )
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'polar mt (Marsaglia, 1964) + shr3 (simple & fast rng)', elapsed_time

      ! Test 64-bit KISS, with Ziggurat method to obtain normal numbers
      kiss_32bits = .FALSE.
      CALL kiss_seed( zseed1_8, zseed2_8, zseed3_8, zseed4_8 )

      zig_rngtype='kiss64'
      call zig_set()

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran4 = zig_normal( )
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'ziggurat (Marsaglia, 2000) + kiss64 (Marsaglia, 2009)', elapsed_time

      ! Test 32-bit KISS, with Ziggurat method to obtain normal numbers
      kiss_32bits = .TRUE.
      CALL kiss_seed( zseed1_4, zseed2_4, zseed3_4, zseed4_4 )

      zig_rngtype='kiss32'
      call zig_set()

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran4 = zig_normal( )
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'ziggurat (Marsaglia, 2000) + kiss32 (Marsaglia, 1992)', elapsed_time


      ! Test the shr3 simple and fast rng, with Ziggurat method to obtain normal numbers
      call shr3_seed(zseed1_4)

      zig_rngtype='shr3'
      call zig_set()

      call CPU_TIME(start_time)
      DO jran = 1, jpran
        zran4 = zig_normal( )
      END DO
      call CPU_TIME(end_time)
      elapsed_time = end_time - start_time
      print *, 'ziggurat (Marsaglia, 2000) + shr3 (simple & fast rng)', elapsed_time

   END SUBROUTINE sto_rng_test
   !!======================================================================
END MODULE storng_check
