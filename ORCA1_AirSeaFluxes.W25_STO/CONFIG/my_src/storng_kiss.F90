MODULE storng_kiss
   !!======================================================================
   !!                       ***  MODULE  storng_kiss  ***
   !! Random number generator, used in NEMO stochastic parameterization
   !!
   !!=====================================================================
   !! History :  3.3  ! 2011-10 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! The module is based on (and includes) the
   !! 64-bit KISS (Keep It Simple Stupid) random number generator
   !! distributed by George Marsaglia :
   !! http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/a85bf5f2a97f5a55
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   kiss64        : 64-bit KISS random number generator (period ~ 2^250)
   !!   kiss32        : 32-bit KISS random number generator (period ~ 2^123)
   !!   kiss_seed     : Define seeds for KISS random number generator
   !!   kiss_state    : Get current state of KISS random number generator
   !!   kiss_save     : Save current state of KISS (for future restart)
   !!   kiss_load     : Load the saved state of KISS
   !!   kiss_reset    : Reset the default seeds
   !!   kiss_check    : Check the KISS pseudo-random sequence
   !!   kiss_uniform  : Real random numbers with uniform distribution in [0,1]
   !!   kiss_normal   : Real random numbers with normal distribution N(0,1)
   !!   kiss_gamma    : Real random numbers with gamma distribution Gamma(k,1)
   !!   kiss_sample   : Select a random sample from a set of integers
   !!
   !!   ---CURRENTLY NOT USED IN NEMO :
   !!   kiss_save, kiss_load, kiss_check, kiss_gamma, kiss_sample
   !!----------------------------------------------------------------------
   USE stoexternal, only : wp, i4, i8

   IMPLICIT NONE
   PRIVATE

   ! Public functions/subroutines
   PUBLIC :: kiss32, kiss64, kiss_seed, kiss_state, kiss_reset
   PUBLIC :: kiss_save, kiss_load, kiss_check
   PUBLIC :: kiss_uniform, kiss_normal, kiss_gamma, kiss_sample

   ! Public variable to choose between 32-bit and 64-bit generators
   LOGICAL, SAVE, PUBLIC :: kiss_32bits=.TRUE.

   ! Default/initial seeds (64-bit kiss)
   INTEGER(KIND=i8), SAVE :: x8=1234567890987654321_i8
   INTEGER(KIND=i8), SAVE :: y8=362436362436362436_i8
   INTEGER(KIND=i8), SAVE :: z8=1066149217761810_i8
   INTEGER(KIND=i8), SAVE :: w8=123456123456123456_i8

   ! Default/initial seeds (32-bit kiss)
   INTEGER(KIND=i4), SAVE :: x4=123456789_i4
   INTEGER(KIND=i4), SAVE :: y4=362436069_i4
   INTEGER(KIND=i4), SAVE :: z4=521288629_i4
   INTEGER(KIND=i4), SAVE :: w4=916191069_i4

   ! Parameters to generate real random variates
   REAL(KIND=wp), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0
   REAL(wp), PARAMETER    ::  scaling_uni1_8 = 0.5_wp/(real(huge(x8),wp)+1.0_wp)
   REAL(wp), PARAMETER    ::  scaling_uni2_8 = 1.0_wp/(real(huge(x8),wp)+1.0_wp)
   REAL(wp), PARAMETER    ::  scaling_uni1_4 = 0.5_wp/(real(huge(x4),wp)+1.0_wp)
   REAL(wp), PARAMETER    ::  scaling_uni2_4 = 1.0_wp/(real(huge(x4),wp)+1.0_wp)

   ! Variables to store 2 Gaussian random numbers with current index (ig)
   INTEGER(KIND=i8), SAVE :: ig=1
   REAL(KIND=wp), SAVE :: gran1, gran2

   ! Interface to 32-bits or 64-bits routines
   INTERFACE kiss_seed
      MODULE PROCEDURE kiss_seed_32, kiss_seed_64
   END INTERFACE

   INTERFACE kiss_state
      MODULE PROCEDURE kiss_state_32, kiss_state_64
   END INTERFACE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynhpg.F90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION kiss64()
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION kiss64  ***
      !!
      !! ** Purpose :   64-bit KISS random number generator
      !!
      !! ** Method  :   combine several random number generators:
      !!                (1) Xorshift (XSH), period 2^64-1,
      !!                (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
      !!                (3) Congruential generator (CNG), period 2^64.
      !!
      !!                overall period:
      !!                (2^250+2^192+2^64-2^186-2^129)/6
      !!                            ~= 2^(247.42) or 10^(74.48)
      !!
      !!                set your own seeds with 'kiss_seed'
      ! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: kiss64, t

      t = ISHFT(x8,58) + w8
      IF (s(x8).eq.s(t)) THEN
         w8 = ISHFT(x8,-6) + s(x8)
      ELSE
         w8 = ISHFT(x8,-6) + 1 - s(x8+t)
      ENDIF
      x8 = t + x8
      y8 = m( m( m(y8,13_i8), -17_i8 ), 43_i8 )
      z8 = 6906969069_i8 * z8 + 1234567_i8

      kiss64 = x8 + y8 + z8

      CONTAINS

         FUNCTION s(k)
            INTEGER(KIND=i8) :: s, k
            s = ISHFT(k,-63)
         END FUNCTION s

         FUNCTION m(k, n)
            INTEGER(KIND=i8) :: m, k, n
            m =  IEOR(k, ISHFT(k, n) )
         END FUNCTION m

   END FUNCTION kiss64


   FUNCTION kiss32()
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION kiss32  ***
      !!
      !! ** Purpose :   32-bit KISS random number generator
      !!                (from https://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90)
      !!
      !! ** Method  :   combine several random number generators:
      !!                (1) Ccongruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
      !!                (2) A 3-shift shift-register generator, period 2^32-1,
      !!                (3) Two 16-bit multiply-with-carry generators,
      !!                       period 597273182964842497>2^59
      !!                Overall period > 2^123
      !!
      !!                set your own seeds with 'kiss_seed'
      ! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(i4)          :: kiss32

      x4 = 69069 * x4 + 1327217885
      y4= ieor (y4, ishft (y4, 13))
      y4= ieor (y4, ishft (y4, -17))
      y4= ieor (y4, ishft (y4, 5))
      z4 = 18000 * iand (z4, 65535) + ishft (z4, - 16)
      w4 = 30903 * iand (w4, 65535) + ishft (w4, - 16)
      kiss32 = ishft(x4 + y4 + ishft (z4, 16) + w4 , -1)

   END FUNCTION kiss32


   SUBROUTINE kiss_seed_64(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_seed  ***
      !!
      !! ** Purpose :   Define seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      x8 = ix
      y8 = iy
      z8 = iz
      w8 = iw

   END SUBROUTINE kiss_seed_64


   SUBROUTINE kiss_seed_32(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_seed  ***
      !!
      !! ** Purpose :   Define seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i4) :: ix, iy, iz, iw

      x4 = ix
      y4 = iy
      z4 = iz
      w4 = iw

   END SUBROUTINE kiss_seed_32


   SUBROUTINE kiss_state_64(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_state  ***
      !!
      !! ** Purpose :   Get current state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      ix = x8 ; iy = y8 ; iz = z8 ; iw = w8

   END SUBROUTINE kiss_state_64


   SUBROUTINE kiss_state_32(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_state  ***
      !!
      !! ** Purpose :   Get current state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i4) :: ix, iy, iz, iw

      ix = x4 ; iy = y4 ; iz = z4 ; iw = w4

   END SUBROUTINE kiss_state_32


   SUBROUTINE kiss_reset()
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_reset  ***
      !!
      !! ** Purpose :   Reset the default seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE

      IF (kiss_32bits) THEN
        x4=123456789_i4
        y4=362436069_i4
        z4=521288629_i4
        w4=916191069_i4
      ELSE
        x8=1234567890987654321_i8
        y8=362436362436362436_i8
        z8=1066149217761810_i8
        w8=123456123456123456_i8
      ENDIF

    END SUBROUTINE kiss_reset


    SUBROUTINE kiss_check(check_type)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_check  ***
      !!
      !! ** Purpose :   Check the KISS pseudo-random sequence
      !!
      !! ** Method  :   Check that it reproduces the correct sequence
      !!                from the default seed
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: iter, niter, correct, iran
      CHARACTER(LEN=*) :: check_type
      LOGICAL :: print_success

      IF (kiss_32bits) STOP 'Cannot check 32-bit kiss'

      ! Save current state of KISS
      CALL kiss_save()
      ! Reset the default seed
      CALL kiss_reset()

      ! Select check type
      SELECT CASE(check_type)
      CASE('short')
         niter = 5_i8
         correct = 542381058189297533_i8
         print_success = .FALSE.
      CASE('long')
         niter = 100000000_i8
         correct = 1666297717051644203_i8 ! Check provided by G. Marsaglia
         print_success = .TRUE.
      CASE('default')
      CASE DEFAULT
         STOP 'Bad check type in kiss_check'
      END SELECT

      ! Run kiss for the required number of iterations (niter)
      DO iter=1,niter
         iran = kiss64()
      ENDDO

      ! Check that last iterate is correct
      IF (iran.NE.correct) THEN
         STOP 'Check failed: KISS internal error !!'
      ELSE
         IF (print_success) PRINT *, 'Check successful: 100 million calls to KISS OK'
      ENDIF

      ! Reload the previous state of KISS
      CALL kiss_load()

   END SUBROUTINE kiss_check


   SUBROUTINE kiss_save
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_save  ***
      !!
      !! ** Purpose :   Save current state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE

      OPEN(UNIT=30,FILE='.kiss_restart')
      IF (kiss_32bits) THEN
        WRITE(30,*) x4
        WRITE(30,*) y4
        WRITE(30,*) z4
        WRITE(30,*) w4
      ELSE
        WRITE(30,*) x8
        WRITE(30,*) y8
        WRITE(30,*) z8
        WRITE(30,*) w8
      ENDIF
      CLOSE(30)

   END SUBROUTINE kiss_save


   SUBROUTINE kiss_load
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_load  ***
      !!
      !! ** Purpose :   Load the saved state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: filexists

      INQUIRE(FILE='.kiss_restart',EXIST=filexists)
      IF (filexists) THEN
         OPEN(UNIT=30,FILE='.kiss_restart')
         IF (kiss_32bits) THEN
           READ(30,*) x4
           READ(30,*) y4
           READ(30,*) z4
           READ(30,*) w4
         ELSE
           READ(30,*) x8
           READ(30,*) y8
           READ(30,*) z8
           READ(30,*) w8
         ENDIF
         CLOSE(30)
      ENDIF

   END SUBROUTINE kiss_load


   FUNCTION kiss_uniform()
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_uniform  ***
      !!
      !! ** Purpose :   Real random numbers with uniform distribution in [0,1]
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: kiss_uniform

      IF (kiss_32bits) THEN
        kiss_uniform = REAL(kiss32(),wp) * scaling_uni2_4
      ELSE
        kiss_uniform = half + REAL(kiss64(),wp) * scaling_uni1_8
      ENDIF

   END FUNCTION kiss_uniform


   FUNCTION kiss_normal()
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_normal  ***
      !!
      !! ** Purpose :   Real random numbers with Gaussian distribution N(0,1)
      !!
      !! ** Method  :   Generate 2 new Gaussian draws (gran1 and gran2)
      !!                from 2 uniform draws on [-1,1] (u1 and u2),
      !!                using the Marsaglia polar method
      !!                (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: kiss_normal
      REAL(KIND=wp) :: u1, u2, rsq, fac

      IF (ig.EQ.1) THEN
         rsq = two
         DO WHILE ( (rsq.GE.one).OR. (rsq.EQ.zero) )
            IF (kiss_32bits) THEN
              u1 = REAL(kiss32(),wp) * scaling_uni2_4
              u2 = REAL(kiss32(),wp) * scaling_uni2_4
            ELSE
              u1 = REAL(kiss64(),wp) * scaling_uni2_8
              u2 = REAL(kiss64(),wp) * scaling_uni2_8
            ENDIF
            rsq = u1*u1 + u2*u2
         ENDDO
         fac = SQRT(-two*LOG(rsq)/rsq)
         gran1 = u1 * fac
         gran2 = u2 * fac
      ENDIF

      ! Output one of the 2 draws
      IF (ig.EQ.1) THEN
         kiss_normal = gran1 ; ig = 2
      ELSE
         kiss_normal = gran2 ; ig = 1
      ENDIF

   END FUNCTION kiss_normal


   FUNCTION kiss_gamma(k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_gamma  ***
      !!
      !! ** Purpose :   Real random numbers with Gamma distribution Gamma(k,1)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: kiss_gamma
      REAL(KIND=wp), PARAMETER :: p1 = 4.5_wp
      REAL(KIND=wp), PARAMETER :: p2 = 2.50407739677627_wp  ! 1+LOG(9/2)
      REAL(KIND=wp), PARAMETER :: p3 = 1.38629436111989_wp  ! LOG(4)
      REAL(KIND=wp) :: k, u1, u2, b, c, d, xx, yy, zz, rr, ee
      LOGICAL :: accepted

      IF (k.GT.one) THEN
         ! Cheng's rejection algorithm
         ! (see Devroye, Non-Uniform Random Variate Generation, p. 413)
         b = k - p3 ; d = SQRT(two*k-one) ; c = k + d

         accepted=.FALSE.
         DO WHILE (.NOT.accepted)
            u1 = kiss_uniform()
            yy = LOG(u1/(one-u1)) / d  ! Mistake in Devroye: "* k" instead of "/ d"
            xx = k * EXP(yy)
            rr = b + c * yy - xx
            u2 = kiss_uniform()
            zz = u1 * u1 * u2

            accepted = rr .GE. (zz*p1-p2)
            IF (.NOT.accepted) accepted =  rr .GE. LOG(zz)
         ENDDO

         kiss_gamma = xx

      ELSEIF (k.LT.one) THEN
        ! Rejection from the Weibull density
        ! (see Devroye, Non-Uniform Random Variate Generation, p. 415)
        c = one/k ; d = (one-k) * EXP( (k/(one-k)) * LOG(k) )

        accepted=.FALSE.
        DO WHILE (.NOT.accepted)
           u1 = kiss_uniform()
           zz = -LOG(u1)
           xx = EXP( c * LOG(zz) )
           u2 = kiss_uniform()
           ee = -LOG(u2)

           accepted = (zz+ee) .GE. (d+xx)  ! Mistake in Devroye: "LE" instead of "GE"
        ENDDO

        kiss_gamma = xx

      ELSE
         ! Exponential distribution
         u1 = kiss_uniform()
         kiss_gamma = -LOG(u1)

      ENDIF

   END FUNCTION kiss_gamma


   SUBROUTINE kiss_sample(a,n,k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_sample  ***
      !!
      !! ** Purpose :   Select a random sample of size k from a set of n integers
      !!
      !! ** Method  :   The sample is output in the first k elements of a
      !!                Set k equal to n to obtain a random permutation
      !!                  of the whole set of integers
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8), DIMENSION(:) :: a
      INTEGER(KIND=i8) :: n, k, i, j, atmp
      REAL(KIND=wp) :: uran

      ! Select the sample using the swapping method
      ! (see Devroye, Non-Uniform Random Variate Generation, p. 612)
      DO i=1,k
         ! Randomly select the swapping element between i and n (inclusive)
         uran = kiss_uniform()
         j = i - 1 + CEILING( REAL(n-i+1,8) * uran )
         ! Swap elements i and j
         atmp = a(i) ; a(i) = a(j) ; a(j) = atmp
      ENDDO

   END SUBROUTINE kiss_sample
END MODULE storng_kiss
