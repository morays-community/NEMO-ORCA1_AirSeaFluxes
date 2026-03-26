MODULE storng_ziggurat
   !!======================================================================
   !!                       ***  MODULE  storng_ziggurat  ***
   !! Random number generator, used in NEMO stochastic parameterization
   !!
   !!=====================================================================
   !! History :  4.0  ! 2023-10 (A. Storto)       Original code
   !!                 ! 2024-01 (J.-M. Brankart)  Include option to use kiss
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! The module applies an efficient method (Ziggurat, Marsaglia and Tsang, 2000)
   !! to compute normal random numbers from a random sequence of 32-bit integers
   !!
   !! It also include a simple and fast 32-bit random generator (shr3 function)
   !!----------------------------------------------------------------------
   !!
   !! Marsaglia & Tsang generator for random normals & random exponentials.
   !! Translated from C by Alan Miller (amiller@bigpond.net.au)
   !!
   !! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
   !! random variables', J. Statist. Software, v5(8).
   !!
   !! This is an electronic journal which can be downloaded from:
   !! http://www.jstatsoft.org/v05/i08
   !!
   !! N.B. It is assumed that all integers are 32-bit.
   !! N.B. The value of M2 has been halved to compensate for the lack of
   !!      unsigned integers in Fortran.
   !!
   !! Latest version - 1 January 2001
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zig_set     : initialize the Ziggurat method
   !!   zig_normal  : get next normal random number with Ziggurat
   !!   zig_exp     : get next exponential random number with Ziggurat
   !!   shr3        : simple 32-bit random number generator (period = 2^32-1)
   !!   shr3_seed   : seed shr3 random number generator
   !!   shr3_reset  : reset shr3 to default seed
   !!   shr3_uni    : get next uniform random number using shr3
   !!   shr3_normal : get next normal random number using shr3 + simple method
   !!----------------------------------------------------------------------
   USE stoexternal, only : sp, dp, wp, i4, i8
   USE storng_kiss

   IMPLICIT NONE
   PRIVATE

   ! Define working precision: single (sp) or double (dp) precision
   INTEGER, PARAMETER :: rp = sp

   ! Public parameter defining the random number generator to use
   CHARACTER(len=6), SAVE, PUBLIC :: zig_rngtype='kiss32'
   INTEGER, SAVE :: rngtype
   LOGICAL, SAVE ::  initialized=.FALSE.

   ! Parameters for the Ziggurat method
   REAL(rp), PARAMETER    ::  m1=2147483648.0_rp,   m2=2147483648.0_rp,      &
                              half=0.5_rp
   REAL(rp)               ::  dn=3.442619855899_rp, tn=3.442619855899_rp,    &
                              vn=0.00991256303526217_rp,                     &
                              q,                    de=7.697117470131487_rp, &
                              te=7.697117470131487_rp,                       &
                              ve=0.003949659822581572_rp
   INTEGER(KIND=i4), SAVE ::  iz, jz, kn(0:127), ke(0:255), hz
   REAL(rp), SAVE         ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)

   ! Parameters for the shr3 random number generator
   INTEGER(KIND=i4), SAVE ::  jsr=123456789  ! default seed
   REAL(rp), PARAMETER    ::  scaling_uni1 = 0.5_rp/(real(huge(jsr),rp)+1.0_rp)
   REAL(rp), PARAMETER    ::  scaling_uni2 = 1.0_rp/(real(huge(jsr),rp)+1.0_rp)

   ! Variables to store 2 Gaussian random numbers with current index (ig)
   INTEGER(KIND=i8), SAVE :: ig=1
   REAL(KIND=wp), SAVE :: gran1, gran2

   ! Public routines
   PUBLIC  :: zig_set, zig_normal, zig_exp
   PUBLIC  :: shr3, shr3_seed, shr3_reset, shr3_uni, shr3_normal, shr3_state

CONTAINS

   SUBROUTINE zig_set()
      !! --------------------------------------------------------------------
      !!                  ***  SUBROUTINE zig_set ***
      !!
      !! ** Purpose : Initialize the Ziggurat method
      !!
      !! ** Method  : Set precomputed tables of parameters  
      !! 
      ! --------------------------------------------------------------------
      INTEGER  :: i

      ! Define type oof random number generator to use
      SELECT CASE(zig_rngtype)
      CASE('shr3')  ! shr3 random number generator
        rngtype = 0
      CASE('kiss32')  ! 32-bit kiss random number generator
        kiss_32bits=.TRUE.
        rngtype = 1
      CASE('kiss64')  ! 64-bit kiss random number generator
        kiss_32bits=.FALSE.
        rngtype = 2
      CASE DEFAULT
        STOP 'Bad random number generator type in storng_ziggurat'
      END SELECT

      ! Reinitialize parameters if zig_set was already called
      IF (initialized) THEN
        dn=3.442619855899_rp
        tn=3.442619855899_rp
        de=7.697117470131487_rp
        te=7.697117470131487_rp
      ENDIF

      !  Tables for RNOR
      q = vn*EXP(half*dn*dn)
      kn(0) = (dn/q)*m1
      kn(1) = 0
      wn(0) = q/m1
      wn(127) = dn/m1
      fn(0) = 1.0_rp
      fn(127) = EXP( -half*dn*dn )
      DO  i = 126, 1, -1
         dn = SQRT( -2.0_rp * LOG( vn/dn + EXP( -half*dn*dn ) ) )
         kn(i+1) = (dn/tn)*m1
         tn = dn
         fn(i) = EXP(-half*dn*dn)
         wn(i) = dn/m1
      END DO

      !  Tables for REXP
      q = ve*EXP( de )
      ke(0) = (de/q)*m2
      ke(1) = 0
      we(0) = q/m2
      we(255) = de/m2
      fe(0) = 1.0_rp
      fe(255) = EXP( -de )
      DO  i = 254, 1, -1
         de = -LOG( ve/de + EXP( -de ) )
         ke(i+1) = m2 * (de/te)
         te = de
         fe(i) = EXP( -de )
         we(i) = de/m2
      END DO

      initialized = .TRUE.

      RETURN
   END SUBROUTINE zig_set


   FUNCTION zig_normal( ) RESULT( fn_val )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION zig_normal ***
      !!
      !! ** Purpose : get next normal random number
      !!
      !! ** Method  : apply the Ziggurat method
      !!
      ! --------------------------------------------------------------------
      REAL(rp)             ::  fn_val

      REAL(rp), PARAMETER  ::  r = 3.442620_rp
      REAL(rp)             ::  x, y

      IF( .NOT. initialized ) CALL zig_set()

      hz = rng_i32( )
      iz = IAND( hz, 127 )
      IF( ABS( hz ) < kn(iz) ) THEN
         fn_val = hz * wn(iz)
      ELSE
         DO
            IF( iz == 0 ) THEN
               DO
                  x = -0.2904764_rp * LOG( rng_uni( ) )
                  y = -LOG( rng_uni( ) )
                  IF( y+y >= x*x ) EXIT
               END DO
               fn_val = r+x
               IF( hz <= 0 ) fn_val = -fn_val
               RETURN
            END IF
            x = hz * wn(iz)
            IF( fn(iz) + rng_uni( )*(fn(iz-1)-fn(iz)) < EXP(-half*x*x) ) THEN
               fn_val = x
               RETURN
            END IF
            hz = rng_i32( )
            iz = IAND( hz, 127 )
            IF( ABS( hz ) < kn(iz) ) THEN
               fn_val = hz * wn(iz)
               RETURN
            END IF
         END DO
      END IF
      RETURN
   END FUNCTION zig_normal


   FUNCTION zig_exp( ) RESULT( fn_val )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION zig_exp ***
      !!
      !! ** Purpose : get next exponential random number
      !!
      !! ** Method  : apply the Ziggurat method
      !!
      ! --------------------------------------------------------------------
      REAL(rp)  ::  fn_val

      REAL(rp)  ::  x

      IF( .NOT. initialized ) CALL zig_set()

      jz = rng_i32( )
      iz = IAND( jz, 255 )
      IF( ABS( jz ) < ke(iz) ) THEN
         fn_val = ABS(jz) * we(iz)
         RETURN
      END IF
      DO
         IF( iz == 0 ) THEN
            fn_val = 7.69711 - LOG( rng_uni( ) )
            RETURN
         END IF
         x = ABS( jz ) * we(iz)
         IF( fe(iz) + rng_uni( )*(fe(iz-1) - fe(iz)) < EXP( -x ) ) THEN
            fn_val = x
            RETURN
         END IF
         jz = rng_i32( )
         iz = IAND( jz, 255 )
         IF( ABS( jz ) < ke(iz) ) THEN
            fn_val = ABS( jz ) * we(iz)
            RETURN
         END IF
      END DO
      RETURN
   END FUNCTION zig_exp


   FUNCTION rng_i32( ) RESULT( ival )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION rng_i32  ***
      !!
      !! ** Purpose :   select random number generator (i32 output)
      !!
      ! --------------------------------------------------------------------
      INTEGER(KIND=i4) ::  ival

      IF (rngtype==0) THEN
         ival = shr3( )
      ELSEIF (rngtype==1) THEN
         ival = kiss32( )
         ival = IOR(ival,ISHFT(ISHFT(ival,-8),31)) ! use 9th bit as sign bit
      ELSEIF (rngtype==2) THEN
         ival = ISHFT( kiss64() , -32 )            ! use last 32 bits
      ENDIF

      RETURN
   END FUNCTION rng_i32


   FUNCTION rng_uni( ) RESULT( fn_val )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION rng_i32  ***
      !!
      !! ** Purpose :   select random number generator (uniform output)
      !!
      ! --------------------------------------------------------------------
      REAL(rp)  ::  fn_val

      IF (rngtype==0) THEN
         fn_val = shr3_uni( )
      ELSEIF (rngtype==1) THEN
         fn_val = kiss_uniform()
      ELSEIF (rngtype==2) THEN
         fn_val = kiss_uniform()
      ENDIF

      RETURN
   END FUNCTION rng_uni


   FUNCTION shr3( ) RESULT( ival )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION shr3  ***
      !!
      !! ** Purpose :   fast and simple 32-bit random number generator
      !!
      !! ** Method  :   Xorshift (XSH), period 2^32-1
      !!
      ! --------------------------------------------------------------------

      INTEGER(KIND=i4) ::  ival
   
      jz = jsr
      jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
      jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
      jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
      ival = jz + jsr
      RETURN
   END FUNCTION shr3


   SUBROUTINE shr3_seed( jsrseed )
      !! --------------------------------------------------------------------
      !!                  ***  SUBROUTINE shr3_seed ***
      !!
      !! ** Purpose : seed shr3 random number genrator
      !!
      !! ** Method  : provide 32-bit integer
      !!
      ! --------------------------------------------------------------------
      INTEGER(KIND=i4), INTENT(IN)  :: jsrseed

      jsr = jsrseed

   END SUBROUTINE shr3_seed


   SUBROUTINE shr3_state( jsrstate )
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE shr3_state  ***
      !!
      !! ** Purpose :   Get current state of shr3 random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i4) :: jsrstate

      jsrstate = jsr

   END SUBROUTINE shr3_state


   SUBROUTINE shr3_reset( )
      !!                  ***  SUBROUTINE shr3_reset ***
      !!
      !! ** Purpose : reset shr3 to default seed
      !!
      ! --------------------------------------------------------------------
      jsr = 123456789

   END SUBROUTINE shr3_reset

   FUNCTION shr3_uni( ) RESULT( fn_val )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION shr3_uni  ***
      !!
      !! ** Purpose : get uniform random number using shr3
      !!
      ! --------------------------------------------------------------------
      REAL(rp)  ::  fn_val
   
      fn_val = half + scaling_uni1 * shr3( )
      RETURN
   END FUNCTION shr3_uni


   FUNCTION shr3_normal( ) RESULT( fn_val )
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION shr3_normal  ***
      !!
      !! ** Purpose :   get normal random number using shr3
      !!
      !! ** Method  :   generate 2 new normal draws (gran1 and gran2)
      !!                from 2 uniform draws on [-1,1] (u1 and u2),
      !!                using the Marsaglia polar method
      !!                (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rp)  ::  fn_val
      REAL(KIND=rp) :: u1, u2, rsq, fac

      IF (ig.EQ.1) THEN
         rsq = 2.0_rp
         DO WHILE ( (rsq.GE.1.0_rp).OR. (rsq.EQ.0.0_rp) )
            u1 = REAL(shr3(),rp) * scaling_uni2
            u2 = REAL(shr3(),rp) * scaling_uni2
            rsq = u1*u1 + u2*u2
         ENDDO
         fac = SQRT(-2.0_rp*LOG(rsq)/rsq)
         gran1 = u1 * fac
         gran2 = u2 * fac
      ENDIF

      ! Output one of the 2 draws
      IF (ig.EQ.1) THEN
         fn_val = gran1 ; ig = 2
      ELSE
         fn_val = gran2 ; ig = 1
      ENDIF

      RETURN
   END FUNCTION shr3_normal

END MODULE storng_ziggurat
