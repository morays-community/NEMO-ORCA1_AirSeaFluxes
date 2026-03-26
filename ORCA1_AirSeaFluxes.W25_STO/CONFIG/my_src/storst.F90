MODULE storst
   !!======================================================================
   !!                       ***  MODULE  storst  ***
   !! Purpose        : read and write restart files for stochastic fields
   !!=====================================================================
   !!   sto_rst_read  : read restart file
   !!   sto_rst_write : write restart file
   !!----------------------------------------------------------------------
   USE netcdf
   USE stoarray         ! stochastic arrays to store in restart
   USE stowhite         ! white noise generator
   USE storng_kiss      ! KISS random number generator
   USE storng_ziggurat  ! Ziggurat algorithm to generate normal numbers
   USE stoexternal, only : wp, lc, jpi, jpj, jpk, cn_member, narea, ln_ensemble, lbc_lnk

   ! NEMO specific modules
   USE iom             ! I/O module

   IMPLICIT NONE
   PRIVATE

   PUBLIC sto_rst_read, sto_rst_write

CONTAINS

   SUBROUTINE sto_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_read  ***
      !!
      !! ** Purpose :   read stochastic parameters from restart file
      !!----------------------------------------------------------------------
      ! Variables to describe seed of random number generator
      INTEGER(KIND=8) :: ziseed8(4) ! RNG seeds in 64-bit integer type
      INTEGER(KIND=4) :: ziseed4(4) ! RNG seeds in 32-bit integer type
      REAL(KIND=8)    :: zrseed8(4) ! RNG seeds in 64-bit real type (with same bits to save in restart)
      REAL(KIND=4)    :: zrseed4(4) ! RNG seeds in 32-bit real type (with same bits to save in restart)
      ! Other variables
      INTEGER             ::   jsto, jidx, jseed
      INTEGER             ::   idg                   ! number of digits
      CHARACTER(LEN=11)   ::   clsto2d='sto2d_000_0' ! stochastic parameter variable name
      CHARACTER(LEN=11)   ::   clsto3d='sto3d_000_0' ! stochastic parameter variable name
      CHARACTER(LEN=15)   ::   clseed='seed0_0000'   ! seed variable name
      CHARACTER(LEN=6)    ::   clfmt                 ! writing format
      ! Type of stochastic parameter
      CHARACTER(len=1) :: sto_typ
      REAL(wp) :: sto_sgn

      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sto_rst_read : read stochastic parameters from restart file'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF

         ! Open the restart file
         CALL iom_open( cn_storst_in, numstor )

         ! Get stochastic parameters from restart file:
         ! 2D stochastic parameters
         DO jsto = 1 , jpsto2d
            sto_typ = sto2d_typ(jsto)
            sto_sgn = sto2d_sgn(jsto)
            DO jidx = 1, jpidx2d + jpidxsup2d
               WRITE(clsto2d(7:9),'(i3.3)') jsto
               WRITE(clsto2d(11:11),'(i1.1)') jidx
               CALL iom_get( numstor, jpdom_auto, clsto2d, sto2d(:,:,jidx,jsto), cd_type = sto_typ, psgn = sto_sgn )
            END DO
         END DO
         ! 3D stochastic parameters
         DO jsto = 1 , jpsto3d
            sto_typ = sto3d_typ(jsto)
            sto_sgn = sto3d_sgn(jsto)
            DO jidx = 1, jpidx3d + jpidxsup3d
               WRITE(clsto3d(7:9),'(i3.3)') jsto
               WRITE(clsto3d(11:11),'(i1.1)') jidx
               CALL iom_get( numstor, jpdom_auto, clsto3d, sto3d(:,:,:,jidx,jsto), cd_type = sto_typ, psgn = sto_sgn )
            END DO
         END DO

         IF (ln_rstseed) THEN
            ! Get saved state of the random number generator
            idg = MAX( INT(LOG10(REAL(jpnij,wp))) + 1, 4 )        ! how many digits to we need to write? min=4, max=9
            WRITE(clfmt, "('(i', i1, '.', i1, ')')") idg, idg     ! "(ix.x)"
            DO jseed = 1 , 4
               WRITE(clseed(5:5)      ,'(i1.1)') jseed
               WRITE(clseed(7:7+idg-1),  clfmt ) narea
               CALL iom_get( numstor, clseed(1:7+idg-1) , zrseed8(jseed) )
            END DO
            SELECT CASE (c_rngtype)
            CASE('kiss64')
              ziseed8 = TRANSFER( zrseed8 , ziseed8)
              CALL kiss_seed( ziseed8(1) , ziseed8(2) , ziseed8(3) , ziseed8(4) )
            CASE('kiss32')
              ziseed8 = TRANSFER( zrseed8 , ziseed8)
              ziseed4 = ziseed8
              CALL kiss_seed( ziseed4(1) , ziseed4(2) , ziseed4(3) , ziseed4(4) )
            CASE('shr3')
              ziseed8 = TRANSFER( zrseed8 , ziseed8)
              ziseed4 = ziseed8
              CALL shr3_seed( ziseed4(1) )
            END SELECT
         ENDIF

         ! Close the restart file
         CALL iom_close( numstor )

      ENDIF

   END SUBROUTINE sto_rst_read


   SUBROUTINE sto_rst_write( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_read  ***
      !!
      !! ** Purpose :   write stochastic parameters in restart file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      ! Variables to describe seed of random number generator
      INTEGER(KIND=8) :: ziseed8(4) ! RNG seeds in 64-bit integer type
      INTEGER(KIND=4) :: ziseed4(4) ! RNG seeds in 32-bit integer type
      REAL(KIND=8)    :: zrseed8(4) ! RNG seeds in 64-bit real type (with same bits to save in restart)
      REAL(KIND=4)    :: zrseed4(4) ! RNG seeds in 32-bit real type (with same bits to save in restart)
      ! Other variables
      INTEGER             ::   jsto, jidx, jseed
      INTEGER             ::   idg                   ! number of digits
      CHARACTER(LEN=20)   ::   clkt                  ! ocean time-step defined as a character
      CHARACTER(LEN=50)   ::   clname                ! restart file name
      CHARACTER(LEN=11)   ::   clsto2d='sto2d_000_0' ! stochastic parameter variable name
      CHARACTER(LEN=11)   ::   clsto3d='sto3d_000_0' ! stochastic parameter variable name
      CHARACTER(LEN=15)   ::   clseed='seed0_0000'   ! seed variable name
      CHARACTER(LEN=6)    ::   clfmt                 ! writing format
      !!----------------------------------------------------------------------

      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN

         IF( kt == nitrst .OR. kt == nitend ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'sto_rst_write : write stochastic parameters in restart file'
               WRITE(numout,*) '~~~~~~~~~~~~~'
            ENDIF
         ENDIF

         ! Put stochastic parameters in restart files
         ! (as opened at previous timestep, see below)
         IF( kt > nit000) THEN
            IF( kt == nitrst .OR. kt == nitend ) THEN
   
               ! Get type and state of random number generator to save in restart file
               SELECT CASE (c_rngtype)
               CASE('kiss64')
                  CALL kiss_state( ziseed8(1) , ziseed8(2) , ziseed8(3) , ziseed8(4) )
                  zrseed8 = TRANSFER( ziseed8 , zrseed8 )
               CASE('kiss32')
                  CALL kiss_state( ziseed4(1) , ziseed4(2) , ziseed4(3) , ziseed4(4) )
                  ziseed8 = ziseed4
                  zrseed8 = TRANSFER( ziseed8 , zrseed8 )
               CASE('shr3')
                  CALL shr3_state( ziseed4(1) )
                  ziseed8 = ziseed4
                  zrseed8 = TRANSFER( ziseed8 , zrseed8 )
               CASE DEFAULT
                  STOP 'Bad type of random number generator in storst'
               END SELECT
      
               ! get and save current state of the random number generator
               idg = MAX( INT(LOG10(REAL(jpnij,wp))) + 1, 4 )    ! how many digits to we need to write? min=4, max=9
               WRITE(clfmt, "('(i', i1, '.', i1, ')')") idg, idg ! "(ix.x)"
               DO jseed = 1 , 4
                  WRITE(clseed(5:5)      ,'(i1.1)') jseed
                  WRITE(clseed(7:7+idg-1),  clfmt ) narea
                  CALL iom_rstput( kt, nitrst, numstow, clseed(1:7+idg-1), zrseed8(jseed) )
               END DO
               ! 2D stochastic parameters
               DO jsto = 1 , jpsto2d
                  DO jidx = 1, jpidx2d + jpidxsup2d
                     WRITE(clsto2d(7:9),'(i3.3)') jsto
                     WRITE(clsto2d(11:11),'(i1.1)') jidx
                     CALL iom_rstput( kt, nitrst, numstow, clsto2d , sto2d(:,:,jidx,jsto) )
                  END DO
               END DO
               ! 3D stochastic parameters
               DO jsto = 1 , jpsto3d
                  DO jidx = 1, jpidx3d + jpidxsup3d
                     WRITE(clsto3d(7:9),'(i3.3)') jsto
                     WRITE(clsto3d(11:11),'(i1.1)') jidx
                     CALL iom_rstput( kt, nitrst, numstow, clsto3d , sto3d(:,:,:,jidx,jsto) )
                  END DO
               END DO
               ! close the restart file
               CALL iom_close( numstow )
   
            ENDIF
         ENDIF

         ! Open the restart file one timestep before writing restart
         IF( kt < nitend) THEN
            IF( kt == nitrst - 1 .OR. nn_stock == 1 .OR. kt == nitend-1 ) THEN
               ! create the filename
               IF( nitrst > 999999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
               ELSE                            ;   WRITE(clkt, '(i8.8)') nitrst
               ENDIF
               clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_storst_out)
               ! print information
               IF(lwp) THEN
                  WRITE(numout,*) '             open stochastic parameters restart file: '//clname
                  IF( kt == nitrst - 1 ) THEN
                     WRITE(numout,*) '             kt = nitrst - 1 = ', kt
                  ELSE
                     WRITE(numout,*) '             kt = '             , kt
                  ENDIF
               ENDIF
               ! open the restart file
               CALL iom_open( clname, numstow, ldwrt = .TRUE. )
            ENDIF
         ENDIF

      ENDIF

   END SUBROUTINE sto_rst_write

END MODULE storst
