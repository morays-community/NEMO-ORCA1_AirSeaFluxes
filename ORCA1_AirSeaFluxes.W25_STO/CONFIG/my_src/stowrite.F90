MODULE stowrite
   !!======================================================================
   !!                       ***  MODULE  stowrite  ***
   !! Purpose        : write stochastic fields in NetCDF file
   !!=====================================================================
   !!   sto_write       : write current stochastic fields in NetCDF file
   !!   sto_write_init  : initialization of the NetCDF file
   !!   sto_write_final : finalize the writing of the file
   !!----------------------------------------------------------------------
   USE stoexternal, only : jpi, jpj, jpk
   USE stoarray
   USE netcdf

   IMPLICIT NONE
   PRIVATE

   ! index of NetCDF file
   INTEGER, SAVE :: ncid
   ! index of dimensions in NetCDF file
   INTEGER, SAVE :: idx, idy, idz, idt
   ! index of stochastic arrays in NetCDF file
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: idsto

   PUBLIC sto_write, sto_write_init, sto_write_final

CONTAINS

   SUBROUTINE sto_write(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_write  ***
      !!
      !! ** Purpose :  write current stochastic fields in NetCDF file
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      INTEGER :: is, jsto
      INTEGER, DIMENSION(4) :: vstart

      ! Write stochastic fields
      print *, 'kt =',kt
      DO jsto=1,jpsto
         vstart = 1
         IF (stofields(jsto)%dim == 2) THEN
            vstart(3) = kt
            is = NF90_PUT_VAR(ncid,idsto(jsto),stofields(jsto)%sto2d,start=vstart)
         ELSEIF (stofields(jsto)%dim == 3) THEN
            vstart(4) = kt
            is = NF90_PUT_VAR(ncid,idsto(jsto),stofields(jsto)%sto3d,start=vstart)
         ELSEIF (stofields(jsto)%dim == 0) THEN
            vstart(1) = kt
            is = NF90_PUT_VAR(ncid,idsto(jsto),stofields(jsto)%sto0d,start=vstart)
         ENDIF
      ENDDO

   END SUBROUTINE sto_write


   SUBROUTINE sto_write_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_write_init  ***
      !!
      !! ** Purpose :   initialization of the NetCDF file
      !!----------------------------------------------------------------------
      INTEGER :: is, jsto

      ALLOCATE(idsto(jpsto))

      ! Create output filr
      is = NF90_CREATE("stofields.nc",NF90_CLOBBER,ncid)

      ! Define dimensions
      is = NF90_DEF_DIM(ncid,'x',jpi,idx)
      is = NF90_DEF_DIM(ncid,'y',jpj,idy)
      is = NF90_DEF_DIM(ncid,'z',jpk,idz)
      is = NF90_DEF_DIM(ncid,'t',nf90_unlimited,idt)

      ! Define variables
      DO jsto=1,jpsto
         IF (stofields(jsto)%dim == 2) THEN
            is = NF90_DEF_VAR(ncid,stofields(jsto)%stoname, &
            &              NF90_FLOAT,(/idx,idy,idt/),idsto(jsto))
         ELSEIF (stofields(jsto)%dim == 3) THEN
            is = NF90_DEF_VAR(ncid,stofields(jsto)%stoname, &
            &              NF90_FLOAT,(/idx,idy,idz,idt/),idsto(jsto))
         ELSEIF (stofields(jsto)%dim == 0) THEN
            is = NF90_DEF_VAR(ncid,stofields(jsto)%stoname, &
            &              NF90_FLOAT,(/idt/),idsto(jsto))
         ENDIF
      ENDDO

      is = NF90_ENDDEF(ncid)
      
   END SUBROUTINE sto_write_init


   SUBROUTINE sto_write_final
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_write_final  ***
      !!
      !! ** Purpose :   finalize the writing of the file
      !!----------------------------------------------------------------------
      INTEGER :: is

      is = NF90_CLOSE(ncid)

   END SUBROUTINE sto_write_final

END MODULE stowrite
