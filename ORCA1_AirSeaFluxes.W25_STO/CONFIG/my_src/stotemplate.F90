MODULE stotemplate
   !!======================================================================
   !!                       ***  MODULE stotemplate  ***
   !!
   !! Purpose : Template module for a new stochastic scheme
   !!           Illustrating how to use the stochastic modules
   !!======================================================================
   USE stoarray

   IMPLICIT NONE
   PRIVATE

   INTEGER :: jstotemplate1   ! index of stochastic field used in this scheme
   INTEGER :: jstotemplate2   ! a possible second index
   INTEGER :: jstotemplate3   ! a third, it could also be an array...
   INTEGER :: jstotemplate4   ! etc.
   INTEGER :: jstotemplate5   ! etc.
   INTEGER :: jstotemplate6   ! etc

   PUBLIC sto_template, sto_template_init

CONTAINS

   SUBROUTINE sto_template(kt)
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_template  ***
      !!
      !! This routine is called at every time step
      !! to make appropriate use of the stochastic field
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      ! Use arrays (for instance with index jstotemplate1) as:
      ! stofields(jstotemplate1)%sto2d(:,:),
      ! stofields(jstotemplate1)%sto3d(:,:,:), or
      ! stofields(jstotemplate1)%sto0d
      ! depending on the requested dimension,
      ! available in stofields(jstotemplate1)%dim

      ! Here we just write some statistics of the current stochastic field
      ! print *, 'min',MINVAL(stofields(jstotemplate1)%sto2d)
      ! print *, 'max',MAXVAL(stofields(jstotemplate1)%sto2d)

   END SUBROUTINE sto_template


   SUBROUTINE sto_template_init
      !!----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE sto_template_init  ***
      !!
      !! This routine is calle at initialization time
      !! to request stochastic field with appropriate features
      !!
      !!----------------------------------------------------------------------
      LOGICAL :: l1, l2, l3, l4, l5, l6

      l1 = .FALSE. ! exclude example 1
      l2 = .FALSE. ! exclude example 2
      l3 = .FALSE. ! exclude example 3
      l4 = .FALSE. ! exclude example 4
      l5 = .FALSE. ! exclude example 5
      l6 = .FALSE. ! exclude example 6

      l3 = .TRUE.  ! include example 3 in output file

      ! Read namelist block corresponding to this stochastic scheme
      ! -> get parameters

      ! Example 1 : white noise
      ! -----------------------
      IF (l1) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate1)
        ! No specification required
      ENDIF

      ! Example 2 : autoregressive processes
      ! ------------------------------------
      IF (l2) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate2)
        ! Set features of the requested stochastic field
        stofields(jstotemplate2)%type_t='arn'
        stofields(jstotemplate2)%corr_t=5.0
        stofields(jstotemplate2)%nar_order=2
        stofields(jstotemplate2)%nar_update=5
      ENDIF

      ! Example 3 : space correlation (with diffusive operator)
      ! -------------------------------------------------------
      IF (l3) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate3)
        ! Set features of the requested stochastic field
        stofields(jstotemplate3)%type_t='arn'
        stofields(jstotemplate3)%corr_t=5.0
        stofields(jstotemplate3)%nar_order=2
        stofields(jstotemplate3)%nar_update=5

        stofields(jstotemplate3)%type_xy='diffusive'
        stofields(jstotemplate3)%diff_passes=50
        stofields(jstotemplate3)%diff_type=1
      ENDIF

      ! Example 4 : space correlation (with kernel convolution)
      ! -------------------------------------------------------
      IF (l4) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate4)
        ! Set features of the requested stochastic field
        stofields(jstotemplate4)%type_t='arn'
        stofields(jstotemplate4)%corr_t=5.0
        stofields(jstotemplate4)%nar_order=2
        stofields(jstotemplate4)%nar_update=5

        stofields(jstotemplate4)%type_xy='kernel'
        stofields(jstotemplate4)%corr_xy=10.0
        stofields(jstotemplate4)%ker_type=0
        stofields(jstotemplate4)%ker_coord=2
      ENDIF

      ! Example 5 : modified marginal distribution: lognormal
      ! -----------------------------------------------------
      IF (l5) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate5)
        ! Set features of the requested stochastic field
        stofields(jstotemplate5)%type_t='arn'
        stofields(jstotemplate5)%corr_t=5.0
        stofields(jstotemplate5)%nar_order=2
        stofields(jstotemplate5)%nar_update=5

        stofields(jstotemplate5)%type_xy='kernel'
        stofields(jstotemplate5)%corr_xy=10.0
        stofields(jstotemplate5)%ker_type=0
        stofields(jstotemplate5)%ker_coord=2

        stofields(jstotemplate5)%type_variate='lognormal'
        stofields(jstotemplate5)%ave=1.0
        stofields(jstotemplate5)%std=0.5
      ENDIF

      ! Example 6 : modified marginal distribution: bounded
      ! ---------------------------------------------------
      IF (l6) THEN
        ! Request index for a new stochastic array
        CALL sto_array_request_new(jstotemplate6)
        ! Set features of the requested stochastic field
        stofields(jstotemplate6)%type_t='arn'
        stofields(jstotemplate6)%corr_t=5.0
        stofields(jstotemplate6)%nar_order=2
        stofields(jstotemplate6)%nar_update=5

        stofields(jstotemplate6)%type_xy='kernel'
        stofields(jstotemplate6)%corr_xy=10.0
        stofields(jstotemplate6)%ker_type=0
        stofields(jstotemplate6)%ker_coord=2

        stofields(jstotemplate6)%type_variate='bounded_atan'
        stofields(jstotemplate6)%min=0.
        stofields(jstotemplate6)%max=1.
        stofields(jstotemplate6)%ave=0.7
        stofields(jstotemplate6)%std=0.1
      ENDIF

   END SUBROUTINE sto_template_init

   !!======================================================================
END MODULE stotemplate
