! File:         module_id_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE id_base

  !***********************************************************
  !
  !# This MODULE contains the definition of TYPE idbase,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (ID) to be set up for |sphincsbssn|.
  !  That is, a binary neutron star system, a rotating
  !  star, a binary black hole system, etc.
  !
  !  PROCEDURES and variables shared by all these types
  !  of ID should belong to TYPE idbase, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific types of ID.
  !
  !  FT 24.09.2021
  !
  !***********************************************************


  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !*************************************************************
  !                                                            *
  !              Definition of TYPE idbase                     *
  !                                                            *
  !   This ABSTRACT TYPE represents a generic ID for           *
  !   |sphincsbssn| (binary neutron star, rotating star...).   *
  !                                                            *
  !*************************************************************

  TYPE, ABSTRACT:: idbase
  !# Represents a generic ID for |sphincsbssn| (binary neutron star, rotating
  !  star, etc.)


    PRIVATE


    INTEGER:: n_matter= 0
    !# Number of matter objects belonging the physical system.
    !  For example, n_objects= 2 for a binary system of stars, and n_objects= 1
    !  for a single star or for a binary system of a black hole and a star.


    LOGICAL:: one_lapse
    !! Logical variable that determines if the lapse function \(\alpha=1\),
    !  i.e., if the geodesic gauge is to be used
    LOGICAL:: zero_shift
    !! Logical variable that determines if the shift \(\beta^i=0\)


    TYPE(timer), PUBLIC:: construction_timer
    !! Timer that times the construction of the appropriate object


    CONTAINS


    !---------------------------!
    !--  DEFERRED PROCEDURES  --!
    !---------------------------!

    !
    !-- PROCEDURES to read the value of a field at a point
    !

    PROCEDURE(read_double_at_pos),        DEFERRED:: read_mass_density
    !# Returns the baryon mass density at the given point

    PROCEDURE(read_integer_at_pos),       DEFERRED:: test_position
    !# Returns 1 if the position has physically acceptable properties,
    !  0 otherwise

    !
    !-- PROCEDURES to read the value of several fields at several points
    !

    PROCEDURE(read_id_full_int),          DEFERRED:: read_id_full
    !# Reads the full ID

    PROCEDURE(read_id_particles_int),     DEFERRED:: read_id_particles
    !! Reads the hydro ID needed to compute the SPH ID

    PROCEDURE(read_id_mass_b_int),        DEFERRED:: read_id_mass_b
    !! Reads the hydro ID needed to compute the baryon mass

    PROCEDURE(read_id_spacetime_int),     DEFERRED:: read_id_spacetime
    !# Reads the spacetime ID needed to compute
    !  the BSSN variables and constraints

    PROCEDURE(read_id_hydro_int),         DEFERRED:: read_id_hydro
    !# Reads the hydro ID needed to compute the constraints on the refined mesh

    PROCEDURE(read_id_k_int),             DEFERRED:: read_id_k
    !! Reads the components of the extrinsic curvature

    !
    !-- PROCEDURES returning the values of some parameters of a matter object
    !

    PROCEDURE(return_spatial_extent_int), DEFERRED:: return_spatial_extent
    !# Returns the spatial extent of the matter objects,
    !  returning the array of 6 numbers
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)

    PROCEDURE(return_double_parameter),   DEFERRED:: return_mass
    !! Returns the masses of the matter objects.

    PROCEDURE(return_position),           DEFERRED:: return_center
    !! Returns the centers of the matter objects.

    PROCEDURE(return_position),           DEFERRED:: return_barycenter
    !! Returns the barycenters (centers of mass) of the matter objects.

    PROCEDURE(return_eos_parameters_int), DEFERRED:: return_eos_parameters
    !# Returns the identification number of the |eos| of the matter objects.
    !  @todo Set up a convention for the identification number

    PROCEDURE(return_string_parameter),   DEFERRED:: return_eos_name
    !! Returns the name of the |eos| of the matter objects.

    !
    !-- PROCEDURES returning the values of some parameters of the entire system
    !

 !   PROCEDURE(return_total_spatial_extent_int), DEFERRED:: &
 !                                                   return_total_spatial_extent
    !# Returns the spatial extent of the physical system considered,
    !  as the array of 6 numbers
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)


    !PROCEDURE(derived_type_constructor_int), NOPASS, DEFERRED:: derived_type_constructor
    !#


    !-------------------------------!
    !--  NON-DEFERRED PROCEDURES  --!
    !-------------------------------!


    !PROCEDURE, NON_OVERRIDABLE:: sanity_check
    !#


    !PROCEDURE, NON_OVERRIDABLE:: initialize
    !#


    PROCEDURE, NON_OVERRIDABLE:: get_total_spatial_extent
    !# Returns the spatial extent of the physical system considered,
    !  as the array of 6 numbers
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)


    PROCEDURE, NON_OVERRIDABLE:: set_n_matter
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to a value


    PROCEDURE, NON_OVERRIDABLE:: get_n_matter
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system


    PROCEDURE, NON_OVERRIDABLE:: set_one_lapse
    !# Sets [[idbase:one_lapse]], the logical variable that determines if
    !  the lapse \(\alpha=1\), i.e., if the geodesic gauge is to be used


    PROCEDURE, NON_OVERRIDABLE:: get_one_lapse
    !# Returns [[idbase:one_lapse]], the logical variable that determines if
    ! the lapse function \(\alpha=1\), i.e., if the geodesic gauge is to be used


    PROCEDURE, NON_OVERRIDABLE:: set_zero_shift
    !# Sets [[idbase:zero_shift]], the logical variable that determines if
    !  the shift \(\beta^i=0\)


    PROCEDURE, NON_OVERRIDABLE:: get_zero_shift
    !# Returns [[idbase:zero_shift]], the logical variable that determines if
    !  the shift \(\beta^i=0\)


    PROCEDURE, NON_OVERRIDABLE:: check_i_matter
    !# Checks that the given index is between 1 and [[idbase:n_matter]],
    !  included. If not, it stops the execution of the program.


    PROCEDURE:: integrate_baryon_mass_density
    !# Integrates the baryon mass density over a matter object, using spherical
    !  coordinates, and computes its radial profile inside the star


  END TYPE idbase


  ABSTRACT INTERFACE


   ! FUNCTION derived_type_constructor_int() RESULT( derived_type )
   ! !#
   !
   !   IMPORT:: idbase
   !   CLASS(idbase), POINTER:: derived_type
   !
   ! END FUNCTION derived_type_constructor_int


    FUNCTION read_double_at_pos( THIS, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION at a given
    !  position

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION:: res
      !! Real number at \((x,y,z)\)

    END FUNCTION read_double_at_pos


    FUNCTION read_integer_at_pos( THIS, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER at a given position

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !! \(z\) coordinate of the desired point
      INTEGER:: res
      !! Integer at \((x,y,z)\)

    END FUNCTION read_integer_at_pos


    FUNCTION return_double_parameter( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_double_parameter


    FUNCTION return_position( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: res
      !# Centers of the matter objects. The first index runs over the matter
      !  objects, the second index over \((x,y,z)\).

    END FUNCTION return_position


    FUNCTION return_integer_parameter( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      INTEGER:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_integer_parameter


    SUBROUTINE return_eos_parameters_int( THIS, i_matter, eos_params )
    !# INTERFACE for a PROCEDURE that returns an array containing the
    !  parametersf the |eos| for the matter objects

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE return_eos_parameters_int


    FUNCTION return_string_parameter( THIS, i_matter ) RESULT( string )
    !# INTERFACE for a PROCEDURE that returns a CHARACTER( LEN= : )

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT ):: THIS
      !! [[idbase]] object which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER( LEN= : ), ALLOCATABLE:: string

    END FUNCTION return_string_parameter


    SUBROUTINE read_id_mass_b_int( THIS, x, y, z, &
                                   g_xx, &
                                   baryon_density, &
                                   gamma_euler )
      !# INTERFACE or the SUBROUTINE reading the hydro ID needed to compute the
      !  baryon mass

      IMPORT:: idbase
      CLASS(idbase),    INTENT( IN OUT ):: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN )    :: z
      DOUBLE PRECISION, INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( IN OUT ):: gamma_euler

    END SUBROUTINE read_id_mass_b_int


    SUBROUTINE read_id_full_int( THIS, n, x, y, z, &
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      g_xx, g_xy, g_xz, &
                                      g_yy, g_yz, g_zz, &
                                      k_xx, k_xy, k_xz, &
                                      k_yy, k_yz, k_zz, &
                                      baryon_density, &
                                      energy_density, &
                                      specific_energy, &
                                      u_euler_x, u_euler_y, u_euler_z )
     !# INTERFACE or the SUBROUTINE reading the full ID

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                  INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE read_id_full_int


    SUBROUTINE read_id_spacetime_int( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
     !# INTERFACE or the SUBROUTINE reading the spacetime ID

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                        INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE read_id_spacetime_int


    SUBROUTINE read_id_hydro_int( THIS, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# INTERFACE or the SUBROUTINE reading the the hydro ID needed to compute
    !  the constraints on the refined mesh

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                        INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE read_id_hydro_int


    SUBROUTINE read_id_particles_int( THIS, n, x, y, z, &
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      g_xx, g_xy, g_xz, &
                                      g_yy, g_yz, g_zz, &
                                      baryon_density, &
                                      energy_density, &
                                      specific_energy, &
                                      pressure, &
                                      u_euler_x, u_euler_y, u_euler_z )
    !# INTERFACE or the SUBROUTINE reading the hydro ID needed to compute the
    !  SPH ID

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                     INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE read_id_particles_int


    SUBROUTINE read_id_k_int( THIS, n, x, y, z,&
                              k_xx, k_xy, k_xz, &
                              k_yy, k_yz, k_zz )
    !# INTERFACE or the SUBROUTINE reading the components of the extrinsic
    !  curvature

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                     INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_zz

    END SUBROUTINE read_id_k_int


    SUBROUTINE integrate_field_int( THIS, center, radius, &
                                    central_density, &
                                    dr, dth, dphi, &
                                    mass, mass_profile, &
                                    mass_profile_idx )
    !# INTERFACE to the SUBROUTINE integrating the baryon mass density to
    !  compute the radial mass profile of a single star.

      IMPORT:: idbase
      !> Object of class [[idbase]] which this PROCEDURE is a member of
      CLASS(idbase), INTENT( IN OUT )      :: THIS
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
      !> Center of the star
      DOUBLE PRECISION, INTENT( IN )    :: center
      !> Central density of the star
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !> Radius of the star
      DOUBLE PRECISION, INTENT( IN )    :: radius
      !> Integration steps
      DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION, INTENT( IN OUT ):: mass
      !> Array storing the radial mass profile of the star
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                       mass_profile

    END SUBROUTINE integrate_field_int


    FUNCTION return_spatial_extent_int( THIS, i_matter ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  matter objects, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN OUT )   :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION return_spatial_extent_int


    FUNCTION return_total_spatial_extent_int( THIS ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  physical system considered, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )   :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION return_total_spatial_extent_int

  END INTERFACE


  INTERFACE

    MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
                                                     central_density, &
                                                     dr, dth, dphi, &
                                                     mass, mass_profile, &
                                                     mass_profile_idx )
    !# INTERFACE to the SUBROUTINE integrating the baryon mass density to
    !  compute the radial mass profile of a single star.

      !> Object of class [[idbase]] which this PROCEDURE is a member of
      CLASS(idbase), INTENT( IN OUT )      :: THIS
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
      !> Center of the star
      DOUBLE PRECISION, INTENT( IN )    :: center
      !> Central density of the star
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !> Radius of the star
      DOUBLE PRECISION, INTENT( IN )    :: radius
      !> Integration steps
      DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION, INTENT( IN OUT ):: mass
      !> Array storing the radial mass profile of the star
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                       mass_profile

    END SUBROUTINE integrate_baryon_mass_density


    MODULE FUNCTION get_n_matter( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      DOUBLE PRECISION:: get_n_matter
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_n_matter


    MODULE SUBROUTINE set_n_matter( THIS, value )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: value
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_n_matter


    MODULE SUBROUTINE check_i_matter( THIS, i_matter )
    !# Checks that the given index `i_matter` is between 1 and
    !  [[idbase:n_matter]], included. If not, it stops the execution of the
    !  program.

      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Value to be checked

    END SUBROUTINE check_i_matter


    MODULE FUNCTION get_total_spatial_extent( THIS ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  physical system considered, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      CLASS(idbase), INTENT( IN OUT )   :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION get_total_spatial_extent


    MODULE FUNCTION get_one_lapse( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      LOGICAL:: get_one_lapse
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_one_lapse


    MODULE SUBROUTINE set_one_lapse( THIS, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      LOGICAL, INTENT( IN ):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_one_lapse


    MODULE FUNCTION get_zero_shift( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      LOGICAL:: get_zero_shift
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_zero_shift


    MODULE SUBROUTINE set_zero_shift( THIS, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      LOGICAL, INTENT( IN ):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_zero_shift


  END INTERFACE


END MODULE id_base
