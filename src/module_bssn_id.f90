! File:         module_bssn_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE formul_bssn_id

  !***********************************************************
  !                                                          *
  !   This module contains the definition of TYPE bssn_id    *
  !                                                          *
  !***********************************************************


  USE utility,          ONLY: ios, err_msg, perc, creturn, run_id, &
                              test_status, compute_g4, &
                              determinant_sym4x4_grid, show_progress
  USE id_base,          ONLY: idbase
  USE formul_3p1_id,    ONLY: formul_3p1
  USE particles_id,     ONLY: particles
  USE timing,           ONLY: timer
  USE mesh_refinement,  ONLY: grid_function_scalar, grid_function


  IMPLICIT NONE


  !********************************************************
  !                                                       *
  !              Definition of TYPE bssn_id               *
  !                                                       *
  ! This class extends the ABSTRACT TYPE formul_3p1 by    *
  ! implementing its deferred methods such that the BSSN  *
  ! variable sare computed on the grid for the LORENE ID, *
  ! stored, exported to a binary file for evolution and   *
  ! to a formatted file. The BSSN constraints can also    *
  ! be computed in different ways, analyzed, and exported *
  ! in different ways.                                    *
  !                                                       *
  !********************************************************

  TYPE, EXTENDS( formul_3p1 ):: bssn_id


    INTEGER:: call_flag= 0
    ! Flag set to a value different than 0 if the SUBROUTINE
    ! compute_and_export_bssn_variables is called

    !
    !-- Arrays storing the BSSN variables for the LORENE ID on the grid
    !

    TYPE(grid_function):: Gamma_u
    !! Conformal connection

    TYPE(grid_function_scalar):: phi
    !! Conformal factor

    TYPE(grid_function_scalar):: trK
    !! Trace of extrinsic curvature

    TYPE(grid_function):: A_BSSN3_ll
    !! Conformal traceless extrinsic curvature

    TYPE(grid_function):: g_BSSN3_ll
    !! Conformal spatial metric

    !
    !-- Connection constraints and its l2 norm and loo norm
    !

    TYPE(grid_function):: GC
    TYPE(grid_function):: GC_parts

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_l2
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_parts_l2
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_loo
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_parts_loo

    LOGICAL, PUBLIC:: export_bin
    LOGICAL, PUBLIC:: export_form_xy, export_form_x

    TYPE(timer):: bssn_computer_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!


    PROCEDURE :: define_allocate_fields => allocate_bssn_fields

    PROCEDURE :: compute_and_export_3p1_variables &
                    => compute_and_export_bssn_variables

    PROCEDURE, PUBLIC :: read_bssn_dump_print_formatted

    PROCEDURE :: print_formatted_lorene_id_3p1_variables &
                    => print_formatted_lorene_id_bssn_variables

    PROCEDURE :: compute_and_export_3p1_constraints_grid &
                    => compute_and_export_bssn_constraints_grid

    PROCEDURE :: compute_and_export_3p1_constraints_particles &
                    => compute_and_export_bssn_constraints_particles

    PROCEDURE :: deallocate_fields => deallocate_bssn_fields


    PROCEDURE :: destruct_bssn_id
    !# Finalizer for members of the extended class bssn_id, not the
    !  primitive class formul_3p1

    FINAL     :: destructor
    !# Destructor; finalizes members from both CLASSES formul_3p1, and bssn_id,
    !  by calling destruct_formul_3p1 and destruct_bssn_id

  END TYPE bssn_id

  !
  !-- Interface of the TYPE bssn_id
  !-- (i.e., declaration of the overloaded constructor)
  !
  INTERFACE bssn_id

    MODULE PROCEDURE:: construct_bssn_id
    !# Constructs the bssn_id object from the number of grid points
    !  along each axis

  END INTERFACE bssn_id

  !
  !-- Interface of the constructor of TYPE bssn_id
  !-- Its implementation is in submodule_BSSN_id_constructor.f90
  !
  INTERFACE

    MODULE FUNCTION construct_bssn_id( id, dx, dy, dz ) RESULT ( bssnid )
    !# Constructs the bssn_id object from the number of grid points
    !  along each axis

      CLASS(idbase), INTENT( IN OUT ):: id
      TYPE(bssn_id)                  :: bssnid
      DOUBLE PRECISION, OPTIONAL     :: dx, dy, dz

    END FUNCTION construct_bssn_id

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE bssn_id
  !-- Their implementations are in submodule_BSSN_id_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE allocate_bssn_fields( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE allocate_bssn_fields

    MODULE SUBROUTINE compute_and_export_bssn_variables( THIS, namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_bssn_variables

    MODULE SUBROUTINE read_bssn_dump_print_formatted( THIS, namefile_bin, &
                                                            namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile_bin
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_bssn_dump_print_formatted

    MODULE SUBROUTINE print_formatted_lorene_id_bssn_variables( THIS, &
                                                                namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_bssn_variables

    MODULE SUBROUTINE compute_and_export_bssn_constraints_grid( THIS, &
                                                           id, &
                                                           namefile, &
                                                           name_logfile )

      CLASS(bssn_id),      INTENT( IN OUT ):: THIS
      CLASS(idbase),      INTENT( IN OUT ):: id
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_grid

    MODULE SUBROUTINE compute_and_export_bssn_constraints_particles( THIS, &
                                                           parts_obj, &
                                                           namefile, &
                                                           name_logfile )

      CLASS(bssn_id),      INTENT( IN OUT ):: THIS
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_particles

    MODULE SUBROUTINE deallocate_bssn_fields( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_bssn_fields

    MODULE SUBROUTINE destruct_bssn_id( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_bssn_id

    MODULE SUBROUTINE destructor( THIS )

      TYPE(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE destructor

  END INTERFACE


END MODULE formul_bssn_id
