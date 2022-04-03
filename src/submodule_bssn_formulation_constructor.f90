! File:         submodule_bssn_formulation_constructor.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020, 2021, 2022 Francesco Torsello                     *
!                                                                       *
! This file is part of SPHINCS_ID                                       *
!                                                                       *
! SPHINCS_ID is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by  *
! the Free Software Foundation, either version 3 of the License, or     *
! (at your option) any later version.                                   *
!                                                                       *
! SPHINCS_ID is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of        *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
! GNU General Public License for more details.                          *
!                                                                       *
! You should have received a copy of the GNU General Public License     *
! along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/>.   *
! The copy of the GNU General Public License should be in the file      *
! 'COPYING'.                                                            *
!************************************************************************

SUBMODULE (bssn_formulation) constructor

  !************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE bssn
  !
  !  FT 23.10.2020
  !
  !  Updated to support mesh refinement
  !
  !  FT 26.03.2021
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE construct_bssn

    !****************************************************
    !
    !# This constructor of TYPE bssn calls the
    !  SUBROUTINES that rely on an idbase object, and
    !  allocates memory. It constructs the grid
    !  using the number of grid points along each axis.
    !
    !  FT 23.10.2020
    !
    !****************************************************

    USE McLachlan_refine, ONLY: initialize_BSSN, deallocate_BSSN
    USE mesh_refinement,  ONLY: levels, allocate_grid_function
    USE Extract_Mass,     ONLY: radius2
    USE constants,        ONLY: zero, one

    IMPLICIT NONE

    ! Initialize the timer
    bssnid% bssn_computer_timer= timer( "bssn_computer_timer" )

    ! Construct the gravity grid and import the LORENE ID on it,
    ! in standard 3+1 formulation
    IF( PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) )THEN

      CALL bssnid% setup_standard_tpo_variables( id, dx, dy, dz )

    ELSE

      CALL bssnid% setup_standard_tpo_variables( id )

    ENDIF

    ! Read and store the BSSN parameters
    CALL initialize_BSSN()
    CALL deallocate_BSSN()

    ! The construct_formul_3p1 SUBROUTINE constructs the grid,
    ! hence the dimensions of the arrays imported from the module BSSN
    ! are know and the arrays can be allocated
    !CALL allocate_bssn_fields( bssnid )

    DEALLOCATE( levels )

    IF( .NOT.ALLOCATED( bssnid% GC_int ))THEN
      ALLOCATE( bssnid% GC_int( bssnid% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array GC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    bssnid% GC_int= HUGE(one)

    IF( .NOT.ALLOCATED( bssnid% GC_parts_int ))THEN
      ALLOCATE( bssnid% GC_parts_int( bssnid% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array GC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    bssnid% GC_parts_int= HUGE(one)

    ! radius2 is the extraction radius. If not set here, then it is 0 by default
    ! and the metric is not interpolated on the particle in
    ! get_metric_on_particles
    radius2= HUGE(DBLE(1.0D0))

    PRINT *
    PRINT *, " * Ready to compute BSSN variables."
    PRINT *

  END PROCEDURE construct_bssn

  !
  !-- Keeping the following two SUBROUTINES separate in case it is needed
  !-- to add other PROCEDURES to the destructor (probably superfluous...)
  !
  MODULE PROCEDURE destruct_bssn

    !**************************************************
    !
    !# Destructor of the EXTENDED TYPE bssn
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    CALL deallocate_bssn_fields( THIS )

  END PROCEDURE destruct_bssn


  MODULE PROCEDURE destructor

    !**************************************************
    !
    !# Destructor of EXTENDED TYPE bssn
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    CALL destruct_bssn( THIS )

#if host == Sunrise

  CALL THIS% deallocate_standard_tpo_variables

#else

#ifdef __INTEL_COMPILER
  CALL deallocate_standard_tpo_variables( THIS )
#endif

#ifdef __GFORTRAN__
  CALL THIS% deallocate_standard_tpo_variables
#endif

#endif

  END PROCEDURE destructor


END SUBMODULE constructor
