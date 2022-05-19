! File:         submodule_bssn_formulation_constraints.f90
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

SUBMODULE (bssn_formulation) constraints

  !************************************************
  !
  !# Implementation of the methods of TYPE bssn
  !  that compute the constraints
  !
  !  FT 9.07.2021
  !
  !************************************************

  USE utility,  ONLY: zero, one, two, three, four, five, ten

  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_export_bssn_constraints_grid

    !***************************************************
    !
    !# Compute, store, analyze and print the BSSN
    !  constraints to a formatted file. The computation
    !  is done by importing the hydro ID on the
    !  gravity grid, without any information on the
    !  particles.
    !
    !  FT 1.02.2021
    !
    !***************************************************

    USE constants,         ONLY: pi
    USE utility,           ONLY: lorene2hydrobase
    USE matrix,            ONLY: invert_4x4_matrix
    USE tensor,            ONLY: itt, itx, ity, itz, ixx, ixy, &
                                 ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                 jyy, jyz, jzz, jx, jy, jz, &
                                 it, ix, iy, iz, n_sym4x4
    USE mesh_refinement,   ONLY: allocate_grid_function, &
                                 levels, nlevels
    USE McLachlan_refine,  ONLY: BSSN_CONSTRAINTS_INTERIOR

    IMPLICIT NONE

    INTEGER:: i, j, k, fd_lim, l
    INTEGER, DIMENSION(3) :: imin, imax
    INTEGER:: unit_logfile, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    TYPE(grid_function_scalar):: baryon_density
    TYPE(grid_function_scalar):: energy_density
    TYPE(grid_function_scalar):: specific_energy
    TYPE(grid_function_scalar):: pressure
    TYPE(grid_function):: v_euler
    TYPE(grid_function):: v_euler_l
    TYPE(grid_function):: u_euler_l
    TYPE(grid_function_scalar):: lorentz_factor
    DOUBLE PRECISION:: u_euler_norm= zero
    DOUBLE PRECISION:: detg4
    ! Spacetime metric
    TYPE(grid_function):: g4
    ! Stress-energy tensor
    TYPE(grid_function):: Tmunu_ll
    ! Spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: g4temp
    ! Inverse spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: ig4

    ! Declaration of debug variables needed to compute the Hamiltonian
    ! constraint directly, without calling the Cactus-bound SUBROUTINE
    ! BSSN_CONSTRAINTS_INTERIOR
    TYPE(grid_function_scalar):: HC_hand
    TYPE(grid_function_scalar):: HC_rho
    TYPE(grid_function_scalar):: HC_trK
    TYPE(grid_function_scalar):: HC_A
    TYPE(grid_function_scalar):: HC_derphi

    CHARACTER( LEN= : ), ALLOCATABLE:: name_constraint
    CHARACTER( LEN= : ), ALLOCATABLE:: name_analysis
    CHARACTER( LEN= : ), ALLOCATABLE:: finalname_logfile
    CHARACTER( LEN= 2 ):: n_reflev

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    ALLOCATE ( levels( this% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    levels = this% levels
    nlevels= this% nlevels

    CALL allocate_grid_function( baryon_density, "baryon_density", 1 )
    CALL allocate_grid_function( energy_density, "energy_density", 1 )
    CALL allocate_grid_function( specific_energy, "specific_energy", 1 )
    CALL allocate_grid_function( pressure, "pressure", 1 )

    CALL allocate_grid_function( v_euler, "v_euler", 3 )
    CALL allocate_grid_function( v_euler_l, "v_euler_l", 3 )
    CALL allocate_grid_function( u_euler_l, "u_euler_l", 4 )
    CALL allocate_grid_function( lorentz_factor, "lorentz_factor", 1 )

    CALL allocate_grid_function( g4, "g4", n_sym4x4 )
    CALL allocate_grid_function( Tmunu_ll, "Tmunu_ll", n_sym4x4 )

    CALL allocate_grid_function( HC_hand, "HC_hand", 1 )
    CALL allocate_grid_function( HC_rho, "HC_rho", 1 )
    CALL allocate_grid_function( HC_trK, "HC_trK", 1 )
    CALL allocate_grid_function( HC_A, "HC_A", 1 )
    CALL allocate_grid_function( HC_derphi, "HC_derphi", 1 )

    CALL allocate_grid_function( this% HC, "HC_id", 1 )
    CALL allocate_grid_function( this% MC, "MC_id", 3 )
    CALL allocate_grid_function( this% GC, "GC_id", 3 )

    CALL allocate_grid_function( this% rho, "HC_rho", 1 )
    CALL allocate_grid_function( this% S, "MC_S", 3 )

    !
    !-- Import the hydro ID on the gravity grid
    !
    PRINT *, "** Importing the hydro ID on the mesh..."
    PRINT *
    ref_levels: DO l= 1, this% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      CALL id% read_id_hydro( this% get_ngrid_x(l), &
                              this% get_ngrid_y(l), &
                              this% get_ngrid_z(l), &
                              this% coords% levels(l)% var, &
                              baryon_density% levels(l)% var, &
                              energy_density% levels(l)% var, &
                              specific_energy% levels(l)% var, &
                              pressure% levels(l)% var, &
                              v_euler% levels(l)% var )

    ENDDO ref_levels
    PRINT *, " * Hydro ID imported."
    PRINT *

    !---------------------------!
    !--  Compute constraints  --!
    !---------------------------!

    !
    !-- Compute the fluid 4-velocity in the coordinate frame
    !

    ref_levels2: DO l= 1, this% nlevels

      PRINT *, "** Computing fluid 4-velocity wrt Eulerian observer ", &
               "on refinement level ", l, "..."

      CALL compute_4velocity_eul()
      PRINT *, " * Fluid 4-velocity wrt Eulerian observer ", &
               "on refinement level", l, "computed."
      PRINT *

      ! Note that the units used in the spacetime part of SPHINCS are the
      ! same units as in the HydroBase thorn in the Einstein Toolkit.
      ! Such units can be found here, https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html
      ! The order of magnitude of the energy density can be found in
      ! https://www.ias.ac.in/article/fulltext/pram/084/05/0927-0941,
      ! and it is 150 MeV fm^{-3} ~ (2.4*10^{-11}J) / (10^{-45}m^3)
      !                           = 2.4*10^34 J m^{-3}

      !
      !-- Compute the stress-energy tensor
      !
      PRINT *, "** Computing stress-energy tensor ", &
               "on refinement level ", l, "..."

      CALL compute_stress_energy()
      PRINT *, " * Stress-energy tensor on refinement level", l, "computed."
      PRINT *

    ENDDO ref_levels2

    ! In debug mode, compute the Hamiltonian constraint by hand
    IF( debug )THEN

      DO l= 1, this% nlevels, 1

#ifdef __INTEL_COMPILER

  ASSOCIATE( HC_rho         => HC_rho% levels(l)%var, &
             HC_trK         => HC_trK% levels(l)%var, &
             HC_A           => HC_A% levels(l)%var, &
             HC_derphi      => HC_derphi% levels(l)%var, &
             HC_hand        => HC_hand% levels(l)%var, &
             phi            => this% phi% levels(l)%var, &
             trK            => this% trK% levels(l)%var, &
             A_BSSN3_ll     => this% A_BSSN3_ll% levels(l)%var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

  HC_rho= zero
  HC_trK= zero
  HC_A= zero
  HC_derphi= zero
  HC_hand= zero

#endif

#ifdef __GFORTRAN__

  HC_rho% levels(l)%var= zero
  HC_trK% levels(l)%var= zero
  HC_A% levels(l)%var= zero
  HC_derphi% levels(l)%var= zero
  HC_hand% levels(l)%var= zero

#endif
          fd_lim= 5

          DO k= fd_lim, this% get_ngrid_z(l) - fd_lim, 1
            DO j= fd_lim, this% get_ngrid_y(l) - fd_lim, 1
              DO i= fd_lim, this% get_ngrid_x(l) - fd_lim, 1

#ifdef __GFORTRAN__

  ASSOCIATE( HC_rho         => HC_rho% levels(l)%var, &
             HC_trK         => HC_trK% levels(l)%var, &
             HC_A           => HC_A% levels(l)%var, &
             HC_derphi      => HC_derphi% levels(l)%var, &
             HC_hand        => HC_hand% levels(l)%var, &
             phi            => this% phi% levels(l)%var, &
             trK            => this% trK% levels(l)%var, &
             A_BSSN3_ll     => this% A_BSSN3_ll% levels(l)%var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

#endif

! The following works with both compilers
!                ASSOCIATE( HC_rho => HC_rho% levels(l)%var( i, j, k ) &
!                )

                HC_rho( i, j, k )= two*pi*EXP(five*phi( i, j, k )) &
                                    *lorene2hydrobase*energy_density( i, j, k )

                HC_trK( i, j, k )= - EXP(five*phi( i, j, k ))/(three*four) &
                                        *trK( i, j, k )**2

                HC_A( i, j, k )= EXP(five*phi( i, j, k ))/two*four &
                 *( A_BSSN3_ll(i, j, k,jxx)*A_BSSN3_ll(i, j, k,jxx) &
                  + A_BSSN3_ll(i, j, k,jxy)*A_BSSN3_ll(i, j, k,jxy) &
                  + A_BSSN3_ll(i, j, k,jxz)*A_BSSN3_ll(i, j, k,jxz) &
                  + A_BSSN3_ll(i, j, k,jxy)*A_BSSN3_ll(i, j, k,jxy) &
                  + A_BSSN3_ll(i, j, k,jyy)*A_BSSN3_ll(i, j, k,jyy) &
                  + A_BSSN3_ll(i, j, k,jyz)*A_BSSN3_ll(i, j, k,jyz) &
                  + A_BSSN3_ll(i, j, k,jxz)*A_BSSN3_ll(i, j, k,jxz) &
                  + A_BSSN3_ll(i, j, k,jyz)*A_BSSN3_ll(i, j, k,jyz) &
                  + A_BSSN3_ll(i, j, k,jzz)*A_BSSN3_ll(i, j, k,jzz) &
                  )

                ! Second derivative of conformal factor with fourth-order FD
                !HC_derphi( ix, iy, iz )= &
                !                   ( -      EXP(this% phi( ix + 2, iy, iz )) &
                !                     + 16.0*EXP(this% phi( ix + 1, iy, iz )) &
                !                     - 30.0*EXP(this% phi( ix    , iy, iz )) &
                !                     + 16.0*EXP(this% phi( ix - 1, iy, iz )) &
                !                     -      EXP(this% phi( ix - 2, iy, iz )) &
                !                     -      EXP(this% phi( ix, iy + 2, iz )) &
                !                     + 16.0*EXP(this% phi( ix, iy + 1, iz )) &
                !                     - 30.0*EXP(this% phi( ix, iy, iz )) &
                !                     + 16.0*EXP(this% phi( ix, iy - 1, iz )) &
                !                     -      EXP(this% phi( ix, iy - 2, iz )) &
                !                     -      EXP(this% phi( ix, iy, iz + 2 )) &
                !                     + 16.0*EXP(this% phi( ix, iy, iz + 1 )) &
                !                     - 30.0*EXP(this% phi( ix, iy, iz )) &
                !                     + 16.0*EXP(this% phi( ix, iy, iz - 1 )) &
                !                     -      EXP(this% phi( ix, iy, iz - 2 )) )&
                !                     /(12.0*this% dx**2)

                ! Second derivative of conformal factor with eighth-order FD
                HC_derphi( i, j, k )= ( &
                                - DBLE(one/560.0)*EXP(phi(i + 4, j, k)) &
                                + DBLE(two*four/315.0)*EXP(phi(i + 3, j, k)) &
                                - DBLE(one/five  )*EXP(phi(i + 2, j, k)) &
                                + DBLE(two*four/five  )*EXP(phi(i + 1, j, k)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(two*four/five  )*EXP(phi(i - 1, j, k)) &
                                - DBLE(one/five  )*EXP(phi(i - 2, j, k)) &
                                + DBLE(two*four/315.0)*EXP(phi(i - 3, j, k)) &
                                - DBLE(one/560.0)*EXP(phi(i - 4, j, k)) &
                                - DBLE(one/560.0)*EXP(phi(i, j + 4, k)) &
                                + DBLE(two*four/315.0)*EXP(phi(i, j + 3, k)) &
                                - DBLE(one/five  )*EXP(phi(i, j + 2, k)) &
                                + DBLE(two*four/five  )*EXP(phi(i, j + 1, k)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(two*four/five  )*EXP(phi(i, j - 1, k)) &
                                - DBLE(one/five  )*EXP(phi(i, j - 2, k)) &
                                + DBLE(two*four/315.0)*EXP(phi(i, j - 3, k)) &
                                - DBLE(one/560.0)*EXP(phi(i, j - 4, k)) &
                                - DBLE(one/560.0)*EXP(phi(i, j, k + 4)) &
                                + DBLE(two*four/315.0)*EXP(phi(i, j, k + 3)) &
                                - DBLE(one/five  )*EXP(phi(i, j, k + 2)) &
                                + DBLE(two*four/five  )*EXP(phi(i, j, k + 1)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(two*four/five  )*EXP(phi(i, j, k - 1)) &
                                - DBLE(one/five  )*EXP(phi(i, j, k - 2)) &
                                + DBLE(two*four/315.0)*EXP(phi(i, j, k - 3)) &
                                - DBLE(one/560.0)*EXP(phi(i, j, k - 4)) )&
                                /(this% levels(l)% dx**2)


                HC_hand( i, j, k )= HC_rho( i, j, k ) + &
                                    HC_trK( i, j, k ) + &
                                    HC_A( i, j, k )   + &
                                    HC_derphi( i, j, k )

#ifdef __GFORTRAN__
  END ASSOCIATE
#endif

              ENDDO
            ENDDO
          ENDDO

#ifdef __INTEL_COMPILER
  END ASSOCIATE
#endif

      ENDDO

    ENDIF

    !
    !-- Compute the BSSN constraints by calling the Cactus-bound procedure
    !-- BSSN_CONSTRAINTS_INTERIOR
    !
    PRINT *, "** Computing contraints..."
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP          SHARED( this, Tmunu_ll ) &
!    !$OMP          PRIVATE( l, imin, imax )
    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse      => this% lapse% levels(l)% var, &
                 shift_u    => this% shift_u% levels(l)% var, &
                 phi        => this% phi% levels(l)% var, &
                 trK        => this% trK% levels(l)% var, &
                 g_BSSN3_ll => this% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => this% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u    => this% Gamma_u% levels(l)% var, &
                 Tmunu_ll   => Tmunu_ll% levels(l)% var, &
                 HC         => this% HC% levels(l)% var, &
                 MC         => this% MC% levels(l)% var, &
                 GC         => this% GC% levels(l)% var, &
                 rho        => this% rho% levels(l)% var, &
                 S          => this% S% levels(l)% var &
      )

        imin(1) = this% levels(l)% nghost_x
        imin(2) = this% levels(l)% nghost_y
        imin(3) = this% levels(l)% nghost_z
        imax(1) = this% get_ngrid_x(l) - this% levels(l)% nghost_x - 1
        imax(2) = this% get_ngrid_y(l) - this% levels(l)% nghost_y - 1
        imax(3) = this% get_ngrid_z(l) - this% levels(l)% nghost_z - 1

        HC = zero
        MC = zero
        GC = zero
        rho= zero
        S  = zero
        CALL bssn_constraint_terms_interior( &
          !
          !-- Input
          !
          this% get_ngrid_x(l), this% get_ngrid_y(l), this% get_ngrid_z(l), &
          imin, imax, &
          this% get_dx(l), this% get_dy(l), this% get_dz(l), &
          g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
          g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
          g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
          A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
          A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
          A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
          trK(:,:,:), phi(:,:,:), &
          Gamma_u(:,:,:,jx), &
          Gamma_u(:,:,:,jy), &
          Gamma_u(:,:,:,jz), &
          Tmunu_ll(:,:,:,itt), &
          Tmunu_ll(:,:,:,itx), &
          Tmunu_ll(:,:,:,ity), &
          Tmunu_ll(:,:,:,itz), &
          Tmunu_ll(:,:,:,ixx), &
          Tmunu_ll(:,:,:,ixy), &
          Tmunu_ll(:,:,:,ixz), &
          Tmunu_ll(:,:,:,iyy), &
          Tmunu_ll(:,:,:,iyz), &
          Tmunu_ll(:,:,:,izz), &
          lapse(:,:,:), &
          shift_u(:,:,:,jx), &
          shift_u(:,:,:,jy), &
          shift_u(:,:,:,jz), &
          !
          !-- Output
          !
          ! Connection constraints
          GC(:,:,:,jx), &
          GC(:,:,:,jy), &
          GC(:,:,:,jz), &
          ! Hamiltonian and momentum constraints
          HC(:,:,:), &
          MC(:,:,:,jx), &
          MC(:,:,:,jy), &
          MC(:,:,:,jz), &
          ! Sources in the Hamiltonian and momentum constraints
          rho(:,:,:), &
          S(:,:,:,jx), &
          S(:,:,:,jy), &
          S(:,:,:,jz) &
        )
       ! CALL BSSN_CONSTRAINTS_INTERIOR( &
       !   !
       !   !-- Input
       !   !
       !   this% get_ngrid_x(l), this% get_ngrid_y(l), this% get_ngrid_z(l), &
       !   imin, imax, &
       !   this% get_dx(l), this% get_dy(l), this% get_dz(l), &
       !   g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
       !   g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
       !   g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
       !   A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
       !   A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
       !   A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
       !   trK(:,:,:), phi(:,:,:), &
       !   Gamma_u(:,:,:,jx), &
       !   Gamma_u(:,:,:,jy), &
       !   Gamma_u(:,:,:,jz), &
       !   Tmunu_ll(:,:,:,itt), &
       !   Tmunu_ll(:,:,:,itx), &
       !   Tmunu_ll(:,:,:,ity), &
       !   Tmunu_ll(:,:,:,itz), &
       !   Tmunu_ll(:,:,:,ixx), &
       !   Tmunu_ll(:,:,:,ixy), &
       !   Tmunu_ll(:,:,:,ixz), &
       !   Tmunu_ll(:,:,:,iyy), &
       !   Tmunu_ll(:,:,:,iyz), &
       !   Tmunu_ll(:,:,:,izz), &
       !   lapse(:,:,:), &
       !   shift_u(:,:,:,jx), &
       !   shift_u(:,:,:,jy), &
       !   shift_u(:,:,:,jz), &
       !   !
       !   !-- Output
       !   !
       !   ! Connection constraints
       !   GC(:,:,:,jx), &
       !   GC(:,:,:,jy), &
       !   GC(:,:,:,jz), &
       !   ! Hamiltonian and momentum constraints
       !   HC(:,:,:), &
       !   MC(:,:,:,jx), &
       !   MC(:,:,:,jy), &
       !   MC(:,:,:,jz) &
       ! )
      END ASSOCIATE
    ENDDO
!    !$OMP END PARALLEL DO
    PRINT *, " * Constraints computed."
    PRINT *

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    DO l= 1, this% nlevels, 1

      ASSOCIATE( HC         => this% HC% levels(l)% var, &
                 MC         => this% MC% levels(l)% var, &
                 GC         => this% GC% levels(l)% var, &
                 rho        => this% rho% levels(l)% var, &
                 S          => this% S% levels(l)% var &
      )
        !
        !-- Export the constraint statistics to a formatted file
        !
        unit_logfile= 2891

        IF( l > 9 )THEN
          WRITE( n_reflev, "(I2)" ) l
        ELSE
          WRITE( n_reflev, "(I1)" ) l
        ENDIF

        finalname_logfile= TRIM(name_logfile)//"-reflev"//TRIM(n_reflev)//".log"

        INQUIRE( FILE= TRIM(finalname_logfile), EXIST= exist )

        IF( exist )THEN
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "REPLACE", &
                  FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "NEW", &
                  FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening ", TRIM(finalname_logfile), &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error when opening " &
        !         // TRIM(name_logfile) )

        IF( .NOT.ALLOCATED( this% HC_l2 ))THEN
          ALLOCATE( this% HC_l2( this% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% MC_l2 ))THEN
          ALLOCATE( this% MC_l2( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% GC_l2 ))THEN
          ALLOCATE( this% GC_l2( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% HC_loo ))THEN
          ALLOCATE( this% HC_loo( this% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% MC_loo ))THEN
          ALLOCATE( this% MC_loo( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% GC_loo ))THEN
          ALLOCATE( this% GC_loo( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
      !  IF( .NOT.ALLOCATED( this% HC_int ))THEN
      !    ALLOCATE( this% HC_int( this% nlevels ), &
      !              STAT= ios, ERRMSG= err_msg )
      !    IF( ios > 0 )THEN
      !      PRINT *, "...allocation error for array HC_loo. ", &
      !               "The error message is", err_msg
      !      STOP
      !    ENDIF
      !    !CALL test_status( ios, err_msg, &
      !    !                "...deallocation error for array HC" )
      !  ENDIF
      !  IF( .NOT.ALLOCATED( this% MC_int ))THEN
      !    ALLOCATE( this% MC_int( this% nlevels, 3 ), &
      !              STAT= ios, ERRMSG= err_msg )
      !    IF( ios > 0 )THEN
      !      PRINT *, "...allocation error for array MC_loo. ", &
      !               "The error message is", err_msg
      !      STOP
      !    ENDIF
      !    !CALL test_status( ios, err_msg, &
      !    !                "...deallocation error for array HC" )
      !  ENDIF
      !  IF( .NOT.ALLOCATED( this% GC_int ))THEN
      !    ALLOCATE( this% GC_int( this% nlevels, 3 ), &
      !              STAT= ios, ERRMSG= err_msg )
      !    IF( ios > 0 )THEN
      !      PRINT *, "...allocation error for array GC_loo. ", &
      !               "The error message is", err_msg
      !      STOP
      !    ENDIF
      !    !CALL test_status( ios, err_msg, &
      !    !                "...deallocation error for array HC" )
      !  ENDIF

        WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

        PRINT *, "** Analyzing constraints on refinement level ", l, "..."

        name_analysis= "bssn-hc-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the Hamiltonian constraint"
        CALL this% analyze_constraint( &
             l, &
             HC, name_constraint, unit_logfile, name_analysis, &
             this% HC_l2(l), this% HC_loo(l), this% HC_int(l), rho )

        name_analysis= "bssn-mc1-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             this% MC_l2(l,jx), this% MC_loo(l,jx), this% MC_int(l,jx), &
             S(:,:,:,jx) )

        name_analysis= "bssn-mc2-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             this% MC_l2(l,jy), this% MC_loo(l,jy), this% MC_int(l,jy), &
             S(:,:,:,jy) )

        name_analysis= "bssn-mc3-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             this% MC_l2(l,jz), this% MC_loo(l,jz), this% MC_int(l,jz), &
             S(:,:,:,jz) )

        name_analysis= "bssn-gc1-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             this% GC_l2(l,jx), this% GC_loo(l,jx), this% GC_int(l,jx) )

        name_analysis= "bssn-gc2-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             this% GC_l2(l,jy), this% GC_loo(l,jy), this% GC_int(l,jy) )

        name_analysis= "bssn-gc3-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             this% GC_l2(l,jz), this% GC_loo(l,jz), this% GC_int(l,jz) )

        CLOSE( UNIT= unit_logfile )

      PRINT *, " * Constraints analyzed. Summary of results saved to ", &
               finalname_logfile
      PRINT *
      END ASSOCIATE
    ENDDO

    IF( this% export_constraints )THEN

      PRINT *, "** Printing constraints to file ", TRIM(namefile), "..."

      !
      !-- Export the constraints to a formatted file
      !
      INQUIRE( FILE= TRIM(namefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
      ELSE
        OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "NEW", &
        FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(namefile) )

      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the stress-energy tensor and the BSSN constraints" &
      // " for the ID " &
      // "on selected grid points"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
      !         // TRIM(namefile) )
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4       5", &
      "       6       7       8       9       10", &
      "       11       12       13       14       15", &
      "       16       17       18       19       20      21"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 2 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
      !        // TRIM(namefile) )
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      refinement level    x   y   z   Stress-energy (10 components)   "&
      // "Hamiltonian constraint       " &
      // "Momentum constraint (three components)       " &
      // "Connection constraint (three components)"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 3 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
      !        // TRIM(namefile) )

      DO l= 1, this% nlevels, 1

        ASSOCIATE( lapse           => this% lapse% levels(l)% var, &
                   shift_u         => this% shift_u% levels(l)% var, &
                   phi             => this% phi% levels(l)% var, &
                   trK             => this% trK% levels(l)% var, &
                   g_BSSN3_ll      => this% g_BSSN3_ll% levels(l)% var, &
                   A_BSSN3_ll      => this% A_BSSN3_ll% levels(l)% var, &
                   g_phys3_ll      => this% g_phys3_ll% levels(l)% var, &
                   k_phys3_ll      => this% k_phys3_ll% levels(l)% var, &
                   Gamma_u         => this% Gamma_u% levels(l)% var, &
                   Tmunu_ll        => Tmunu_ll% levels(l)% var, &
                   v_euler_l       => v_euler_l% levels(l)% var, &
                   u_euler_l       => u_euler_l% levels(l)% var, &
                   v_euler         => v_euler% levels(l)% var, &
                   lorentz_factor  => lorentz_factor% levels(l)% var, &
                   HC              => this% HC% levels(l)% var, &
                   MC              => this% MC% levels(l)% var, &
                   GC              => this% GC% levels(l)% var, &
                   HC_rho          => HC_rho% levels(l)%var, &
                   HC_trK          => HC_trK% levels(l)%var, &
                   HC_A            => HC_A% levels(l)%var, &
                   HC_derphi       => HC_derphi% levels(l)%var, &
                   HC_hand         => HC_hand% levels(l)%var, &
                   g4              => g4% levels(l)% var, &
                   baryon_density  => baryon_density% levels(l)% var, &
                   specific_energy => specific_energy% levels(l)% var, &
                   energy_density  => energy_density% levels(l)% var, &
                   pressure        => pressure% levels(l)% var &
        )

          ! Being abs_grid a local array, it is good practice to allocate it on
          ! the heap, otherwise it will be stored on the stack which has a very
          ! limited size. This results in a segmentation fault.
          IF( ALLOCATED( abs_grid ) )THEN
            DEALLOCATE( abs_grid )
          ENDIF
          ALLOCATE( abs_grid( this% get_ngrid_x(l), this% get_ngrid_y(l), &
                              this% get_ngrid_z(l), 3 ) )

          DO k= 1, this% get_ngrid_z(l), 1
            DO j= 1, this% get_ngrid_y(l), 1
              DO i= 1, this% get_ngrid_x(l), 1

                abs_grid( i, j, k, jx )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jx ) )
                abs_grid( i, j, k, jy )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                abs_grid( i, j, k, jz )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jz ) )

              ENDDO
            ENDDO
          ENDDO

          min_abs_y= 1D+20
          min_abs_z= 1D+20
          DO k= 1, this% get_ngrid_z(l), 1
            DO j= 1, this% get_ngrid_y(l), 1
              DO i= 1, this% get_ngrid_x(l), 1

                IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                    < min_abs_y )THEN
                  min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                  min_ix_y= i
                  min_iy_y= j
                  min_iz_y= k
                ENDIF

                IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                    < min_abs_z )THEN
                  min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                  min_ix_z= i
                  min_iy_z= j
                  min_iz_z= k
                ENDIF

              ENDDO
            ENDDO
          ENDDO

          DO k= 1, this% get_ngrid_z(l), 1

            IF( MOD( k, this% cons_step ) /= 0 ) CYCLE

            DO j= 1, this% get_ngrid_y(l), 1

              IF( MOD( j, this% cons_step ) /= 0 ) CYCLE

              DO i= 1, this% get_ngrid_x(l), 1

                IF( MOD( i, this% cons_step ) /= 0 ) CYCLE

                IF( this% export_constraints_xy .AND. &
                    ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                      this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) ) )THEN
                  CYCLE
                ENDIF
                IF( this% export_constraints_x .AND. &
                    ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                      this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) &
                      .OR. &
                      this% coords% levels(l)% var( i, j, k, jy ) /= &
                      this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                    min_iz_y, jy ) ) )THEN
                  CYCLE
                ENDIF

                IF( debug )THEN
                  WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
                           FMT = * )&
                    l, &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    baryon_density( i, j, k ), &
                    energy_density( i, j, k ), &
                    specific_energy( i, j, k ), &
                    pressure( ix, iy, iz ), &
                    u_euler_l( i, j, k, it ), &
                    u_euler_l( i, j, k, ix ), &
                    u_euler_l( i, j, k, iy ), &
                    u_euler_l( i, j, k, iz ), &
                    u_euler_l( i, j, k, it ), &
                    u_euler_l( i, j, k, ix ), &
                    u_euler_l( i, j, k, iy ), &
                    u_euler_l( i, j, k, iz ), &
                    v_euler( i, j, k, jx ), &
                    v_euler( i, j, k, jy ), &
                    v_euler( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC( i, j, k ), &
                    HC_hand( i, j, k ), &
                    HC_rho( i, j, k ), &
                    HC_trK( i, j, k ), &
                    HC_A( i, j, k ), &
                    HC_derphi( i, j, k ), &
                    lorentz_factor( i, j, k ), &
                    lapse( ix, iy, iz ), &
                    shift_u( i, j, k, jx ), &
                    shift_u( i, j, k, jy ), &
                    shift_u( i, j, k, jz ), &
                    g4( i, j, k, ixx ), &
                    g4( i, j, k, ixy ), &
                    g4( i, j, k, ixz ), &
                    g4( i, j, k, iyy ), &
                    g4( i, j, k, iyz ), &
                    g4( i, j, k, izz ), &
                    !g_BSSN3_ll( i, j, k, jxx ), &
                    !g_BSSN3_ll( i, j, k, jxy ), &
                    !g_BSSN3_ll( i, j, k, jxz ), &
                    !g_BSSN3_ll( i, j, k, jyy ), &
                    !g_BSSN3_ll( i, j, k, jyz ), &
                    !g_BSSN3_ll( i, j, k, jzz ), &
                    k_phys3_ll( i, j, k, jxx ), &
                    k_phys3_ll( i, j, k, jxy ), &
                    k_phys3_ll( i, j, k, jxz ), &
                    k_phys3_ll( i, j, k, jyy ), &
                    k_phys3_ll( i, j, k, jyz ), &
                    k_phys3_ll( i, j, k, jzz ), &
                    A_BSSN3_ll( i, j, k, jxx ), &
                    A_BSSN3_ll( i, j, k, jxy ), &
                    A_BSSN3_ll( i, j, k, jxz ), &
                    A_BSSN3_ll( i, j, k, jyy ), &
                    A_BSSN3_ll( i, j, k, jyz ), &
                    A_BSSN3_ll( i, j, k, jzz ), &
                    trK( i, j, k ), &
                    phi( i, j, k ), &
                    Gamma_u( i, j, k, 1 ), &
                    Gamma_u( i, j, k, 2 ), &
                    Gamma_u( i, j, k, 3 )
                ELSE
                  WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC( i, j, k ), &
                    MC( i, j, k, jx ), &
                    MC( i, j, k, jy ), &
                    MC( i, j, k, jz ), &
                    GC( i, j, k, jx ), &
                    GC( i, j, k, jy ), &
                    GC( i, j, k, jz )
                ENDIF

                IF( ios > 0 )THEN
                  PRINT *, "...error when writing the arrays in ", &
                           TRIM(namefile), &
                           ". The error message is", err_msg
                  STOP
                ENDIF
                !CALL test_status( ios, err_msg, &
                !          "...error in writing " &
                !          // "the arrays in " // TRIM(namefile) )
              ENDDO
            ENDDO
          ENDDO
        END ASSOCIATE
      ENDDO

      CLOSE( UNIT= 20 )

      PRINT *, " * Printed."
      PRINT *

    ENDIF

    !DEALLOCATE( baryon_density )
    !DEALLOCATE( energy_density )
    !DEALLOCATE( specific_energy )
    !DEALLOCATE( pressure )
    !DEALLOCATE( v_euler )
    !!DEALLOCATE( u_coord )
    !!DEALLOCATE( u_coord_l )
    !DEALLOCATE( g4 )
    !DEALLOCATE( g4temp )
    !DEALLOCATE( ig4 )
    !DEALLOCATE( Tmunu_ll )
    DEALLOCATE( levels )


    CONTAINS


    SUBROUTINE compute_4velocity_eul

      !**************************************************
      !
      !# Compute the components of the fluid \(4\)-velocity
      !  wrt the Eulerian observer
      !
      !  FT 25.04.2022
      !
      !**************************************************

      IMPLICIT NONE

#ifdef __INTEL_COMPILER

  ASSOCIATE( v_euler_l      => v_euler_l% levels(l)% var, &
             u_euler_l      => u_euler_l% levels(l)% var, &
             v_euler        => v_euler% levels(l)% var, &
             lorentz_factor => lorentz_factor% levels(l)% var, &
             lapse          => this% lapse% levels(l)% var, &
             shift_u        => this% shift_u% levels(l)% var, &
             g_phys3_ll     => this% g_phys3_ll% levels(l)% var, &
             g4             => g4% levels(l)% var, &
             Tmunu_ll       => Tmunu_ll% levels(l)% var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP          SHARED( this, v_euler_l, u_euler_l, lorentz_factor, &
  !$OMP                  v_euler, Tmunu_ll, energy_density, pressure, &
  !$OMP                  show_progress, l ) &
  !$OMP          PRIVATE( i, j, k, g4, detg4, g4temp, ig4, u_euler_norm, &
  !$OMP                   perc )

#endif
        DO k= 1, this% get_ngrid_z(l), 1
          DO j= 1, this% get_ngrid_y(l), 1
            DO i= 1, this% get_ngrid_x(l), 1

#ifdef __GFORTRAN__

  ASSOCIATE( v_euler_l      => v_euler_l% levels(l)% var, &
             u_euler_l      => u_euler_l% levels(l)% var, &
             v_euler        => v_euler% levels(l)% var, &
             lorentz_factor => lorentz_factor% levels(l)% var, &
             lapse          => this% lapse% levels(l)% var, &
             shift_u        => this% shift_u% levels(l)% var, &
             g_phys3_ll     => this% g_phys3_ll% levels(l)% var, &
             g4             => g4% levels(l)% var, &
             Tmunu_ll       => Tmunu_ll% levels(l)% var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

#endif

              !energy_density( i, j, k )= baryon_density( i, j, k ) &
              !                            + ( specific_energy(i,j,k) + 1.0 ) &
              !                                 *baryon_density( i, j, k )

              v_euler_l(i,j,k,jx)= g_phys3_ll(i,j,k,jxx)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jxy)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jxz)*v_euler(i,j,k,jz)
              v_euler_l(i,j,k,jy)= g_phys3_ll(i,j,k,jxy)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jyy)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jyz)*v_euler(i,j,k,jz)
              v_euler_l(i,j,k,jz)= g_phys3_ll(i,j,k,jxz)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jyz)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jzz)*v_euler(i,j,k,jz)

              lorentz_factor( i, j, k )= one/SQRT( one &
                              - ( v_euler_l(i,j,k,jx)*v_euler(i,j,k,jx) &
                                + v_euler_l(i,j,k,jy)*v_euler(i,j,k,jy) &
                                + v_euler_l(i,j,k,jz)*v_euler(i,j,k,jz) ) )


              u_euler_l(i,j,k,it)= lorentz_factor( i, j, k ) &
                 *( - lapse( i, j, k ) &
                    + v_euler_l( i, j, k, jx )*shift_u( i, j, k, jx ) &
                    + v_euler_l( i, j, k, jy )*shift_u( i, j, k, jy ) &
                    + v_euler_l( i, j, k, jz )*shift_u( i, j, k, jz ) )
              u_euler_l(i,j,k,ix)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jx )
              u_euler_l(i,j,k,iy)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jy )
              u_euler_l(i,j,k,iz)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jz )

              CALL compute_g4( lapse(i,j,k), shift_u(i,j,k,:), &
                               g_phys3_ll(i,j,k,:), g4(i,j,k,:) )

              CALL determinant_sym4x4( g4(i,j,k,:), detg4 )

              IF( ABS( detg4 ) < 1.0D-10 )THEN
                  PRINT *, "The determinant of the spacetime metric "&
                           // "is effectively 0 at the grid point " &
                           // "(i,j,k)= (", i, ",", j, ",", k, &
                              ")."
                  PRINT *, "detg4=", detg4
                  PRINT *
                  STOP
              ELSEIF( detg4 > zero )THEN
                  PRINT *, "The determinant of the spacetime metric "&
                           // "is positive at the grid point " &
                           // "(i,j,k)= (", i, ",", j, ",", k, &
                              ")."
                  PRINT *, "detg4=", detg4
                  PRINT *
                  STOP
              ENDIF

              g4temp(1,1)= g4(i,j,k,itt)
              g4temp(1,2)= g4(i,j,k,itx)
              g4temp(1,3)= g4(i,j,k,ity)
              g4temp(1,4)= g4(i,j,k,itz)

              g4temp(2,1)= g4(i,j,k,itx)
              g4temp(2,2)= g4(i,j,k,ixx)
              g4temp(2,3)= g4(i,j,k,ixy)
              g4temp(2,4)= g4(i,j,k,ixz)

              g4temp(3,1)= g4(i,j,k,ity)
              g4temp(3,2)= g4(i,j,k,ixy)
              g4temp(3,3)= g4(i,j,k,iyy)
              g4temp(3,4)= g4(i,j,k,iyz)

              g4temp(4,1)= g4(i,j,k,itz)
              g4temp(4,2)= g4(i,j,k,ixz)
              g4temp(4,3)= g4(i,j,k,iyz)
              g4temp(4,4)= g4(i,j,k,izz)

              CALL invert_4x4_matrix( g4temp, ig4 )

              u_euler_norm= ig4(it,it)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,it) &
                          + two*ig4(it,ix)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,ix) &
                          + two*ig4(it,iy)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iy) &
                          + two*ig4(it,iz)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iz) &
                          + ig4(ix,ix)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,ix) &
                          + two*ig4(ix,iy)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iy) &
                          + two*ig4(ix,iz)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iz) &
                          + ig4(iy,iy)* &
                            u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iy) &
                          + two*ig4(iy,iz)* &
                            u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iz) &
                          + two*ig4(iz,iz)* &
                            u_euler_l(i,j,k,iz)*u_euler_l(i,j,k,iz)

              IF( ABS( u_euler_norm + one ) > 1.0D-4 )THEN
                  PRINT *, "** ERROR! The fluid 4-velocity in the " &
                           // "coordinate frame does not have norm -1. " &
                           // "The norm is", u_euler_norm
                  STOP
              ENDIF

#ifdef __GFORTRAN__
  END ASSOCIATE
#endif

            ENDDO
          ENDDO
        ENDDO

#ifdef __INTEL_COMPILER
  !$OMP END PARALLEL DO
  END ASSOCIATE
#endif


    END SUBROUTINE compute_4velocity_eul


    SUBROUTINE compute_stress_energy

      !**************************************************
      !
      !# Compute the components of the stress-energy tensor
      !
      !  FT 25.04.2022
      !
      !**************************************************

      IMPLICIT NONE

      Tmunu_ll% levels(l)% var= zero

#ifdef __INTEL_COMPILER

  ASSOCIATE( v_euler_l      => v_euler_l% levels(l)% var, &
             u_euler_l      => u_euler_l% levels(l)% var, &
             v_euler        => v_euler% levels(l)% var, &
             lorentz_factor => lorentz_factor% levels(l)% var, &
             lapse          => this% lapse% levels(l)% var, &
             shift_u        => this% shift_u% levels(l)% var, &
             g_phys3_ll     => this% g_phys3_ll% levels(l)% var, &
             g4             => g4% levels(l)% var, &
             Tmunu_ll       => Tmunu_ll% levels(l)% var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP          SHARED( this, v_euler_l, u_euler_l, lorentz_factor, &
  !$OMP                  v_euler, Tmunu_ll, energy_density, pressure, &
  !$OMP                  show_progress, l ) &
  !$OMP          PRIVATE( i, j, k, g4, detg4, g4temp, ig4, u_euler_norm, &
  !$OMP                   perc )

#endif
              DO k= 1, this% get_ngrid_z(l), 1
                DO j= 1, this% get_ngrid_y(l), 1
                  DO i= 1, this% get_ngrid_x(l), 1

#ifdef __GFORTRAN__

  ASSOCIATE( v_euler_l      => v_euler_l% levels(l)% var, &
             u_euler_l      => u_euler_l% levels(l)% var, &
             v_euler        => v_euler% levels(l)% var, &
             lorentz_factor => lorentz_factor% levels(l)% var, &
             lapse          => this% lapse% levels(l)% var, &
             shift_u        => this% shift_u% levels(l)% var, &
             g_phys3_ll     => this% g_phys3_ll% levels(l)% var, &
             g4             => g4% levels(l)% var, &
             Tmunu_ll       => Tmunu_ll% levels(l)% var, &
             energy_density => energy_density% levels(l)% var, &
             pressure       => pressure% levels(l)% var &
  )

#endif

                    Tmunu_ll(i,j,k,itt)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,it) &
                            + pressure(i,j,k)*g4(i,j,k,itt) &
                             )

                    Tmunu_ll(i,j,k,itx)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,ix) &
                            + pressure(i,j,k)*g4(i,j,k,itx) &
                             )

                    Tmunu_ll(i,j,k,ity)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iy) &
                            + pressure(i,j,k)*g4(i,j,k,ity) &
                             )

                    Tmunu_ll(i,j,k,itz)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iz) &
                            + pressure(i,j,k)*g4(i,j,k,itz) &
                             )

                    Tmunu_ll(i,j,k,ixx)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,ix) &
                            + pressure(i,j,k)*g4(i,j,k,ixx) &
                             )

                    Tmunu_ll(i,j,k,ixy)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iy) &
                            + pressure(i,j,k)*g4(i,j,k,ixy) &
                             )

                    Tmunu_ll(i,j,k,ixz)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iz) &
                            + pressure(i,j,k)*g4(i,j,k,ixz) &
                             )

                    Tmunu_ll(i,j,k,iyy)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iy) &
                            + pressure(i,j,k)*g4(i,j,k,iyy)  &
                             )

                    Tmunu_ll(i,j,k,iyz)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iz) &
                            + pressure(i,j,k)*g4(i,j,k,iyz) &
                             )

                    Tmunu_ll(i,j,k,izz)= lorene2hydrobase*( &
                            ( energy_density(i,j,k) + pressure(i,j,k) ) &
                            *u_euler_l(i,j,k,iz)*u_euler_l(i,j,k,iz) &
                            + pressure(i,j,k)*g4(i,j,k,izz) &
                             )

#ifdef __GFORTRAN__
  END ASSOCIATE
#endif

                  ENDDO
                ENDDO
              ENDDO

#ifdef __INTEL_COMPILER
  !$OMP END PARALLEL DO
  END ASSOCIATE
#endif

    END SUBROUTINE compute_stress_energy


  END PROCEDURE compute_and_export_bssn_constraints_grid


  MODULE PROCEDURE compute_and_export_bssn_constraints_particles

    !**************************************************
    !
    !# Compute, store and print the BSSN constraints
    !  to a formatted file. The computaton is done
    !  mapping the physical metric from the gravity
    !  to the particles, computing e stress-energy
    !  tensor on the particles, and mapping it to the
    !  gravity grid.
    !  @todo use the SPH density to compute the
    !       stress-energy tensor, rather than the
    !       density from the |id|
    !
    !  FT 1.02.2021
    !
    !**************************************************

    USE units,                ONLY: set_units
    USE tensor,               ONLY: itt, itx, ity, itz, ixx, ixy, &
                                    ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                    jyy, jyz, jzz, jx, jy, jz

    USE mesh_refinement,             ONLY: allocate_grid_function, levels, &
                                           rad_coord, nlevels, &
                                           deallocate_grid_function, coords
    USE ADM_refine,                  ONLY: lapse, shift_u, &
                                           g_phys3_ll, &
                                           allocate_ADM, deallocate_ADM
    USE BSSN_refine,                 ONLY: allocate_BSSN, deallocate_BSSN
    USE Tmunu_refine,                ONLY: Tmunu_ll, allocate_Tmunu, &
                                           deallocate_Tmunu
    USE McLachlan_refine,            ONLY: BSSN_CONSTRAINTS_INTERIOR, &
                                           allocate_Ztmp, deallocate_Ztmp
    USE GravityAcceleration_refine,  ONLY: allocate_GravityAcceleration, &
                                           deallocate_GravityAcceleration


    USE input_output,         ONLY: read_options
    USE options,              ONLY: ndes
    USE sph_variables,        ONLY: npart, &  ! particle number
                                    pos_u, &  ! particle positions
                                    vel_u, &  ! particle velocities in
                                              ! coordinate frame
                                    nlrf,  &  ! baryon number density in
                                              ! local rest frame
                                    !ehat,  &  ! canonical energy per baryon
                                    nu,    &  ! canonical baryon number per
                                              ! particle
                                    Theta, &  ! Generalized Lorentz factor
                                    h,     &  ! Smoothing length
                                    Pr,    &  ! Pressure
                                    u,     &  ! Internal energy in local rest
                                              ! frame (no kinetic energy)
                                    !temp,  &  ! Temperature
                                    !av,    &  ! Dissipation
                                    !Ye,    &  ! Electron fraction
                                    !divv,  &  ! Divergence of velocity vel_u
                                    !Nstar, &  ! Comput.frame baryon number
                                              ! density
                                    allocate_SPH_memory, &
                                    deallocate_SPH_memory
    USE RCB_tree_3D,          ONLY: allocate_RCB_tree_memory_3D,&
                                    deallocate_RCB_tree_memory_3D, iorig
    USE kernel_table,         ONLY: ktable
    USE gradient,             ONLY: allocate_gradient, deallocate_gradient
    USE sphincs_sph,          ONLY: density, flag_dead_ll_cells
    USE set_h,                ONLY: exact_nei_tree_update
    USE alive_flag,           ONLY: alive

    USE map_particles_2_grid, ONLY: map_2_grid_hash
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    get_metric_on_particles
    !USE particle_mesh,        ONLY: all_lists, flag_nei_cell, pp_g
    USE particle_mesh_hash,   ONLY: deallocate_hash_memory

    IMPLICIT NONE

    INTEGER:: i, j, k, l, a, allocation_status
    INTEGER, DIMENSION(3) :: imin, imax
    INTEGER:: unit_logfile, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z
    INTEGER, SAVE:: counter= 1


    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_loc
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_loc
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: theta_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sph_density

    CHARACTER( LEN= : ), ALLOCATABLE:: name_constraint
    CHARACTER( LEN= : ), ALLOCATABLE:: name_analysis
    CHARACTER( LEN= : ), ALLOCATABLE:: finalname_logfile
    CHARACTER( 2 ):: n_reflev

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    ALLOCATE ( levels( this% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    nlevels= this% nlevels
    levels = this% levels
    coords = this% coords

    DO l= 1, this% nlevels, 1
      levels(l)% ngrid_x= this% levels(l)% ngrid_x
      levels(l)% ngrid_x= this% levels(l)% ngrid_x
      levels(l)% ngrid_x= this% levels(l)% ngrid_x
    ENDDO

    IF( debug ) PRINT *, "ngrid_x=", this% levels(1)%ngrid_x
    IF( debug ) PRINT *, "ngrid_y=", this% levels(1)%ngrid_y
    IF( debug ) PRINT *, "ngrid_z=", this% levels(1)%ngrid_z
    IF( debug ) PRINT *

    CALL allocate_grid_function( this% HC_parts, "HC_parts_ID", 1 )
    CALL allocate_grid_function( this% MC_parts, "MC_parts_ID", 3 )
    CALL allocate_grid_function( this% GC_parts, "GC_parts_ID", 3 )

    CALL allocate_grid_function( this% rho_parts, "rho_parts", 1 )
    CALL allocate_grid_function( this% S_parts, "S_parts", 3 )

    PRINT *, "Mapping hydro fields from particles to grid..."

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    CALL allocate_grid_function( rad_coord, 'rad_coord', 1 )

    ! Initialize the stress-energy tensor to 0
    DO l= 1, this% nlevels, 1
      Tmunu_ll%   levels(l)% var= zero
      rad_coord%  levels(l)% var= this% rad_coord%  levels(l)% var
      g_phys3_ll% levels(l)% var= this% g_phys3_ll% levels(l)% var
      shift_u%    levels(l)% var= this% shift_u%    levels(l)% var
      lapse%      levels(l)% var= this% lapse%      levels(l)% var
    ENDDO

    npart= parts_obj% get_npart()

    IF( .NOT. ALLOCATED( nlrf_loc ) )THEN
        ALLOCATE( nlrf_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for nlrf_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( nu_loc ) )THEN
        ALLOCATE( nu_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for nu_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( u_loc ) )THEN
        ALLOCATE( u_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for u_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( pressure_loc ) )THEN
        ALLOCATE( pressure_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for pressure_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( theta_loc ) )THEN
        ALLOCATE( theta_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for theta_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( pos_loc ) )THEN
        ALLOCATE( pos_loc( 3, npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for pos_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( vel_loc ) )THEN
        ALLOCATE( vel_loc( 3, npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for vel_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( sph_density ) )THEN
        ALLOCATE( sph_density( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for sph_density'
       STOP
    ENDIF

    ! Set the SPH density to 0 by default
    sph_density= zero

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory

    IF( debug ) PRINT *, "-2"

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    IF( debug ) PRINT *, "-1"

    h= parts_obj% get_h()

    !IF( counter == 1 )THEN
    !  ! tabulate kernel, get ndes
    !  CALL ktable(ikernel,ndes)
    !ENDIF

    IF( debug ) PRINT *, "0"

    nu_loc      = parts_obj% get_nu()
    pos_loc     = parts_obj% get_pos()
    vel_loc     = parts_obj% get_vel()
    u_loc       = parts_obj% get_u()
    nlrf_loc    = parts_obj% get_nlrf()
    theta_loc   = parts_obj% get_theta()
    pressure_loc= parts_obj% get_pressure_cu()

    IF( debug ) PRINT *, "1"

    PRINT *, " * Allocating needed memory..."
    PRINT *

    ! flag that particles are 'alive'
    IF( .NOT.ALLOCATED( alive ) ) ALLOCATE( alive( npart ) )
    alive( 1:npart )= 1

    CALL allocate_gradient( npart )

    IF( debug ) PRINT *, "2"

    CALL allocate_metric_on_particles( npart )

    IF( debug ) PRINT *, "3"

    !---------------------------!
    !--  Compute constraints  --!
    !---------------------------!

    PRINT *, " * Mapping metric from the grid to the particles..."
    PRINT *
    CALL get_metric_on_particles( npart, &
                                  pos_loc )

    IF( debug ) PRINT *, "4"

    !
    !-- Seems like computing neighbors and SPH density is not needed to map
    !-- the stress-energy tensor from the particles to the grid
    !

    PRINT *, " * Computing neighbours..."
    PRINT *
    CALL exact_nei_tree_update( ndes,    &
                                npart,   &
                                pos_loc, &
                                nu_loc )

    !IF( debug ) PRINT *, "5"
    !
    !PRINT *, " * Computing SPH density..."
    !PRINT *
    nu   = nu_loc
    pos_u= pos_loc
    vel_u= vel_loc
    u    = u_loc
    nlrf = nlrf_loc
    Theta= theta_loc
    Pr   = pressure_loc
    !CALL density( npart,   &
    !              pos_loc, &
    !              sph_density )

    IF( debug ) PRINT *, "6"

    IF( debug .AND. .TRUE. ) PRINT *, "npart= ", npart
    IF( debug .AND. .TRUE. ) PRINT *, "nu_loc= ", nu_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "pos_loc= ", pos_loc(2,npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "vel_loc= ", vel_loc(2,npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "u_loc= ", u_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "nlrf_loc= ", nlrf_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "theta_loc= ", theta_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "pressure_loc= ", pressure_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *

    !IF( counter == 2 ) STOP

    PRINT *, " * Mapping stress-energy tensor from the particles to the grid..."
    PRINT *
    CALL map_2_grid_hash( npart        , &
                          nu_loc       , &
                          pos_loc      , &
                          vel_loc      , &
                          u_loc        , &
                          nlrf_loc     , &
                          theta_loc    , &
                          pressure_loc )

    !IF( counter == 2 )THEN
    !  STOP
    !ENDIF

    IF( debug ) PRINT *, "6.5"

    !IF( counter == 2 ) STOP

    IF( debug ) PRINT *, "7"

    !
    !-- Deallocate SPH MODULE variables
    !
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )
    !IF( ALLOCATED( flag_nei_cell% levels ) )THEN
    !  CALL deallocate_grid_function( flag_nei_cell, 'flag_nei_cell' )
    !ENDIF
    !IF( ALLOCATED( pp_g% levels ) )THEN
    !  CALL deallocate_grid_function( pp_g, 'pp_g' )
    !ENDIF
    !IF( ALLOCATED(all_lists% levels) )THEN
    !  CALL deallocate_grid_function( all_lists, 'all_lists' )
    !ENDIF
    CALL deallocate_hash_memory
    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE( alive )
    !DEALLOCATE(W_no_norm)
    !DEALLOCATE(dWdv_no_norm)
    !DEALLOCATE(fmass)
    !DEALLOCATE(fpoten)
    !DEALLOCATE(dphidh)
    CALL deallocate_RCB_tree_memory_3D
    CALL deallocate_SPH_memory

    IF( debug ) PRINT *, "8.1"

    !
    !-- Compute the BSSN constraints by calling the Cactus-bound procedure
    !-- BSSN_CONSTRAINTS_INTERIOR
    !
    PRINT *, " * Computing constraints using particle data..."
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP          SHARED( this, Tmunu_ll ) &
!    !$OMP          PRIVATE( l, imin, imax )
    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse      => this% lapse% levels(l)% var, &
                 shift_u    => this% shift_u% levels(l)% var, &
                 phi        => this% phi% levels(l)% var, &
                 trK        => this% trK% levels(l)% var, &
                 g_BSSN3_ll => this% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => this% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u    => this% Gamma_u% levels(l)% var, &
                 Tmunu_ll   => Tmunu_ll% levels(l)% var, &
                 HC_parts   => this% HC_parts % levels(l)% var, &
                 MC_parts   => this% MC_parts % levels(l)% var, &
                 GC_parts   => this% GC_parts % levels(l)% var, &
                 rho_parts  => this% rho_parts% levels(l)% var, &
                 S_parts    => this% S_parts% levels(l)% var &
      )

        imin(1) = this% levels(l)% nghost_x
        imin(2) = this% levels(l)% nghost_y
        imin(3) = this% levels(l)% nghost_z
        imax(1) = this% get_ngrid_x(l) - this% levels(l)% nghost_x - 1
        imax(2) = this% get_ngrid_y(l) - this% levels(l)% nghost_y - 1
        imax(3) = this% get_ngrid_z(l) - this% levels(l)% nghost_z - 1

        HC_parts    = zero
        MC_parts    = zero
        GC_parts    = zero
        rho_parts   = zero
        S_parts     = zero
        CALL bssn_constraint_terms_interior( &
          !
          !-- Input
          !
          this% get_ngrid_x(l), this% get_ngrid_y(l), this% get_ngrid_z(l), &
          imin, imax, &
          this% get_dx(l), this% get_dy(l), this% get_dz(l), &
          g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
          g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
          g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
          A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
          A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
          A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
          trK(:,:,:), phi(:,:,:), &
          Gamma_u(:,:,:,jx), &
          Gamma_u(:,:,:,jy), &
          Gamma_u(:,:,:,jz), &
          Tmunu_ll(:,:,:,itt), &
          Tmunu_ll(:,:,:,itx), &
          Tmunu_ll(:,:,:,ity), &
          Tmunu_ll(:,:,:,itz), &
          Tmunu_ll(:,:,:,ixx), &
          Tmunu_ll(:,:,:,ixy), &
          Tmunu_ll(:,:,:,ixz), &
          Tmunu_ll(:,:,:,iyy), &
          Tmunu_ll(:,:,:,iyz), &
          Tmunu_ll(:,:,:,izz), &
          lapse(:,:,:), &
          shift_u(:,:,:,jx), &
          shift_u(:,:,:,jy), &
          shift_u(:,:,:,jz), &
          !
          !-- Output
          !
          ! Connection constraints
          GC_parts(:,:,:,jx), &
          GC_parts(:,:,:,jy), &
          GC_parts(:,:,:,jz), &
          ! Hamiltonian and momentum constraints
          HC_parts(:,:,:), &
          MC_parts(:,:,:,jx), &
          MC_parts(:,:,:,jy), &
          MC_parts(:,:,:,jz), &
          ! Sources in the Hamiltonian and momentum constraints
          rho_parts(:,:,:), &
          S_parts(:,:,:,jx), &
          S_parts(:,:,:,jy), &
          S_parts(:,:,:,jz) &
        )
        !CALL BSSN_CONSTRAINTS_INTERIOR( &
        !  !
        !  !-- Input
        !  !
        !  this% get_ngrid_x(l), this% get_ngrid_y(l), this% get_ngrid_z(l), &
        !  imin, imax, &
        !  this% get_dx(l), this% get_dy(l), this% get_dz(l), &
        !  g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
        !  g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
        !  g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
        !  A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
        !  A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
        !  A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
        !  trK(:,:,:), phi(:,:,:), &
        !  Gamma_u(:,:,:,jx), &
        !  Gamma_u(:,:,:,jy), &
        !  Gamma_u(:,:,:,jz), &
        !  Tmunu_ll(:,:,:,itt), &
        !  Tmunu_ll(:,:,:,itx), &
        !  Tmunu_ll(:,:,:,ity), &
        !  Tmunu_ll(:,:,:,itz), &
        !  Tmunu_ll(:,:,:,ixx), &
        !  Tmunu_ll(:,:,:,ixy), &
        !  Tmunu_ll(:,:,:,ixz), &
        !  Tmunu_ll(:,:,:,iyy), &
        !  Tmunu_ll(:,:,:,iyz), &
        !  Tmunu_ll(:,:,:,izz), &
        !  lapse(:,:,:), &
        !  shift_u(:,:,:,jx), &
        !  shift_u(:,:,:,jy), &
        !  shift_u(:,:,:,jz), &
        !  !
        !  !-- Output
        !  !
        !  ! Connection constraints
        !  GC_parts(:,:,:,jx), &
        !  GC_parts(:,:,:,jy), &
        !  GC_parts(:,:,:,jz), &
        !  ! Hamiltonian and momentum constraints
        !  HC_parts(:,:,:), &
        !  MC_parts(:,:,:,jx), &
        !  MC_parts(:,:,:,jy), &
        !  MC_parts(:,:,:,jz) &
        !)
      END ASSOCIATE
    ENDDO
!    !$OMP END PARALLEL DO
    PRINT *, " * Constraints computed."
    PRINT *

    IF( debug ) PRINT *, "0"

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    DO l= 1, this% nlevels, 1

      ASSOCIATE( HC_parts  => this% HC_parts% levels(l)% var, &
                 MC_parts  => this% MC_parts% levels(l)% var, &
                 GC_parts  => this% GC_parts% levels(l)% var, &
                 rho_parts => this% rho_parts% levels(l)% var, &
                 S_parts   => this% S_parts% levels(l)% var &
      )

        unit_logfile= 2791

        IF( l > 9 )THEN
          WRITE( n_reflev, "(I2)" ) l
        ELSE
          WRITE( n_reflev, "(I1)" ) l
        ENDIF

        finalname_logfile= TRIM(name_logfile)//"-reflev"//TRIM(n_reflev)//".log"

        INQUIRE( FILE= TRIM(finalname_logfile), EXIST= exist )

        IF( debug ) PRINT *, "1"

        IF( exist )THEN
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "REPLACE", &
                  FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "NEW", &
                  FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening ", TRIM(finalname_logfile), &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error when opening " &
        !                  // TRIM(name_logfile) )

        IF( debug ) PRINT *, "2"

        IF( .NOT.ALLOCATED( this% HC_parts_l2 ))THEN
          ALLOCATE( this% HC_parts_l2( this% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% MC_parts_l2 ))THEN
          ALLOCATE( this% MC_parts_l2( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% GC_parts_l2 ))THEN
          ALLOCATE( this% GC_parts_l2( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% HC_parts_loo ))THEN
          ALLOCATE( this% HC_parts_loo( this% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% MC_parts_loo ))THEN
          ALLOCATE( this% MC_parts_loo( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( this% GC_parts_loo ))THEN
          ALLOCATE( this% GC_parts_loo( this% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
     !   IF( .NOT.ALLOCATED( this% HC_parts_int ))THEN
     !     ALLOCATE( this% HC_parts_int( this% nlevels ), &
     !               STAT= ios, ERRMSG= err_msg )
     !     IF( ios > 0 )THEN
     !       PRINT *, "...allocation error for array MC_loo. ", &
     !                "The error message is", err_msg
     !       STOP
     !     ENDIF
     !     !CALL test_status( ios, err_msg, &
     !     !                "...deallocation error for array HC" )
     !   ENDIF
     !   IF( .NOT.ALLOCATED( this% MC_parts_int ))THEN
     !     ALLOCATE( this% MC_parts_int( this% nlevels, 3 ), &
     !               STAT= ios, ERRMSG= err_msg )
     !     IF( ios > 0 )THEN
     !       PRINT *, "...allocation error for array MC_loo. ", &
     !                "The error message is", err_msg
     !       STOP
     !     ENDIF
     !     !CALL test_status( ios, err_msg, &
     !     !                "...deallocation error for array HC" )
     !   ENDIF
     !   IF( .NOT.ALLOCATED( this% GC_parts_int ))THEN
     !     ALLOCATE( this% GC_parts_int( this% nlevels, 3 ), &
     !               STAT= ios, ERRMSG= err_msg )
     !     IF( ios > 0 )THEN
     !       PRINT *, "...allocation error for array GC_loo. ", &
     !                "The error message is", err_msg
     !       STOP
     !     ENDIF
     !     !CALL test_status( ios, err_msg, &
     !     !                "...deallocation error for array HC" )
     !   ENDIF

        WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

        PRINT *, "** Analyzing constraints on refinement level ", l, "..."

        name_analysis= "bssn-hc-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the Hamiltonian constraint"
        CALL this% analyze_constraint( &
             l, &
             HC_parts, name_constraint, unit_logfile, name_analysis, &
             this% HC_parts_l2(l), this% HC_parts_loo(l), &
             this% HC_parts_int(l), rho_parts )


        name_analysis= "bssn-mc1-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             this% MC_parts_l2(l,jx), this% MC_parts_loo(l,jx), &
             this% MC_parts_int(l,jx), S_parts(:,:,:,jx) )

        name_analysis= "bssn-mc2-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             this% MC_parts_l2(l,jy), this% MC_parts_loo(l,jy), &
             this% MC_parts_int(l,jy), S_parts(:,:,:,jy) )

        name_analysis= "bssn-mc3-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the momentum constraint"
        CALL this% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             this% MC_parts_l2(l,jz), this% MC_parts_loo(l,jz), &
             this% MC_parts_int(l,jz), S_parts(:,:,:,jz) )

        name_analysis= "bssn-gc1-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             this% GC_parts_l2(l,jx), this% GC_parts_loo(l,jx), &
             this% GC_parts_int(l,jx) )

        name_analysis= "bssn-gc2-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             this% GC_parts_l2(l,jy), this% GC_parts_loo(l,jy), &
             this% GC_parts_int(l,jy) )

        name_analysis= "bssn-gc3-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the connection constraint"
        CALL this% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             this% GC_parts_l2(l,jz), this% GC_parts_loo(l,jz), &
             this% GC_parts_int(l,jz) )

        CLOSE( UNIT= unit_logfile )

        PRINT *, " * Constraints analyzed. Summary of results saved to ", &
                 finalname_logfile
        PRINT *
      END ASSOCIATE
    ENDDO

    IF( this% export_constraints )THEN

      PRINT *, " * Printing constraints to file ", TRIM(namefile), "..."

      !
      !-- Export the constraints to a formatted file
      !

      INQUIRE( FILE= TRIM(namefile), EXIST= exist )

      IF( debug ) PRINT *, "1"

      IF( exist )THEN
          OPEN( UNIT= 21, FILE= TRIM(namefile), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= 21, FILE= TRIM(namefile), STATUS= "NEW", &
          FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(namefile) )

      IF( debug ) PRINT *, "2"

      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the BSSN constraints computed with the mapping routines ", &
      "for the ID on selected grid points"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
      !         // TRIM(namefile) )
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4       5", &
      "       6       7       8       9       10", &
      "       11       12       13       14       15", &
      "       16       17       18       19       20"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 2 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
      !        // TRIM(namefile) )
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      refinement level    x   y   z   Stress-energy (10 components)   "&
      // "Hamiltonian constraint       " &
      // "Momentum constraint (three components)       " &
      // "Connection constraint (three components)"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 3 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
      !        // TRIM(namefile) )

      IF( debug ) PRINT *, "3"

      DO l= 1, this% nlevels, 1

        ASSOCIATE( lapse           => this% lapse% levels(l)% var, &
                   shift_u         => this% shift_u% levels(l)% var, &
                   phi             => this% phi% levels(l)% var, &
                   trK             => this% trK% levels(l)% var, &
                   g_BSSN3_ll      => this% g_BSSN3_ll% levels(l)% var, &
                   A_BSSN3_ll      => this% A_BSSN3_ll% levels(l)% var, &
                   g_phys3_ll      => this% g_phys3_ll% levels(l)% var, &
                   k_phys3_ll      => this% k_phys3_ll% levels(l)% var, &
                   Gamma_u         => this% Gamma_u% levels(l)% var, &
                   Tmunu_ll        => Tmunu_ll% levels(l)% var, &
                   HC_parts        => this% HC_parts% levels(l)% var, &
                   MC_parts        => this% MC_parts% levels(l)% var, &
                   GC_parts        => this% GC_parts% levels(l)% var &
        )

          ! Being abs_grid a local array, it is good practice to allocate it on
          ! the heap, otherwise it will be stored on the stack which has a very
          ! limited size. This results in a segmentation fault.
          IF( ALLOCATED( abs_grid ) )THEN
            DEALLOCATE( abs_grid )
          ENDIF
          ALLOCATE( abs_grid( this% get_ngrid_x(l), this% get_ngrid_y(l), &
                              this% get_ngrid_z(l), 3 ) )

          DO k= 1, this% get_ngrid_z(l), 1
            DO j= 1, this% get_ngrid_y(l), 1
              DO i= 1, this% get_ngrid_x(l), 1

                abs_grid( i, j, k, jx )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jx ) )
                abs_grid( i, j, k, jy )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                abs_grid( i, j, k, jz )= &
                            ABS( this% coords% levels(l)% var( i, j, k, jz ) )

              ENDDO
            ENDDO
          ENDDO

          min_abs_y= 1D+20
          min_abs_z= 1D+20
          DO k= 1, this% get_ngrid_z(l), 1
            DO j= 1, this% get_ngrid_y(l), 1
              DO i= 1, this% get_ngrid_x(l), 1

                IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                    < min_abs_y )THEN
                  min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                  min_ix_y= i
                  min_iy_y= j
                  min_iz_y= k
                ENDIF

                IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                    < min_abs_z )THEN
                  min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                  min_ix_z= i
                  min_iy_z= j
                  min_iz_z= k
                ENDIF

              ENDDO
            ENDDO
          ENDDO

          DO k= 1, this% get_ngrid_z(l), 1

            IF( MOD( k, this% cons_step ) /= 0 ) CYCLE

            DO j= 1, this% get_ngrid_y(l), 1

              IF( MOD( j, this% cons_step ) /= 0 ) CYCLE

              DO i= 1, this% get_ngrid_x(l), 1

                IF( MOD( i, this% cons_step ) /= 0 ) CYCLE

                IF( this% export_constraints_xy .AND. &
                    ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                      this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) ) )THEN
                  CYCLE
                ENDIF
                IF( this% export_constraints_x .AND. &
                    ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                      this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) &
                      .OR. &
                      this% coords% levels(l)% var( i, j, k, jy ) /= &
                      this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                    min_iz_y, jy ) ) )THEN
                  CYCLE
                ENDIF

                IF( debug )THEN
                  WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), & ! columns 18
                    !pos_loc( 1, ix, iy, iz ), &
                    !pos_loc( 2, ix, iy, iz ), &
                    !pos_loc( 3, ix, iy, iz ), &
                    !vel_loc( 1, ix, iy, iz ), &
                    !vel_loc( 2, ix, iy, iz ), &
                    !vel_loc( 3, ix, iy, iz ), &
                    !nu_loc( ix, iy, iz ), &
                    !u_loc( ix, iy, iz ), &
                    !nlrf_loc( ix, iy, iz ), &
                    !theta_loc( ix, iy, iz ), &
                    !pressure_loc( ix, iy, iz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC_parts( i, j, k ), &
                    MC_parts( i, j, k, jx ), &
                    MC_parts( i, j, k, jy ), &
                    MC_parts( i, j, k, jz ), &
                    GC_parts( i, j, k, jx ), &
                    GC_parts( i, j, k, jy ), &
                    GC_parts( i, j, k, jz ), &
                    lapse( i, j, k ), &
                    shift_u( i, j, k, jx ), &
                    shift_u( i, j, k, jy ), &
                    shift_u( i, j, k, jz ), &
                    g_BSSN3_ll( i, j, k, jxx ), &
                    g_BSSN3_ll( i, j, k, jxy ), &
                    g_BSSN3_ll( i, j, k, jxz ), &
                    g_BSSN3_ll( i, j, k, jyy ), &
                    g_BSSN3_ll( i, j, k, jyz ), &
                    g_BSSN3_ll( i, j, k, jzz ), &
                    k_phys3_ll( i, j, k, jxx ), &
                    k_phys3_ll( i, j, k, jxy ), &
                    k_phys3_ll( i, j, k, jxz ), &
                    k_phys3_ll( i, j, k, jyy ), &
                    k_phys3_ll( i, j, k, jyz ), &
                    k_phys3_ll( i, j, k, jzz ), &
                    A_BSSN3_ll( i, j, k, jxx ), &
                    A_BSSN3_ll( i, j, k, jxy ), &
                    A_BSSN3_ll( i, j, k, jxz ), &
                    A_BSSN3_ll( i, j, k, jyy ), &
                    A_BSSN3_ll( i, j, k, jyz ), &
                    A_BSSN3_ll( i, j, k, jzz ), &
                    trK( i, j, k ), &
                    phi( i, j, k ), &
                    Gamma_u( i, j, k, 1 ), &
                    Gamma_u( i, j, k, 2 ), &
                    Gamma_u( i, j, k, 3 )
                ELSE
                  WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    this% coords% levels(l)% var( i, j, k, jx ), &
                    this% coords% levels(l)% var( i, j, k, jy ), &
                    this% coords% levels(l)% var( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC_parts( i, j, k ), &
                    MC_parts( i, j, k, jx ), &
                    MC_parts( i, j, k, jy ), &
                    MC_parts( i, j, k, jz ), &
                    GC_parts( i, j, k, jx ), &
                    GC_parts( i, j, k, jy ), &
                    GC_parts( i, j, k, jz )
                ENDIF

                IF( ios > 0 )THEN
                  PRINT *, "...error when writing the arrays in ", &
                           TRIM(namefile), &
                           ". The error message is", err_msg
                  STOP
                ENDIF
                !CALL test_status( ios, err_msg, &
                !                  "...error in writing " &
                !                  // "the arrays in " // TRIM(namefile) )
              ENDDO
            ENDDO
          ENDDO
        END ASSOCIATE
      ENDDO

      IF( debug ) PRINT *, "4"

      CLOSE( UNIT= 21 )

      PRINT *, " * Printed."
      PRINT *

    ENDIF

    !
    !-- Deallocate spacetime MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    !CALL deallocate_gravity_grid()
    DEALLOCATE( levels )

    ! Count the number of times that this SUBROUTINE is called, since the
    ! kernel has to be tabulated only once in the present implementation
    counter= counter+ 1

  END PROCEDURE compute_and_export_bssn_constraints_particles


END SUBMODULE constraints
