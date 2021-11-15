! File:         submodule_formul_3p1_io.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_3p1_id) formul_3p1_io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE formul_3p1 that handle I/O (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary

    !************************************************
    !
    !# Prints a summary of the properties of the refined mesh,
    !  and optionally, to a formatted file whose name
    !  is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: l, last_level, i_matter

    last_level= THIS% get_nlevels()

    PRINT *, " * Spacetime:"
    PRINT *
    PRINT *, "   Number of refinement levels= ", last_level
    PRINT *
    PRINT *, "   Number of grid points on each level= ", &
             THIS% get_ngrid_x( 1 ), "**3"
    PRINT *
    DO l= 1, last_level, 1
      PRINT *, "   Resolution on level ", l, "= ", THIS% get_dx(l)
    ENDDO
    PRINT *
    DO l= 1, last_level, 1
      PRINT *, "   x boundary of level ", l, "= ", THIS% get_xR(l)
      PRINT *, "   y boundary of level ", l, "= ", THIS% get_yR(l)
      PRINT *, "   z boundary of level ", l, "= ", THIS% get_zR(l)
    ENDDO
    PRINT *
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Number of grid points across the x-axis-diameter of ", &
               "matter object ", i_matter, "=", THIS% npoints_xaxis(i_matter)
    ENDDO
    PRINT *

  END PROCEDURE print_summary


END SUBMODULE formul_3p1_io
