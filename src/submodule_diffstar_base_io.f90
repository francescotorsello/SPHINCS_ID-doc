! File:         submodule_diffstar_base_io.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_base) diffstar_base_io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE diffstarbase that handle I/O (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !------------------------------!
  !--  OVERRIDING SUBROUTINES  --!
  !------------------------------!


  MODULE PROCEDURE print_summary_drs

    !************************************************
    !
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted
    !  file whose name is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    USE constants,      ONLY: lorene2hydrobase, MSun_geo, kg2g, m2cm

    IMPLICIT NONE

    PRINT *, " * Differentially rotating star (DRS):"
    PRINT *
    PRINT *, "   Baryon mass of the DRS=", THIS% mass, "Msun"
    PRINT *, "   Gravitational mass of the DRS=", THIS% mass_grav, "Msun"
    PRINT *
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=0 = ", &
                 THIS% r_eq, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=pi/2 = ", &
                 THIS% r_eq_pi2, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=pi = ", &
                 THIS% r_eq_pi, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=3pi/2 = ", &
                 THIS% r_eq_3pi2, "Msun_geo"
    PRINT *, "   Polar radius of the DRS= ", THIS% r_pole, "Msun_geo"
    PRINT *, "   Ratio between polar radius at equatiorial radius at phi=0= ", &
                 THIS% r_ratio, "Msun_geo"
    PRINT *
    PRINT *, "   EOS for the DRS= ", THIS% eos
    PRINT *
    PRINT *, "   Central baryon mass density for the DRS= ", &
                 THIS% rho_center, "Msun/Msun_geo**3= ", &
                 THIS% rho_center/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
    PRINT *


  END PROCEDURE print_summary_drs


END SUBMODULE diffstar_base_io
