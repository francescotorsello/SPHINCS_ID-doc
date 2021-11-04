!& File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) spherical_surfaces

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method of TYPE particles that places
  !  particles on spherical surfaces inside
  !  a star
  !
  !  FT 19.04.2021
  !
  !************************************************


  IMPLICIT NONE

  ! Be careful! if you define quantities here, they will be global
  ! If you call the SUBROUTINES multiple times, they will use the SAME variables

  !PRIVATE


  !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii


  CONTAINS


  MODULE PROCEDURE place_particles_spherical_surfaces

    !**********************************************
    !
    !# Places particles on spherical surfaces
    !  inside a star
    !
    !  FT 19.04.2021
    !
    !**********************************************

    !$ USE OMP_LIB
    USE constants, ONLY: pi, MSun, MSun_geo, km2m, kg2g, lorene2hydrobase, &
                         golden_ratio, third, half, amu, g2kg, sixth
    USE matrix,    ONLY: determinant_4x4_matrix
    USE NR,        ONLY: indexx
    USE APM,       ONLY: assign_h

    IMPLICIT NONE

    INTEGER:: n_shells, itr2, cnt, &
              r, th, phi, i_shell, npart_test, npart_shell_tmp, &
              cnt2, rel_sign, dim_seed, r_cnt, prev_shell, &
              npart_discard, npart_shell_cnt, size_pos_shell
    !INTEGER, PARAMETER:: max_length= 5D+6
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx, seed
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq

    DOUBLE PRECISION:: xtemp, ytemp, ztemp, m_p, &
                       dr, dth, dphi, phase, phase_th, mass, &
                       dr_shells, dth_shells, dphi_shells, col, long, rad, &
                       proper_volume, mass_test, mass_test2,&
                       proper_volume_test, npart_shell_kept, &
                       rand_num, rand_num2, delta_r, shell_thickness, &
                       upper_bound_tmp, lower_bound_tmp, col_tmp

    DOUBLE PRECISION, PARAMETER:: huge_real= 1.0D30!ABS( HUGE(0.0D0) )

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii, shell_masses, &
                                                  alpha, m_parts, vol_shell, &
                                                  vol_shell2, mass_shell, &
                                                  mass_shell2, shell_scales

    LOGICAL:: exist, high_mass, low_mass, kept_all

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile, finalnamefile2

    TYPE:: colatitude_pos_shell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: colatitudes
    END TYPE

    TYPE(colatitude_pos_shell), DIMENSION(:), ALLOCATABLE:: colatitude_pos

    TYPE:: pos_on_shells
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_shell
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_th
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_phi
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_shell
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_shell2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xx
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: gamma_euler
    END TYPE

    TYPE(pos_on_shells), DIMENSION(:), ALLOCATABLE:: pos_shells

    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: pos_shell_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: g_xx_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: bar_density_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: gam_euler_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pvol_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: npart_surface_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: npart_discarded

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, "** Executing the place_particles_spherical_surfaces SUBROUTINE..."
    PRINT *

    CALL RANDOM_SEED( SIZE= dim_seed )
    ALLOCATE( seed( dim_seed ) )
    seed( 1 )= 0
    seed( 2 )= 1
    DO itr= 3, dim_seed
      seed( itr )= seed( itr - 1 ) + seed( itr - 2 )
    ENDDO
    CALL RANDOM_SEED( PUT= seed )

    !-----------------------------------!
    !-- Compute desired particle mass --!
    !-----------------------------------!

    m_p= mass_star/npart_approx

    !------------------------------------------!
    !-- Compute number of spherical surfaces --!
    !------------------------------------------!

    n_shells= number_surfaces( m_p, center, radius, get_density )

    !------------------------------------------------!
    !-- Allocate memory for the spherical surfaces --!
    !------------------------------------------------!

    IF(.NOT.ALLOCATED( shell_radii ))THEN
      ALLOCATE( shell_radii( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_radii in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( shell_masses ))THEN
      ALLOCATE( shell_masses( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_masses in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( shell_scales ))THEN
      ALLOCATE( shell_scales( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_scales in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( vol_shell ))THEN
      ALLOCATE( vol_shell( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( vol_shell2 ))THEN
      ALLOCATE( vol_shell2( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell2 in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( mass_shell ))THEN
      ALLOCATE( mass_shell( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( mass_shell2 ))THEN
      ALLOCATE( mass_shell2( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell2 in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( m_parts ))THEN
      ALLOCATE( m_parts( 1:n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    !-------------------------------------------------------!
    !-- Place surfaces based on mass density a that point --!
    !-------------------------------------------------------!

    CALL place_surfaces( central_density, center, radius, m_p, n_shells, &
                         shell_radii, last_r, get_density )

    ! Printout
    PRINT *, " * Number of the spherical surfaces= ", n_shells
    PRINT *, " * Radii of the surfaces in units of the equatorial radius", &
             " of the star, towards the companion= "
    PRINT *, shell_radii/radius
    PRINT *

    !---------------------------------!
    !-- Compute radial mass profile --!
    !---------------------------------!

    PRINT *, " * Integrating the baryon mass density to get the mass profile..."
    PRINT *

    dr             = radius/300.0D0
    dth            = pi/2.0D0/200.0D0
    dphi           = 2.0D0*pi/200.0D0
    CALL integrate_density( center, radius, &
                            central_density, &
                            dr, dth, dphi, &
                            mass, mass_profile, &
                            mass_profile_idx )

    mass_profile( 2:3, : )= mass_profile( 2:3, : )*mass_star/mass

    !---------------------------------------------!
    !-- Assign masses to each spherical surface --!
    !---------------------------------------------!

    CALL assign_surfaces_mass( shell_masses, shell_radii, radius, dr, &
                               n_shells, mass_profile_idx, mass_profile, &
                               mass_star )

    !----------------------------------------------------!
    !-- Print mass profile and surfaces' radii to file --!
    !----------------------------------------------------!

    IF( PRESENT(filename_mass_profile) )THEN
      finalnamefile= filename_mass_profile
    ELSE
      finalnamefile= "mass_profile.dat"
    ENDIF

    IF( PRESENT(filename_shells_radii) )THEN
      finalnamefile2= filename_shells_radii
    ELSE
      finalnamefile2= "shell_radii.dat"
    ENDIF

    CALL print_mass_profile_surface_radii( mass_profile, mass_profile_idx, &
                                           shell_radii, radius, dr, &
                                           n_shells, &
                                           filename_mass_profile, &
                                           filename_shells_radii )

    !---------------------------------------------------------!
    !-- Initialize quantities before starting the iteration --!
    !---------------------------------------------------------!

    PRINT *, " * Initializing quantities before starting the iteration..."
    PRINT *

    ALLOCATE( npart_shell( n_shells ) )
    ALLOCATE( npart_shelleq( n_shells ) )
    ALLOCATE( alpha( n_shells ) )
    ALLOCATE( colatitude_pos( n_shells ) )
    ALLOCATE( pos_shells( n_shells ) )

    initialization: DO r= 1, n_shells, 1

      IF( ALLOCATED( pos_shells( r )% pos_shell ) ) &
        DEALLOCATE( pos_shells( r )% pos_shell )

      IF( ALLOCATED( pos_shells( r )% pvol_shell ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell )

      IF( ALLOCATED( pos_shells( r )% pvol_shell2 ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell2 )

      IF( ALLOCATED( pos_shells( r )% g_xx ) )&
        DEALLOCATE( pos_shells( r )% g_xx )

      IF( ALLOCATED( pos_shells( r )% baryon_density ) ) &
        DEALLOCATE( pos_shells( r )% baryon_density )

      IF( ALLOCATED( pos_shells( r )% gamma_euler ) ) &
        DEALLOCATE( pos_shells( r )% gamma_euler )

      IF( ALLOCATED( pos_shells( r )% pos_th ) ) &
        DEALLOCATE( pos_shells( r )% pos_th )

      IF( ALLOCATED( pos_shells( r )% pos_phi ) ) &
        DEALLOCATE( pos_shells( r )% pos_phi )

      ALLOCATE( pos_shells( r )% pos_shell     ( 3, npart_approx ) )
      ALLOCATE( pos_shells( r )% pvol_shell    (    npart_approx ) )
      ALLOCATE( pos_shells( r )% pvol_shell2   (    npart_approx ) )
      ALLOCATE( pos_shells( r )% g_xx          (    npart_approx ) )
      ALLOCATE( pos_shells( r )% baryon_density(    npart_approx ) )
      ALLOCATE( pos_shells( r )% gamma_euler   (    npart_approx ) )
      ALLOCATE( pos_shells( r )% pos_th        (    npart_approx ) )
      ALLOCATE( pos_shells( r )% pos_phi       (    npart_approx ) )

      pos_shells(r)% pos_shell= 0.0D0
      pos_shells(r)% pos_phi= -1.0D0
      pos_shells(r)% pos_th= -1.0D0
      pos_shells(r)% pvol_shell= 0.0D0
      pos_shells(r)% pvol_shell2= 0.0D0
      pos_shells(r)% g_xx= 0.0D0
      pos_shells(r)% baryon_density= 0.0D0
      pos_shells(r)% gamma_euler= 0.0D0
      m_parts( r )= m_p
      npart_shelleq( r )= CEILING( SQRT(DBLE(2*shell_masses( r )/m_parts( r ))))

    ENDDO initialization

    pos  = 0.0D0
    pmass = 0.0D0
    phase= 0.0D0
    proper_volume= 0.0D0
    vol_shell= 0.0D0
    vol_shell2= 0.0D0
    dr_shells= radius/n_shells
    npart_out= 0
    upper_bound_tmp= upper_bound
    lower_bound_tmp= lower_bound
    r= CEILING(DBLE(n_shells)/2.0D0)
    cnt2= 0
    r_cnt= 1

    ! These array are needed to be able to parallelize the loops on each surface
    ALLOCATE( pos_shell_tmp  ( 3, 5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( g_xx_tmp       (    5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( bar_density_tmp(    5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( gam_euler_tmp  (    5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( pvol_tmp       (    5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( npart_discarded(    5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( npart_surface_tmp(  5*CEILING(SQRT(DBLE(2*npart_approx))), &
                                  5*CEILING(SQRT(DBLE(2*npart_approx))) ) )

    !--------------------------------------------------!
    !--  Main iteration over the spherical surfaces  --!
    !--------------------------------------------------!

    PRINT *, " * Assigning first half of particle positions..."
    PRINT *

    place_particles_on_northern_emispheres: DO

      ! Correct npart_shelleq to be divisible by 4
      IF( MOD( npart_shelleq( r ), 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq( r )= npart_shelleq( r ) + rel_sign
      ENDIF
      IF( MOD( npart_shelleq( r )/2, 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq( r )= 2*( npart_shelleq( r )/2 + rel_sign )
      ENDIF
      IF( MOD( npart_shelleq( r ), 4 ) /= 0 )THEN
        PRINT *, " * ERROR! npart_shelleq(", r, ")=", npart_shelleq( r ), &
                 " is not divisible by 4 in the main iteration. ", &
                 " Check the algorithm. Stopping..."
        STOP
      ENDIF

      ! Compute number of particles on the spherical surface
      npart_shell( r )= ( npart_shelleq( r )**2.0D0 )/2.0D0

      ! Compute angular step in azimuth phi (constant on each shell)
      IF( npart_shelleq( r ) == 0 )THEN
        alpha( r )= 0.0D0
      ELSE
        alpha( r )= 2.0D0*pi/DBLE(npart_shelleq( r ))
      ENDIF

      ! Compute angular positions in colatitude theta,
      ! according to https://mathworld.wolfram.com/SpherePointPicking.html
      IF( ALLOCATED( colatitude_pos( r )% colatitudes ) ) &
        DEALLOCATE( colatitude_pos( r )% colatitudes )
      ALLOCATE( colatitude_pos( r )% colatitudes( npart_shelleq( r )/4 ) )

      IF( shell_radii(r) < 0.95D0*last_r*radius )THEN

        CALL compute_colatitudes_uniformly_in( pi/2.0D0, 9.5D0/10.0D0*pi, &
                                      colatitude_pos( r )% colatitudes( : ) )

      ELSE

        CALL compute_colatitudes_uniformly_in( pi/2.0D0, 9.5D0/10.0D0*pi, &
                                      colatitude_pos( r )% colatitudes( : ) )
        !CALL compute_colatitudes_uniformly_in( pi/2.0D0, 2.0D0/3.0D0*pi, &
        !                              colatitude_pos( r )% colatitudes( : ) )

      ENDIF

          !            alpha( r )*1.0D0/2.0D0 + ( itr2 - 1 )*alpha( r )

        !  ACOS( 2.0D0*( 1.0D0 - COS( pi/3.0D0*( 2.0D0/3.0D0 + DBLE(itr2 - 1)*DBLE(npart_shelleq( r )/4 + 1.0D0 -(1.0D0/2.0D0)-(2.0D0/3.0D0) )/DBLE(npart_shelleq( r )/4 - 1.0D0 ) ) &
        !                   /DBLE(npart_shelleq( r )/4 + 1.0D0 ) ) ) &
        !      - 1.0D0 )
              !5.0D0/12.0D0

        !colatitude_pos( r )% colatitudes( itr2 )= &
        !              colatitude_pos( r )% colatitudes( itr2 ) &
        !              *( 1 + rel_sign*0.05D0*phase_th )
      DO itr2= 1, npart_shelleq( r )/4, 1

        IF( colatitude_pos( r )% colatitudes( itr2 ) <= pi/2.0D0 .OR. &
            colatitude_pos( r )% colatitudes( itr2 ) >= pi &
        )THEN
          PRINT *, "** ERROR! ", &
                   "The colatitudes are not in the OPEN interval (pi/2,pi). ", &
                   "Stopping..."
          STOP
        ENDIF

      ENDDO

      npart_discard    = 0
      npart_shell_cnt  = 0
      npart_shell_tmp  = npart_shell( r )
      ! Initialize te  mporary arrays
      pos_shell_tmp    = huge_real
      g_xx_tmp         = 0.0D0
      bar_density_tmp  = 0.0D0
      gam_euler_tmp    = 0.0D0
      pvol_tmp         = 0.0D0
      npart_discarded  = 0.0D0
      npart_surface_tmp= 0.0D0

      IF( debug ) PRINT *, "Right before OMP, shell ", r, "iteration ", cnt2 + 1

      dphi_shells= alpha(r)

      !$OMP PARALLEL DO DEFAULT(NONE), &
      !$OMP             PRIVATE( phase, col, col_tmp, xtemp, ytemp, ztemp, &
      !$OMP                      dth_shells, delta_r, long, &
      !$OMP                      th, phi, rand_num2, phase_th, rel_sign ), &
      !$OMP             SHARED( r, npart_shelleq, center, rad, alpha, &
      !$OMP                     pos_shells, colatitude_pos, n_shells, &
      !$OMP                     dr_shells, shell_radii, shell_thickness, THIS, &
      !$OMP                     g_xx_tmp, bar_density_tmp, gam_euler_tmp, &
      !$OMP                     pos_shell_tmp, pvol_tmp, dphi_shells, radius, &
      !$OMP                     npart_discarded, npart_surface_tmp, last_r )
      DO phi= 1, npart_shelleq( r ), 1

        IF( debug ) PRINT *, "Right before loop over phi"

        DO th= 1, npart_shelleq( r )/4, 1 !npart_shelleq( r ) is even, see above

          !
          !-- Randomize positions, if specified by the user in the
          !-- parameter file lorene_bns_id_particles.par
          !
          IF( THIS% randomize_phi )THEN

            CALL RANDOM_NUMBER( phase )
            phase= phase*alpha(r)

          ENDIF

        !  IF( shell_radii(r) < 0.95D0*last_r*radius )THEN
        !
        !    long= phase + phi*alpha(r)
        !
        !  ELSE
        !
        !    long= phase + phi*alpha(r)/3.0D0 - pi/3.0D0
        !
        !  ENDIF
          long= phase + phi*alpha(r)


          col= colatitude_pos(r)% colatitudes(th)
          IF( THIS% randomize_theta )THEN

            CALL RANDOM_NUMBER( phase_th )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            col_tmp= col*( 1.0D0 + rel_sign*0.05D0*phase_th )

            IF( col_tmp < pi .AND. col_tmp > pi/2.0D0 )THEN

              col= col_tmp

            ENDIF

          ENDIF

          rad= shell_radii(r)
          IF( THIS% randomize_r )THEN

            CALL RANDOM_NUMBER( delta_r )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            IF( r/n_shells < 0.95D0 )THEN
              rad= rad + rel_sign*delta_r*0.35D0*dr_shells
            ELSE
              !rad= rad - ( 1.0D0 + delta_r )*0.35D0*dr_shells
              rad= rad + ( - delta_r*0.35D0 - 0.5D0 )*dr_shells
            ENDIF

          ENDIF

          IF( rad < 0 )THEN
            PRINT *, " * ERROR! rad < 0. Check the computation of the radial", &
                     " coordinates of the particles. Stopping.."
            STOP
          ENDIF

          !
          !-- Compute Cartesian coordinates of the candidate particle positions
          !
          xtemp= rad*COS(long)*SIN(col) + center
          ytemp= rad*SIN(long)*SIN(col)
          ztemp= rad*COS(col)

          IF( ISNAN( xtemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "xtemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( ISNAN( ytemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ytemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( ISNAN( ztemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ztemp is a NaN. Stopping.."
            STOP
          ENDIF

          ! Import ID needed to compute the particle masses
          CALL get_id( xtemp, ytemp, ztemp, &
                       g_xx_tmp( th, phi ), &
                       bar_density_tmp( th, phi ), &
                       gam_euler_tmp( th, phi ) )

          ! Place a particle at a given position only if the hydro
          ! computed by LORENE is acceptable
          IF( bar_density_tmp( th, phi ) > 0.0D0 &
              !pos_shells(r)% baryon_density( itr + 1 ) > 0.0D0 &
              .AND. &
              validate_position_final( xtemp, ytemp, ztemp ) == 0 )THEN

            !npart_shell_cnt= npart_shell_cnt + 1
            npart_surface_tmp( th, phi )= 1
            pos_shell_tmp( 1, th, phi )= xtemp
            pos_shell_tmp( 2, th, phi )= ytemp
            pos_shell_tmp( 3, th, phi )= ztemp

            ! Compute particle volume
            pvol_tmp( th, phi )= particle_volume( rad, col, dr_shells, &
                                               dth_shells, dphi_shells, th, &
                                               colatitude_pos(r)% colatitudes, &
                                               npart_shelleq(r) )

            ! Safety check
            IF( pvol_tmp( th, phi ) <= 0 )THEN
                ! pos_shells(r)% pvol_shell2( itr + 1 ) <= 0 )THEN
              PRINT *, "When placing first half of particles"
              PRINT *, "pvol_tmp( ", r, ",", th, ",", phi, " ) =", &
                       pvol_tmp( th, phi )
              PRINT *, "dr_shells=", dr_shells, &
                       "dth_shells=", dth_shells, &
                       "dphi_shells=", dphi_shells
              STOP
            ENDIF

          ELSE

            ! If the hydro is not positive, or the position is outside the star,
            ! discard the position and count the number of discarded positions
            !npart_discard= npart_discard + 2
            npart_discarded( th, phi )= 2

          ENDIF

        ENDDO

      ENDDO
      !$OMP END PARALLEL DO

      npart_discard   = SUM( SUM( npart_discarded, DIM= 1 ), DIM= 1 )
      npart_shell_cnt = SUM( SUM( npart_surface_tmp, DIM= 1 ), DIM= 1 )
      npart_shell( r )= npart_shell( r ) - npart_discard
      npart_out       = npart_out + npart_shell( r )/2

      IF( debug ) PRINT *, "Right after OMP"

      ! Safety check
      IF( npart_shell_cnt /= npart_shell( r )/2 )THEN
        PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
        PRINT *, " * npart_shell_cnt=", npart_shell_cnt, &
                 "npart_shell( r )/2=", npart_shell( r )/2
        PRINT *, " * npart_shell_cnt should be equal to npart_shell( r )/2. " &
                 // "Stopping..."
        PRINT *
        STOP
      ENDIF

      ! Set up the next step in pathological cases
      IF( npart_shell( r ) < 0 ) npart_shell( r )= 0
      IF( npart_shell( r ) == 0 )THEN
        m_parts( r )= m_parts( prev_shell )
        PRINT *, " * Placed", npart_shell( r )/2, &
                 " particles on one emisphere of spherical shell ", r, &
                 " out of ", n_shells
        IF( r == 1 )THEN
          EXIT
        ELSEIF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
          !PRINT *, "r=", r
          r= r - 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r == n_shells )THEN
          !PRINT *, "r=", r
          r= CEILING(DBLE(n_shells)/2.0D0) - 1
          r_cnt= r_cnt + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r >= CEILING(DBLE(n_shells)/2.0D0) )THEN
          !PRINT *, "r=", r
          r= r + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ENDIF
      ELSE
        m_parts( r )= shell_masses( r )/DBLE(npart_shell( r ))
      ENDIF

      IF( debug ) PRINT *, " * Before storing the particles"
      IF( debug ) PRINT *, "11"

      IF( debug ) PRINT *, "Right before correction of particle number"
      IF( debug ) PRINT *, "npart_out=", npart_out

      ! If it's not the first populated surface
      not_first_populated_surface: IF( r /= CEILING(DBLE(n_shells)/2.0D0) )THEN

        ! Identify the previous surface
        IF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
          prev_shell= r + 1
        ELSEIF( r > CEILING(DBLE(n_shells)/2.0D0) )THEN
          prev_shell= r - 1
        ELSEIF( r == 1 )THEN
          EXIT
        ENDIF

        ! Logical variables that steer the iteration

        ! This speeds up the iteration considerably, since the inner layers
        ! seem to always want a larger tolerance
        IF( r <= 0.45D0*n_shells )THEN
          upper_bound_tmp= upper_bound_tmp*1.1D0
          lower_bound_tmp= lower_bound_tmp*0.9D0
        ENDIF

        ! Is the particle mass too high?
        high_mass= m_parts( r )/m_parts( prev_shell ) > upper_bound_tmp
        ! Is the particle mass too low?
        low_mass = m_parts( r )/m_parts( prev_shell ) < lower_bound_tmp
        ! How many positions were kept, placing a particle on them?
        npart_shell_kept= DBLE(npart_shell( r ))/DBLE(npart_shell_tmp)
        ! Were all the positions kept?
        kept_all = npart_shell_kept == 1.0D0

        ! If the particle mass is too high and all positions were kept
        adjust_particle_number_surface: IF( high_mass .AND. kept_all )THEN

          IF( debug ) PRINT *, "Case 1"

          cnt2= cnt2 + 1

          ! If this is the (max_steps + 1)st step
          too_many_steps: IF( cnt2 > max_steps )THEN

            ! Allow for a bit more different particle mass
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor

            !
            !-- Special treatment for the positions near the surface
            !

            ! If the range of particle masses is getting too generous
            ! near the surface
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) > 1.1D0*upper_bound &
            )THEN

              ! Increase the number of positions on this surface
              ! by a factor between 5 and 10
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= rand_num*npart_shelleq( r )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )

              ! Reset the particle mass tolerance range
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound

              ! Reset total particle number
              npart_out= npart_out - npart_shell( r )/2

              ! Reset counter
              cnt2= 1

              ! Replace particles on this surface
              CYCLE

            ENDIF

            ! Reset counter
            cnt2= 1

          ENDIF too_many_steps

          ! If this is not yet the (max_steps + 1)st step

          ! Reset total particle number
          npart_out= npart_out - npart_shell( r )/2

          ! Increase randomly particle number on the equator which determines
          ! the particle number on the surface
          ! More particles = lower particle mass, with fixed surface mass
          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) + 1*NINT( 1 + 1.0*rand_num ) &
                                                 + 1*NINT( 1 + 1.0*rand_num2 )

          ! Treat pathological cases
          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          ! Replace particles on this surace
          CYCLE

        ! The cases below do similar things to the one above, so there are no
        ! comments on each line. See the case above for explanations
        ELSEIF( low_mass .AND. kept_all )THEN

          IF( debug ) PRINT *, "Case 2"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) < 0.9D0*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= npart_shelleq( r )/rand_num
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) - 1*NINT( 1 + 1.0*rand_num ) &
                                                 - 1*NINT( 1 + 1.0*rand_num2 )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( high_mass .AND. .NOT.kept_all ) THEN

          IF( debug ) PRINT *, "Case 3"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) > 1.1D0*upper_bound &
                !upper_bound_tmp > 1.1D0*upper_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= rand_num*npart_shelleq( r )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1

          ! If x% of the positions were kept, divide the old particle number
          ! on the equator by x, and adjust with some other random and
          ! non-random factors which turn out to work well a posteriori
          npart_shelleq( r )= CEILING( SQRT( &
                                2*(shell_masses( r )/m_parts( prev_shell )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( prev_shell ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( low_mass .AND. .NOT.kept_all ) THEN

          IF( debug ) PRINT *, "Case 4"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) < 0.9D0*lower_bound &
                !lower_bound_tmp < 0.9D0*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= npart_shelleq( r )/rand_num
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1

          ! If x% of the positions were kept, divide the old particle number
          ! on the equator by x, and adjust with some other random and
          ! non-random factors which turn out to work well a posteriori
          npart_shelleq( r )= CEILING( SQRT( &
                                2*(shell_masses( r )/m_parts( prev_shell )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( prev_shell ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          IF( debug ) PRINT *, "Right after correction of particle number"

          ! Replace particle on this surface
          CYCLE

        ENDIF adjust_particle_number_surface

      ENDIF not_first_populated_surface

      ! >TODO: Safety check to be updated
      !IF( r_cnt > 1 )THEN
      !  npart_test= 0
      !  DO itr= 1, r, 1
      !    npart_test= npart_test + npart_shell( itr )
      !  ENDDO
      !  IF( npart_test/2 /= npart_out )THEN
      !    PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
      !             "equal to the total number of particles. Stopping.."
      !    PRINT *, "npart_test=", npart_test/2, ", npart_out=", npart_out
      !    PRINT *
      !    STOP
      !  ENDIF
      !ENDIF

      IF( debug ) PRINT *, "10"

      ! At this point, the particles are placed on this surface
      ! Print out the result
      PRINT *, " * Placed", npart_shell( r )/2, &
               " particles on one emisphere of spherical shell ", r, &
               " out of ", n_shells
      PRINT *, "   Shell radius= ", shell_radii( r )/radius*100.0D0, &
              "% of the radius of the star"
      PRINT *, "   Placed", npart_out, " particles overall, so far."
      IF( r /= CEILING(DBLE(n_shells)/2.0D0) ) PRINT *, &
               "   Ratio of particle masses on last 2 shells: ", &
               "   m_parts(", r, ")/m_parts(", prev_shell, ")= ",  &
               m_parts( r )/m_parts( prev_shell )

      ! Save particles to non-temporary variables
      size_pos_shell= SIZE( pos_shells(r)% pos_shell( 1, : ) )

      IF( npart_shelleq(r)*npart_shelleq(r)/4 > size_pos_shell )THEN

        CALL reallocate_array_2d( pos_shells(r)% pos_shell, 3, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_shells(r)% g_xx, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_shells(r)% baryon_density, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_shells(r)% gamma_euler, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_shells(r)% pvol_shell2, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

      ENDIF

      itr= 0
      DO th= 1, npart_shelleq(r)/4, 1
        DO phi= 1, npart_shelleq(r), 1

          IF( pos_shell_tmp( 1, th, phi ) /= huge_real &
              !pos_shell_tmp( 1, th, phi ) < center + 1.2D0*radius &
              !.AND. &
              !pos_shell_tmp( 1, th, phi ) > center - 1.2D0*radius
          )THEN

            itr= itr + 1

            pos_shells(r)% pos_shell( 1, itr )= pos_shell_tmp( 1, th, phi )
            pos_shells(r)% pos_shell( 2, itr )= pos_shell_tmp( 2, th, phi )
            pos_shells(r)% pos_shell( 3, itr )= pos_shell_tmp( 3, th, phi )

            pos_shells(r)% g_xx( itr )          = g_xx_tmp( th, phi )
            pos_shells(r)% baryon_density( itr )= bar_density_tmp( th, phi )
            pos_shells(r)% gamma_euler( itr )   = gam_euler_tmp( th, phi )
            pos_shells(r)% pvol_shell2( itr )   = pvol_tmp( th, phi )


          ENDIF

        ENDDO
      ENDDO
      ! Safety check
      IF( npart_shell_cnt /= itr )THEN
        PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
        PRINT *, " * npart_shell_cnt=", npart_shell_cnt, &
                 ", itr=", itr, ". npart_shell_cnt should be equal to itr. "
        PRINT *, " * npart_shell( r )/2=", npart_shell( r )/2
        STOP
      ENDIF

      ! Set up next step
      IF( r == n_shells )THEN
        r= CEILING(DBLE(n_shells)/2.0D0) - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "last shell"
      ELSEIF( r == 1 )THEN
        IF( debug ) PRINT *, "exit"
        EXIT
      ELSEIF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
        r= r - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "inner layers"
      ELSEIF( r >= CEILING(DBLE(n_shells)/2.0D0) )THEN
        r= r + 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "outer layers"
      ENDIF

      IF( debug ) PRINT *, "12"

    ENDDO place_particles_on_northern_emispheres

    !-----------------------------!
    !--  End of main iteration  --!
    !-----------------------------!

    ! Print out the total number of particles on the northern emispheres,
    ! and the final mass ratio
    PRINT *, " * Particles on the northern emispheres=", npart_out
    PRINT *, " * Particle mass ratio= ", MAXVAL(m_parts)/MINVAL(m_parts)
    PRINT *

    ! Safety check
    npart_test= SUM( npart_shell, DIM= 1 )
    IF( npart_test/2 /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * npart_test/2=", npart_test/2, ", npart_out=", npart_out
      PRINT *, " * Array npart_shell=", npart_shell
      PRINT *
      STOP
    ENDIF

    IF( debug ) PRINT *, "13"

    ! Deallocate temporary arrays
    IF( ALLOCATED(pos_shell_tmp) )THEN
      DEALLOCATE( pos_shell_tmp  , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(g_xx_tmp) )THEN
      DEALLOCATE( g_xx_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(bar_density_tmp) )THEN
      DEALLOCATE( bar_density_tmp, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(gam_euler_tmp) )THEN
      DEALLOCATE( gam_euler_tmp  , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(pvol_tmp) )THEN
      DEALLOCATE( pvol_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(npart_discarded) )THEN
      DEALLOCATE( npart_discarded       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array npart_discarded in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(npart_surface_tmp) )THEN
      DEALLOCATE( npart_surface_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "14"

    !
    !-- Mirror particles from the northern emispheres to the southern ones
    !
    PRINT *, " * Mirroring particles from the northern emispheres to the", &
             " southern ones..."
    PRINT *

    IF( debug ) PRINT *, " * npart/2=", npart_out

    DO r= 1, n_shells, 1

      DO itr= 1, npart_shell( r )/2, 1

        npart_out= npart_out + 1

        pos_shells(r)% pos_shell( 1, npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% pos_shell( 1, itr )
        pos_shells(r)% pos_shell( 2, npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% pos_shell( 2, itr )
        pos_shells(r)% pos_shell( 3, npart_shell( r )/2 + itr )= &
                                        - pos_shells(r)% pos_shell( 3, itr )
        pos_shells(r)% g_xx( npart_shell( r )/2 + itr )= &
                                                  pos_shells(r)% g_xx( itr )
        pos_shells(r)% baryon_density( npart_shell( r )/2 + itr )= &
                                        pos_shells(r)% baryon_density( itr )
        pos_shells(r)% gamma_euler( npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% gamma_euler( itr )
        pos_shells(r)% pvol_shell2( npart_shell( r )/2 + itr )= &
                                           pos_shells(r)% pvol_shell2( itr )

        ! Safety checks
        IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, r, itr, pos_shells(r)% pos_shell( 1, itr ), &
                   pos_shells(r)% pos_shell( 2, itr ), &
                   pos_shells(r)% pos_shell( 3, itr ), &
                   pos_shells(r)% baryon_density( itr )
          STOP
        ENDIF
        IF( pos_shells(r)% pvol_shell2( itr ) < 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
                   pos_shells(r)% pvol_shell2( itr )
          STOP
        ENDIF

      ENDDO

    ENDDO
    PRINT *, " * Final number of particles=", npart_out
    PRINT *

    ! Safety checks (maybe redundant at this point, but better to be paranoid)
    npart_test= 0
    DO r= 1, n_shells, 1
      npart_test= npart_test + npart_shell( r )
    ENDDO
    IF( npart_test /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * npart_test", npart_test, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF
    IF( SUM( npart_shell, DIM=1 ) /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * SUM( npart_shell )", SUM( npart_shell, DIM=1 ), &
               ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

    debug_if: IF( debug )THEN

      mass_test= 0.0D0
      mass_test2= 0.0D0
      proper_volume_test= 0.0D0
      proper_volume= 0.0D0
      i_shell= 1
      DO r= 1, n_shells, 1
        !DO itr= i_shell, (i_shell - 1) + (npart_shelleq(r)**2.0D0)/2.0D0, 1
        !  CALL bns_obj% import_id( &
        !           pos( 1, itr ), pos( 2, itr ), pos( 3, itr ), &
        !           g_xx, baryon_density, gamma_euler )
        !
        !  pvol( itr )= m_parts( r )/( baryon_density*g_xx*SQRT(g_xx)*gamma_euler )
        !
        !  proper_volume_test= proper_volume_test + 2.0D0*pvol( itr )*g_xx*SQRT(g_xx)
        !  mass_test= mass_test + &
        !              2.0D0*baryon_density*pvol( itr )*g_xx*SQRT(g_xx)*gamma_euler
        !ENDDO
        !i_shell= i_shell + (npart_shelleq(r)**2.0D0)/2.0D0
        DO itr= 1, npart_shell( r ), 1

          IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
            PRINT *, "When computing particle volume"
            PRINT *, r, itr, pos_shells(r)% pos_shell( 1, itr ), &
                     pos_shells(r)% pos_shell( 2, itr ), &
                     pos_shells(r)% pos_shell( 3, itr ), &
                     pos_shells(r)% baryon_density( itr )
          ENDIF
          !IF( pos_shells(r)% pvol_shell2( itr ) <= 0.0D0 )THEN
          !  PRINT *, "When computing particle volume"
          !  PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
          !           pos_shells(r)% pvol_shell2( itr )
          !  STOP
          !ENDIF

          pos_shells(r)% pvol_shell( itr )= m_parts( r ) &
                            /( pos_shells(r)% baryon_density( itr ) &
                              *pos_shells(r)% g_xx( itr ) &
                              *SQRT(pos_shells(r)% g_xx( itr )) &
                              *pos_shells(r)% gamma_euler( itr ) )

          proper_volume_test= proper_volume_test + &
                              pos_shells(r)% pvol_shell( itr )! &
                              !*pos_shells(r)% g_xx( itr ) &
                              !*SQRT(pos_shells(r)% g_xx( itr ))

          proper_volume= proper_volume + pos_shells(r)% pvol_shell2( itr )
          !PRINT *, proper_volume

          mass_test= mass_test + pos_shells(r)% baryon_density( itr ) &
                    *pos_shells(r)% pvol_shell( itr ) &
                    *pos_shells(r)% g_xx( itr ) &
                    *SQRT(pos_shells(r)% g_xx( itr )) &
                    *pos_shells(r)% gamma_euler( itr )
          mass_test2= mass_test2 + pos_shells(r)% baryon_density( itr ) &
                    *pos_shells(r)% pvol_shell2( itr ) &
                    *pos_shells(r)% g_xx( itr ) &
                    *SQRT(pos_shells(r)% g_xx( itr )) &
                    *pos_shells(r)% gamma_euler( itr )

        ENDDO
      ENDDO

      PRINT *, mass_test, mass_test2, mass_star, proper_volume_test, proper_volume
      PRINT *

      mass_test = 0.0D0
      mass_test2= 0.0D0
      proper_volume_test = 0.0D0
      proper_volume= 0.0D0
      vol_shell( r )  = 0.0D0
      vol_shell2( r ) = 0.0D0
      mass_shell( r ) = 0.0D0
      mass_shell2( r )= 0.0D0
      DO r= 1, n_shells, 1
        DO itr= 1, npart_shell( r ), 1
          !IF( pos_shells(r)% pvol_shell2( itr ) <= 0 )THEN
          !  PRINT *, "When computing shell volumes and masses"
          !  PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
          !           pos_shells(r)% pvol_shell2( itr )
          !  STOP
          !ENDIF
          vol_shell( r )  = vol_shell( r )  + pos_shells(r)% pvol_shell( itr )
          vol_shell2( r ) = vol_shell2( r ) + pos_shells(r)% pvol_shell2( itr )
          mass_shell( r ) = mass_shell( r ) + &
                            pos_shells(r)% baryon_density( itr ) &
                           *pos_shells(r)% pvol_shell( itr ) &
                           *pos_shells(r)% g_xx( itr ) &
                           *SQRT(pos_shells(r)% g_xx( itr )) &
                           *pos_shells(r)% gamma_euler( itr )
          mass_shell2( r )= mass_shell2( r ) + &
                            pos_shells(r)% baryon_density( itr ) &
                           *pos_shells(r)% pvol_shell2( itr ) &
                           *pos_shells(r)% g_xx( itr ) &
                           *SQRT(pos_shells(r)% g_xx( itr )) &
                           *pos_shells(r)% gamma_euler( itr )
        ENDDO
        mass_test= mass_test + mass_shell( r )
        mass_test2= mass_test2 + mass_shell2( r )
        proper_volume= proper_volume + vol_shell( r )
        proper_volume_test= proper_volume_test + vol_shell2( r )
        IF( r > 1 )THEN
          PRINT *, "shell", r
          PRINT *, "  shell volumes:", vol_shell( r ), vol_shell2( r ), &
                   4.0D0/3.0D0*pi* &
                   ( shell_radii( r )**3.0D0 - shell_radii( r - 1 )**3.0D0 )
          PRINT *, "  shell masses:", mass_shell( r ), mass_shell2( r ), &
                   shell_masses( r )
          PRINT *
        ENDIF
      ENDDO
      PRINT *
      PRINT *, "masses of the star:", mass_test, mass_test2, mass_star
      PRINT *, "volumes of the star:", proper_volume, proper_volume_test, &
                                       4.0D0/3.0D0*pi*radius**3.0D0
      PRINT *

      !STOP

      !DO r= 1, n_shells, 1
      !  DO itr= 1, npart_shell( r ), 1
      !
      !    PRINT*, (m_parts( r )*MSun/amu)/pos_shells(r)% pvol_shell( itr ) &
      !            /(pos_shells(r)% g_xx( itr ) &
      !              *SQRT(pos_shells(r)% g_xx( itr )) &
      !              *pos_shells(r)% gamma_euler( itr )), &
      !            pos_shells(r)% baryon_density( itr )*MSun/amu
      !  ENDDO
      !ENDDO
      !STOP

    ENDIF debug_if

    !---------------------------------------!
    !--  Save particles to output arrays  --!
    !---------------------------------------!

    IF(.NOT.ALLOCATED( pos ))THEN
      ALLOCATE( pos( 3, npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( pvol ))THEN
      ALLOCATE( pvol( npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( pmass ))THEN
      ALLOCATE( pmass( npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    cnt= 0
    npart_shell_tmp= 0
    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r ), 1
        pos( 1, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 1, itr )
        pos( 2, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 2, itr )
        pos( 3, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 3, itr )
        pvol( itr + npart_shell_tmp )  = pos_shells(r)% pvol_shell2( itr )
        pmass( itr + npart_shell_tmp ) = m_parts( r )
        cnt= cnt + 1
      ENDDO
      npart_shell_tmp= cnt
    ENDDO
    ! Safety check
    IF( cnt /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "cnt", cnt, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

    !-------------------------------------------------------------------!
    !--  Print particle positions to file (TODO: make this optional)  --!
    !-------------------------------------------------------------------!

    PRINT *, " * Printing particle positions to file..."
    PRINT *

    IF( PRESENT(filename_shells_pos) )THEN
      finalnamefile= filename_shells_pos
    ELSE
      finalnamefile= "spherical_surfaces_pos.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
              ". The error message is", err_msg
      STOP
    ENDIF

    !DO itr = 1, npart_out, 1
    !
    !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    !    pos( 1, itr ), pos( 2, itr ), pos( 3, itr )
    !
    !  IF( ios > 0 )THEN
    !    PRINT *, "...error when writing the arrays in " &
    !             // TRIM(finalnamefile), ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !
    !ENDDO

    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r ), 1

        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          r, pos_shells(r)% pos_shell( 1, itr ), &
          pos_shells(r)% pos_shell( 2, itr ), &
          pos_shells(r)% pos_shell( 3, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in " &
                   // TRIM(finalnamefile), ". The error message is", err_msg
          STOP
        ENDIF

      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

    PRINT *, " * SUBROUTINE place_particles_spherical_surfaces executed."
    PRINT *

    IF( debug ) PRINT *, "20"



    CONTAINS



    FUNCTION validate_position_final( x, y, z ) RESULT( answer )

      !*******************************************************
      !
      !# Returns validate_position( x, y, z ) if the latter
      !  is present, 0 otherwise
      !
      !  FT 22.09.2021
      !
      !*******************************************************

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: z
      !! \(z\) coordinate of the desired point
      INTEGER:: answer
      !! validate_position( x, y, z ) if the latter is present, 0 otherwise

      IF( PRESENT(validate_position) )THEN

        answer= validate_position( x, y, z )

      ELSE

        answer= 0

      ENDIF

    END FUNCTION validate_position_final

  END PROCEDURE place_particles_spherical_surfaces


  FUNCTION number_surfaces( m_p, center, radius, get_dens ) &
           RESULT( n_shells_tmp )

    !************************************************
    !
    !# Compute the number of spherical surfaces
    !  by integrating the linear particle density
    !  along the larger equatorial radius
    !
    !  FT 22.07.2021
    !
    !************************************************

    USE constants, ONLY: third

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT( IN ):: m_p, center, radius
    INTERFACE
      FUNCTION get_dens( x, y, z ) RESULT( density )
        !! Returns the baryon mass density at the desired point
        DOUBLE PRECISION, INTENT(IN):: x
        !! \(x\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: y
        !! \(y\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: z
        !! \(z\) coordinate of the desired point
        DOUBLE PRECISION:: density
        !> Baryon mass density at \((x,y,z)\)
      END FUNCTION get_dens
    END INTERFACE


    DOUBLE PRECISION:: n_shells_tmp

    INTEGER:: r
  !  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: particle_profile
  !
  !  IF(.NOT.ALLOCATED( particle_profile ))THEN
  !    ALLOCATE( particle_profile( 2, 500 ), STAT= ios, &
  !              ERRMSG= err_msg )
  !    IF( ios > 0 )THEN
  !       PRINT *, "...allocation error for array particle_profile in" &
  !                // "FUNCTION number_surfaces. ", &
  !                "The error message is", err_msg
  !       STOP
  !    ENDIF
  !  ENDIF

    n_shells_tmp= 0.0D0
  !  particle_profile= 0.0D0

    DO r= 1, 500, 1

      n_shells_tmp= n_shells_tmp + &
                      radius/500*( ( get_dens( &
                                     center + r*radius/500, 0.0D0, 0.0D0 ) &
                                     )/m_p )**third
      !particle_profile( 1, r )= r*radius/500
      !particle_profile( 2, r )= n_shells_tmp

    ENDDO

    n_shells_tmp= NINT( n_shells_tmp )

  END FUNCTION number_surfaces


  SUBROUTINE reallocate_array_1d( array, new_dim )

    !************************************
    !
    !# Reallocate a 1-dimensional array
    !
    !  FT 22.07.2021
    !
    !************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(IN OUT):: array
    INTEGER, INTENT(IN):: new_dim

    DEALLOCATE( array, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...deallocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

    ALLOCATE( array( new_dim ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

  END SUBROUTINE reallocate_array_1d


  SUBROUTINE reallocate_array_2d( array, new_dim, new_dim2 )

    !************************************
    !
    !# Reallocate a 2-dimensional array
    !
    !  FT 22.07.2021
    !
    !************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(IN OUT):: array
    INTEGER, INTENT(IN):: new_dim, new_dim2

    DEALLOCATE( array, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...deallocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

    ALLOCATE( array( new_dim, new_dim2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

  END SUBROUTINE reallocate_array_2d


  SUBROUTINE place_surfaces( central_dens, center, radius, m_p, n_shells, &
                             shell_radii, last_r, get_dens )

    !************************************************
    !
    !# Place the spherical surface, according to
    !  the baryon mass density of the star
    !  along the larger equatorial radius
    !
    !  FT 23.07.2021
    !
    !************************************************

    USE constants,  ONLY: third

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_shells
    DOUBLE PRECISION, INTENT( IN ):: central_dens, center, radius, m_p, last_r
    DOUBLE PRECISION, DIMENSION( n_shells ), INTENT( IN OUT ):: shell_radii
    INTERFACE
      FUNCTION get_dens( x, y, z ) RESULT( density )
        !! Returns the baryon mass density at the desired point
        DOUBLE PRECISION, INTENT(IN):: x
        !! \(x\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: y
        !! \(y\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: z
        !! \(z\) coordinate of the desired point
        DOUBLE PRECISION:: density
        !> Baryon mass density at \((x,y,z)\)
      END FUNCTION get_dens
    END INTERFACE

    DOUBLE PRECISION:: rho_tmp

    !central_density= bns_obj% get_rho_center1()

    shell_radii= 0.0D0

    shell_radii(1)= ( central_dens/m_p )**(-third)
    DO itr= 2, n_shells, 1

      rho_tmp= get_dens( center + shell_radii( itr - 1 ), &
                                             0.0D0, 0.0D0 )

      IF( rho_tmp == 0 )THEN
        shell_radii= shell_radii*itr/n_shells
      ENDIF

      shell_radii( itr )= shell_radii( itr - 1 ) + ( rho_tmp/m_p )**(-third)

    ENDDO
    shell_radii= shell_radii*(radius*last_r/shell_radii(n_shells))

  END SUBROUTINE place_surfaces


  SUBROUTINE assign_surfaces_mass( shell_masses, shell_radii, radius, dr, &
                                   n_shells, mass_profile_idx, mass_profile, &
                                   mass_star )

    !*************************************************
    !
    !# Assign a mass to each spherical surface,
    !  based on the radial mass profile of the star
    !  (computed along the larger equatorial radius)
    !
    !  FT 23.07.2021
    !
    !*************************************************

    USE constants,  ONLY: third

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_shells
    DOUBLE PRECISION, INTENT( IN ):: radius, dr, mass_star

    INTEGER, DIMENSION( : ),                 INTENT( IN ):: mass_profile_idx
    DOUBLE PRECISION, DIMENSION( n_shells ), INTENT( IN ):: shell_radii
    DOUBLE PRECISION, DIMENSION( :, : ),     INTENT( IN ):: mass_profile
    DOUBLE PRECISION, DIMENSION( n_shells ), INTENT( IN OUT ):: shell_masses

    INTEGER shell_index, itr2

    shell_index= 1
    itr2= 1
    shell_masses= 0.0D0
    assign_masses_to_surfaces: DO itr= 1, NINT(radius/dr), 1

      IF( shell_index == n_shells )THEN

        shell_masses( shell_index )= SUM( mass_profile( 2, &
         mass_profile_idx(itr2):mass_profile_idx(NINT(radius/dr)-1) ), DIM= 1 )

        EXIT

      ENDIF

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= shell_radii( shell_index ) &!+ radius/DBLE(2*n_shells)
      )THEN

       shell_masses( shell_index )= SUM( mass_profile( 2, &
                     mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

       itr2= itr + 1
       shell_index= shell_index + 1

      ENDIF

    ENDDO assign_masses_to_surfaces

    ! Safety check
    IF( ABS( SUM( shell_masses, DIM= 1 ) - mass_star )/mass_star > 5.0D-3 )THEN
      PRINT *, " ** The masses of the shells do not add up to the ", &
               "mass of the star. Stopping..."
      PRINT *, " * SUM( shell_masses )= ", SUM( shell_masses, DIM=1 )
      PRINT *, " * Baryon mass of the star= ", mass_star
      PRINT *, " * Array shell_masses=", shell_masses
      PRINT *
      STOP
    ENDIF

  END SUBROUTINE assign_surfaces_mass


  SUBROUTINE print_mass_profile_surface_radii( mass_profile, mass_profile_idx, &
                                               shell_radii, radius, dr, &
                                               n_shells, &
                                               filename_mass_profile, &
                                               filename_shells_radii )

    !*************************************************
    !
    !# Print star's radial mass profile and radii of
    !  spherical surfaces to different ASCII files
    !
    !  FT 23.07.2021
    !
    !*************************************************

    USE constants,  ONLY: third

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_shells
    DOUBLE PRECISION, INTENT( IN ):: radius, dr

    INTEGER, DIMENSION( : ),                 INTENT( IN ):: mass_profile_idx
    DOUBLE PRECISION, DIMENSION( n_shells ), INTENT( IN ):: shell_radii
    DOUBLE PRECISION, DIMENSION( :, : ),     INTENT( IN ):: mass_profile

    CHARACTER( LEN= * ), INTENT( IN ):: filename_mass_profile, &
                                        filename_shells_radii

    LOGICAL:: exist

    PRINT *, " * Print mass profile to file..."
    PRINT *

    INQUIRE( FILE= TRIM(filename_mass_profile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(filename_mass_profile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 2, FILE= TRIM(filename_mass_profile), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(filename_mass_profile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    write_data_loop: DO itr = 1, NINT(radius/dr) - 1, 1

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        mass_profile( 1, mass_profile_idx(itr) ), &
        mass_profile( 2, mass_profile_idx(itr) ), &
        mass_profile( 3, mass_profile_idx(itr) )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(filename_mass_profile), ". The error message is", &
                 err_msg
        STOP
      ENDIF

    ENDDO write_data_loop

    CLOSE( UNIT= 2 )

    PRINT *, " * Print surfaces' radii to file..."
    PRINT *

    INQUIRE( FILE= TRIM(filename_shells_radii), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(filename_shells_radii), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 2, FILE= TRIM(filename_shells_radii), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(filename_shells_radii), &
              ". The error message is", err_msg
      STOP
    ENDIF

    DO itr = 1, n_shells, 1

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        shell_radii( itr )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(filename_shells_radii), ". The error message is", &
                 err_msg
        STOP
      ENDIF

    ENDDO

    CLOSE( UNIT= 2 )

  END SUBROUTINE print_mass_profile_surface_radii


  FUNCTION particle_volume( rad, col, dr_shells, dth_shells, dphi_shells, th, &
                            colatitudes, npart_equator ) RESULT( pvol )

    !*******************************************
    !
    !# Compute the geometrical particle volume
    !  not the proper particle volume.
    !
    !  FT 23.07.2021
    !
    !*******************************************

    USE constants,  ONLY: pi

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: th, npart_equator
    DOUBLE PRECISION, INTENT( IN ):: rad, col, dr_shells, dphi_shells
    DOUBLE PRECISION, INTENT( IN OUT ):: dth_shells
    DOUBLE PRECISION, DIMENSION(:), INTENT( IN ):: colatitudes

    DOUBLE PRECISION:: pvol

    IF( th == 1 )THEN

    !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
      IF( npart_equator == 4 )THEN

        dth_shells= pi

      ELSE

        dth_shells= 2.0D0*ABS( col - &
                ( col + colatitudes(th + 1) )/2.0D0 )
      ENDIF

    ELSEIF( th == npart_equator/4 )THEN

    !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
      dth_shells= 2.0D0*ABS( ( colatitudes(th - 1) &
                + col )/2.0D0 - col )

    ELSE

      dth_shells= ABS( &
              ( colatitudes(th + 1) + col )/2.0D0 &
            - ( col + colatitudes(th - 1) )/2.0D0 )

    ENDIF

    pvol= rad**2.0D0*SIN(col)*dr_shells*dth_shells*dphi_shells! &

  END FUNCTION particle_volume


  SUBROUTINE compute_colatitudes_uniformly_in( alpha, beta, colatitudes )

    !**************************************************
    !
    !# Compute the colatitudes according to a
    !  uniform distribution over a spherical
    !  surface, between alpha and beta, with
    !  pi/2 < alpha < beta < pi.
    !  The values are stored in the array colatitudes
    !  See https://mathworld.wolfram.com/SpherePointPicking.html
    !
    !  FT 6.10.2021
    !
    !**************************************************

    USE constants, ONLY: pi

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: alpha, beta
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: colatitudes

    INTEGER:: n, size_col, i

    IF( alpha < pi/2 .OR. alpha > pi )THEN
      PRINT *, "ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument alpha should lie in [pi/2,pi]. Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( beta < pi/2 .OR. beta > pi )THEN
      PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument beta should lie in [pi/2,pi]. Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( alpha > beta )THEN
      PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument alpha should be less than argument beta. Stopping..."
      PRINT *
      STOP
    ENDIF

    size_col= SIZE(colatitudes)
    n= 4*size_col

    DO i= 1, size_col, 1

      colatitudes(i)= ACOS( 2.0D0*( DBLE( i + 1 )*( COS(alpha) + 1.0D0 )/2.0D0 &
            + ( ( COS(beta) + 1.0D0 )/2.0D0 &
              - (DBLE(n + 1)/4.0D0 + 1.0D0)*( COS(alpha) + 1.0D0 )/2.0D0 ) &
                      *4.0D0*DBLE(i)/DBLE(n + 1) ) - 1.0D0 )

      !PRINT *, "colatitudes(", i, ")", colatitudes(i)/pi, "pi"

      IF( ISNAN( colatitudes(i) ) )THEN
        PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in! ", &
                 "colatitudes(", i, ") is a NaN! Stopping.."
        PRINT *, DBLE( i + 1 )*( COS(alpha) + 1.0D0 )/2.0D0
        PRINT *, ( COS(beta) + 1.0D0 )/2.0D0
        PRINT *, DBLE(size_col) + 1.0D0
        PRINT *, DBLE(i)/DBLE(size_col)
        PRINT *
        STOP
      ENDIF

    ENDDO

  END SUBROUTINE


END SUBMODULE spherical_surfaces
