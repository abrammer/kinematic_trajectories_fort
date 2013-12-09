!###############################################################################
!  calctraj.wpet.lint_cubep.f90
!  Author: Matthew Janiga (matthew.janiga@gmail.com)
!  Last Updated: Mar. 25, 2011
!
!  Description: Calculation program using linear time and cubic pressure.
!  All calctraj programs use bicubic x and y.
!###############################################################################

MODULE mod_calctraj
  USE mod_fileio
  USE mod_grid
  USE netcdf
  USE mod_traj
  USE mod_checkinp

  IMPLICIT NONE
  SAVE

  ! #### Derived types ##### 

  ! Contains the 3D wind.
  TYPE vec_wind
     REAL :: u,v,w,sp
  END TYPE vec_wind


  ! Contains the position as described by the descriptor arrays (time,lev,lat,lon).
  ! The position in geographic space.
  TYPE geo_pos
     REAL :: time,lev,lat,lon
  END TYPE geo_pos


  ! Interpolation requries two sets of values the interpolation weighting (iw*)
  ! and the neighboorhood to read in (nb*). In addition, real pressure values
  ! are stored for the vertical coordinate to weight closer values more in the
  ! vertical polynomial interpolation.

  ! i, j, k, and l contain the full fractional element position (i.e. 1:nlat)
  ! iwi, iwj, iwk, and iwl contain the fractional position in the neigborhood (-1:2)
  ! p is the current pressure (equivalent to geo_pos%lev)
  ! ploc is the local array of pressures in the neighboorhood

  TYPE interp_coord
     REAL :: i, j, k, l
     REAL :: iwi, iwj, iwk, iwl
     INTEGER, DIMENSION(4) :: nbj, nbk, nbl
     INTEGER, DIMENSION(2) :: nbi
     REAL :: p
     REAL, DIMENSION(4) :: ploc
  END TYPE interp_coord


CONTAINS

  !###########################################################################
  ! calctraj()                                                            ####
  ! Calls individual subroutines which perform calculate the trajectories ####
  !###########################################################################

  SUBROUTINE calctraj()
    IMPLICIT NONE

    ! Local variables
    INTEGER :: i,t          ! Time step and trajectory counters


    ! Open the grid file and auxiliary file.
    CALL check( nf90_open( gridfile, nf90_nowrite, gncid) )
    CALL check( nf90_open( auxfile, nf90_nowrite, ancid) )

    ! Initializes the full trajectory arrays with initial positions
    c_step = 1
    CALL init_traj() 

    ! Begin loop of time steps
    DO i=2,nstep
       c_step = i
       IF( ANY(t_on) )THEN
          PRINT *, 'Time step: ',i
          DO t=1,ntraj 
             IF( t_on(t) )THEN
                CALL petterson(t)
             END IF
          END DO
       ELSE
          PRINT *, 'All trajectories have left the domain, time steping halted'
          EXIT
       END IF
    END DO

    ! Close the grid file.
    CALL check( nf90_close(gncid) )
    CALL check( nf90_close(ancid) )

  END SUBROUTINE calctraj


  !###########################################################################
  ! init_traj()                                                           ####
  ! Initializes the trajectory arrays finds initial values of wind / aux  ####
  !###########################################################################

  SUBROUTINE init_traj()
    IMPLICIT NONE

    ! Local variables
    INTEGER :: t,a                       ! Counters
    TYPE(interp_coord) :: ic             ! Interpolation coordinates
    TYPE(vec_wind) :: wind               ! 3D wind vector
    TYPE(geo_pos) :: gpos                ! Geographic position
    LOGICAL :: in_grid,no_miss,no_grnd   ! Success flags


    ! Initialize arrays to missing
    t_time = mv
    t_lev = mv
    t_lat = mv
    t_lon = mv
    t_u = mv
    t_v = mv
    t_w = mv
    t_sp = mv

    ! Initialize last valid to last total
    last_valid = nstep

    PRINT *, "Trajectory arrays initialized to missing."

    ! Move over initial coordinates
    DO t=1,ntraj
       IF( t_on(t) )THEN

          in_grid = .TRUE.
          no_miss = .TRUE.
          no_grnd = .TRUE.

          t_time(t,c_step) = i_time(t)
          t_lev(t,c_step) = i_lev(t)
          t_lat(t,c_step) = i_lat(t)
          t_lon(t,c_step) = i_lon(t)

          gpos%time = i_time(t)
          gpos%lev = i_lev(t)
          gpos%lat = i_lat(t)
          gpos%lon = i_lon(t)

          CALL get_int_coords(t,gpos,ic,in_grid)   

          IF(in_grid)THEN
             CALL get_grid_vals(t,ic,wind,no_miss)
             IF(no_miss)THEN
                t_sp(t,c_step) = wind%sp

                ! Get Auxiliary values
                IF(naux >= 1)THEN
                   CALL get_aux_vals(t,ic)
                END IF

                ! Underground check
                IF(t_sp(t,c_step) <= t_lev(t,c_step))THEN
                   PRINT *, 'Trajectory encountered ground'
                   t_on(t) = .FALSE.
                   no_grnd = .FALSE.
                ELSE
                   t_u(t,c_step) = wind%u
                   t_v(t,c_step) = wind%v
                   t_w(t,c_step) = wind%w
 
                   ! Get Auxiliary values 
                   IF(naux >= 1)THEN
                      CALL get_aux_vals(t,ic)
                   END IF
                END IF
             END IF
          END IF

          IF(in_grid.EQV..FALSE. .OR. no_grnd.EQV..FALSE. .OR. no_miss.EQV..FALSE.)THEN
             t_u(t,c_step) = mv
             t_v(t,c_step) = mv
             t_w(t,c_step) = mv

             ! Set Auxiliary values to missing 
             IF(naux >= 1)THEN
                DO a = 1,naux
                   t_aux(t,c_step,a) = mv
                END DO
             END IF

             last_valid(t) = mvi  ! None valid 
          END IF

       END IF

    END DO

  END SUBROUTINE init_traj


  !###########################################################################
  ! petterson()                                                           ####
  ! Performs iterative time-step following Petterson (1940)               ####
  !###########################################################################

  SUBROUTINE petterson(t)
    IMPLICIT NONE

    ! Input/output variables
    INTEGER, INTENT(in) :: t                      ! Trajectory to compute

    ! Local variables
    TYPE(interp_coord) :: ic                      ! The output interpolation coordinates
    TYPE(geo_pos) :: gpos0,gpos1                  ! Geographic position (initial and next)
    TYPE(vec_wind) :: wind0, wind1, winda         ! 3D Wind (initial = 0, next = 1, average = a)
    LOGICAL :: in_grid,no_miss,no_grnd            ! Success flags
    INTEGER :: it                                 ! For iteration
    INTEGER :: ac                                 ! Auxiliary loop


    in_grid = .TRUE.
    no_miss = .TRUE.
    no_grnd = .TRUE.

    ! First Guess of X(0)
    wind0%u = t_u(t,c_step-1)
    wind0%v = t_v(t,c_step-1)
    wind0%w = t_w(t,c_step-1)

    gpos0%time = t_time(t,c_step-1)
    gpos0%lev = t_lev(t,c_step-1)
    gpos0%lat = t_lat(t,c_step-1)
    gpos0%lon = t_lon(t,c_step-1)

    CALL get_pos(wind0,gpos0,gpos1)           ! Find new gpos given initial vector wind
    CALL get_int_coords(t,gpos1,ic,in_grid)   ! Find fractional element position given gpos0
    IF(in_grid)THEN
       CALL get_grid_vals(t,ic,winda,no_miss) ! Given interpolation coordiantes get vector wind
       IF(no_miss)THEN
          IF(winda%sp <= gpos1%lev)THEN
             no_grnd = .FALSE.
          END IF
       END IF
    END IF

    ! Solve for X(1) iteratively
    DO it=1,3
       IF(in_grid .AND. no_grnd .AND. no_miss)THEN

          CALL get_pos(winda,gpos0,gpos1)               ! Find new pos given average of initial and guess
          CALL get_int_coords(t,gpos1,ic,in_grid)       ! Find interpolation coordinates at iteration i
          IF(in_grid)THEN
             CALL get_grid_vals(t,ic,wind1,no_miss)     ! Given interpolation coordiantes get vector wind at i
             IF(no_miss)THEN
               IF(winda%sp <= gpos1%lev)THEN
                  no_grnd = .FALSE.
               ELSE
                  ! Find average of initial and current guess             
                  winda%u = 0.5*(wind0%u+wind1%u)
                  winda%v = 0.5*(wind0%v+wind1%v)
                  winda%w = 0.5*(wind0%w+wind1%w)
               END IF
             END IF
          END IF

          ! IF(converged)THEN
          !   EXIT
          ! END IF

       END IF
    END DO

    ! Above steps have been successfully completed. Store values
    IF(in_grid .AND. no_grnd .AND. no_miss)THEN
       t_time(t,c_step) = gpos1%time
       t_lev(t,c_step) = gpos1%lev
       t_lat(t,c_step) = gpos1%lat
       t_lon(t,c_step) = gpos1%lon
       t_u(t,c_step) = wind1%u
       t_v(t,c_step) = wind1%v
       t_w(t,c_step) = wind1%w
       t_sp(t,c_step) = wind1%sp
       IF(naux >= 1)THEN
          CALL get_aux_vals(t,ic)
       END IF
    ELSE
       t_time(t,c_step) = mv
       t_lev(t,c_step) = mv
       t_lat(t,c_step) = mv
       t_lon(t,c_step) = mv
       t_u(t,c_step) = mv
       t_v(t,c_step) = mv
       t_w(t,c_step) = mv
       t_sp(t,c_step) = wind1%sp
       IF(naux >= 1)THEN
          DO ac = 1,naux
             t_aux(t,c_step,ac) = mv
          END DO
       END IF
       last_valid(t) = c_step - 1 
    END IF

  END SUBROUTINE petterson


  !###########################################################################
  ! get_int_coords()                                                      ####
  ! Use linear interpolation to convert the geographic postions to values ####
  ! used to read in a neigboorhood and the interpolation weighting to get ####
  ! the interpolated value from this neighboorhood fo values.
  !###########################################################################

  SUBROUTINE get_int_coords(t,gpos,ic,valid)
    IMPLICIT NONE

    ! Input/output variables
    INTEGER, INTENT(in) :: t                    ! Trajectory
    TYPE(geo_pos), INTENT(in) :: gpos           ! Geographic position
    TYPE(interp_coord), INTENT(out) :: ic       ! The output interpolation coordinates
    LOGICAL, INTENT(inout) :: valid             ! Success flag

    ! Local variables
    INTEGER :: p                                ! Counter 


    ! Interpolate time to fractional element postion
    IF( gpos%time < gtime(1) .OR. gpos%time > gtime(gntime) )THEN
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE 
       DO p=1,gntime-1
          IF( gpos%time >= gtime(p) .AND. gpos%time <= gtime(p+1)) THEN
             ic%i = 1.0*p+(gpos%time-gtime(p))/(gtime(p+1)-gtime(p))
             EXIT
          END IF
       END DO
    END IF

    ! Interpolate level to fractional element postion
    IF( gpos%lev > glev(1) .OR. gpos%lev < glev(gnlev) )THEN
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE 
       DO p=1,gnlev-1
          IF( gpos%lev <= glev(p) .AND. gpos%lev >= glev(p+1)) THEN
             ic%j = 1.0*p+(glev(p)-gpos%lev)/(glev(p)-glev(p+1))
             ic%p = gpos%lev  ! Actual pressure level
             EXIT
          END IF
       END DO
    END IF

    ! Interpolate latitude to fractional element postion
    IF( gpos%lat < glat(1) .OR. gpos%lat > glat(gnlat) )THEN
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE
       DO p=1,gnlat-1
          IF( gpos%lat >= glat(p) .AND. gpos%lat <= glat(p+1)) THEN
             ic%k = 1.0*p+(gpos%lat-glat(p))/(glat(p+1)-glat(p))
             EXIT
          END IF
       END DO
    END IF

    ! Interpolate longtude to fractional element postion
    IF( gpos%lon < glon(1) .OR. gpos%lon > glon(gnlon) )THEN
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE
       DO p=1,gnlon-1
          IF( gpos%lon >= glon(p) .AND. gpos%lon <= glon(p+1)) THEN
             ic%l = 1.0*p+(gpos%lon-glon(p))/(glon(p+1)-glon(p))
             EXIT
          END IF
       END DO
    END IF
    
    ! Determine neighborhood and weigthting for time
    IF( ic%i < 1.0 .OR. ic%i > gntime )THEN
       ! Out of bounds
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE IF( ic%i >= 1.0 .AND. ic%i <= (gntime) )THEN
       ! Central
       ic%nbi = (/0,1/) + INT(ic%i)
       ic%iwi = ic%i-INT(ic%i)
    END IF

   ! Determine neighborhood and weigthting for level
    IF( ic%j < 1.0 .OR. ic%j > gnlev )THEN
       ! Out of bounds
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE IF( ic%j >= 1.0 .AND. ic%j < 2.0)THEN
       ! Right end
       ic%nbj = (/1,2,3,4/)
       ic%iwj = ic%j-INT(ic%j)-1.0
       ic%ploc = glev( ic%nbj )  ! Actual pressure levels
    ELSE IF( ic%j >= (gnlev-1) .AND. ic%j <= gnlev)THEN
       ! Left end
       ic%nbj = (/gnlev-3,gnlev-2,gnlev-1,gnlev/)
       ic%iwj = ic%j-INT(ic%j)+1.0
       ic%ploc = glev( ic%nbj )  ! Actual pressure levels
    ELSE IF( ic%j >= 2.0 .AND. ic%j <= (gnlev-1) )THEN
       ! Central
       ic%nbj = (/-1,0,1,2/) + INT(ic%j)
       ic%iwj = ic%j-INT(ic%j)
       ic%ploc = glev( ic%nbj )  ! Actual pressure levels
    END IF

   ! Determine neighborhood and weigthting for latitude
    IF( ic%k < 1.0 .OR. ic%k > gnlat )THEN
       ! Out of bounds
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE IF( ic%k >= 1.0 .AND. ic%k < 2.0)THEN
       ! Right end
       ic%nbk = (/1,2,3,4/)
       ic%iwk = ic%k-INT(ic%k)-1.0
    ELSE IF( ic%k >= (gnlat-1) .AND. ic%k <= gnlat)THEN
       ! Left end
       ic%nbk = (/gnlat-3,gnlat-2,gnlat-1,gnlat/)
       ic%iwk = ic%k-INT(ic%k)+1.0
    ELSE IF( ic%k >= 2.0 .AND. ic%k <= (gnlat-1) )THEN
       ! Central
       ic%nbk = (/-1,0,1,2/) + INT(ic%k)
       ic%iwk = ic%k-INT(ic%k)
    END IF

   ! Determine neighborhood and weigthting for longitude
    IF( ic%l < 1.0 .OR. ic%l > gnlon )THEN
       ! Out of bounds
       valid = .FALSE.
       t_on(t) = .FALSE.
    ELSE IF( ic%l >= 1.0 .AND. ic%l < 2.0)THEN
       ! Right end
       ic%nbl = (/1,2,3,4/)
       ic%iwl = ic%l-INT(ic%l)-1.0
    ELSE IF( ic%l >= (gnlon-1) .AND. ic%l <= gnlon)THEN
       ! Left end
       ic%nbl = (/gnlon-3,gnlon-2,gnlon-1,gnlon/)
       ic%iwl = ic%l-INT(ic%l)+1.0
    ELSE IF( ic%l >= 2.0 .AND. ic%l <= (gnlon-1) )THEN
       ! Central
       ic%nbl = (/-1,0,1,2/) + INT(ic%l)
       ic%iwl = ic%l-INT(ic%l)
    END IF


  END SUBROUTINE get_int_coords


!###########################################################################
! get_pos()                                                             ####
! Calculates gpos1 given vector wind and gpos0 using Haversine's        ####
! formula.                                                              ####
!###########################################################################

  SUBROUTINE get_pos(wind,gpos0,gpos1)
    IMPLICIT NONE

    ! Input/output variables
    TYPE(geo_pos), INTENT(in) :: gpos0        ! Initial geographic position
    TYPE(geo_pos), INTENT(out) :: gpos1       ! Next geographic position
    TYPE(vec_wind), INTENT(in) :: wind        ! 3D vector wind (initial or average depending on call)

    ! Local Variables
    REAL, PARAMETER :: pi = 3.14159265358979  ! pi
    REAL, PARAMETER :: radius = 6371220.0     ! Earth radius (m)
    REAL :: latr0, lonr0, latr1, lonr1        ! Initial and guess lat/lon (radians)
    REAL :: dir, rdist, adj                   ! Direction, radial distance, backward/forward switch


    ! Reverse direction if time_step < 0
    IF(time_step < 0)THEN
       adj = -1.0
    ELSE
       adj = 1.0
    END IF

    ! Current lat and lon in radians
    latr0 = gpos0%lat*pi/180.0
    lonr0 = gpos0%lon*pi/180.0

    ! Solve for radial distance of the motion over 1 time step and the bearing
    ! of the wind based on u and v. Bearing is reversed for back trajectories.
    rdist = SQRT( wind%u**2.0+wind%v**2.0 )*(adj*time_step/radius)
    dir = pi+ATAN2( adj*-1.0*wind%u, adj*-1.0*wind%v )

    ! Calculate new lat and lon using Haversine's forumula
    latr1 = ASIN( SIN(latr0)*COS(rdist)+COS(latr0)*SIN(rdist)*COS(dir) )
    lonr1 = lonr0+ATAN2( SIN(dir)*SIN(rdist)*COS(latr0), COS(rdist)-SIN(latr0)*SIN(latr1) )

    ! Convert to degrees
    gpos1%lat = latr1*180.0/pi
    gpos1%lon= lonr1*180.0/pi

    ! New time
    gpos1%time = gpos0%time+time_step

    ! New vertical position
    gpos1%lev =  gpos0%lev + wind%w*time_step

  END SUBROUTINE get_pos


!###########################################################################
! get_grid_vals()                                                       ####
! Given the interpolation coordinates a neighborhood of values is read  ####
! in and used to find an interpolated value given the input weightings. ####
!###########################################################################

  SUBROUTINE get_grid_vals(t,ic,wind,no_miss)                     
    IMPLICIT NONE

    ! Input/Output Variables
    INTEGER, INTENT(in) :: t                           ! Trajectory number
    TYPE(interp_coord), INTENT(in) :: ic               ! Interpolation coordinates
    TYPE(vec_wind), INTENT(inout) :: wind              ! 3D wind interpolated from wind_nbor and fep
    LOGICAL, INTENT(inout) :: no_miss                  ! Success flag

    ! Local Variables
    TYPE(vec_wind), DIMENSION(4,4,4,2) :: wind_nbor    ! Temporarily stores wind neighboorhood to get interpolated wind
    REAL, DIMENSION(4,4,2) :: sp_nbor


    ! Get 4D neighborhood values
    wind_nbor%u = get_nbor_xyzt(ic%nbi,ic%nbj,ic%nbk,ic%nbl,gncid,u_varid)
    wind_nbor%v = get_nbor_xyzt(ic%nbi,ic%nbj,ic%nbk,ic%nbl,gncid,v_varid)
    wind_nbor%w = get_nbor_xyzt(ic%nbi,ic%nbj,ic%nbk,ic%nbl,gncid,w_varid)
    sp_nbor = get_nbor_xyt(ic%nbi,ic%nbk,ic%nbl,gncid,sp_varid)

    ! Check for missing values. If there are missing values interpolation
    ! is not performed and all values are set to missing. Also the trajectory
    ! is turned off.

    IF( ANY(wind_nbor%u == mv) .OR. ANY(wind_nbor%v == mv) .OR. ANY(wind_nbor%w == mv) .OR. ANY(wind_nbor%sp == mv) )THEN
       PRINT *, 'Missing values encountered'
       no_miss = .FALSE.
       t_on(t) = .FALSE.
    ELSE
       ! Get interpolated value given weightings
       wind%u = interp_4d(ic, wind_nbor%u)
       wind%v = interp_4d(ic, wind_nbor%v)
       wind%w = interp_4d(ic, wind_nbor%w)
       wind%sp = interp_3d(ic, sp_nbor)
    END IF

  END SUBROUTINE get_grid_vals


!###########################################################################
! get_aux_vals()                                                        ####
! Given the interpolation coordinates a neighborhood of values is read  ####
! in and used to find an interpolated value given the input weightings. ####
!###########################################################################

  SUBROUTINE get_aux_vals(t,ic)                     
    IMPLICIT NONE

    ! Input/Output Variables
    INTEGER, INTENT(in) :: t                 ! Trajectory number
    TYPE(interp_coord), INTENT(in) :: ic     ! Interpolation coordinates

    ! Local Variables
    REAL, DIMENSION(4,4,4,2) :: a4d_xyzt     ! Temporarily stores 4D auxiliary (e.g. pv)
    REAL, DIMENSION(4,4,2) :: a3d_xyt        ! Temporaily stores 3D auxiliary (e.g. mslp)
    INTEGER :: ac                            ! Loop


    ! Loop through the auxililiary values
    DO ac = 1,naux
       ! If 4D auxiliary value
       IF(aux_dim(ac) == 4)THEN
          ! Get neighborhood
          a4d_xyzt = get_nbor_xyzt(ic%nbi,ic%nbj,ic%nbk,ic%nbl,ancid,aux_varid(ac))
          IF( ANY(a4d_xyzt == mv) )THEN
             t_aux(t,c_step,ac) = mv
          ELSE
             ! Get interpolated value given weightings
             t_aux(t,c_step,ac) = interp_4d(ic,a4d_xyzt)
          END IF
       END IF

       ! If 3D auxiliary value
       IF(aux_dim(ac) == 3)THEN
          ! Get neighborhood
          a3d_xyt = get_nbor_xyt(ic%nbi,ic%nbk,ic%nbl,ancid,aux_varid(ac))
          IF( ANY(a3d_xyt == mv) )THEN
             t_aux(t,c_step,ac) = mv
          ELSE
             ! Get interpolated value given weightings
             t_aux(t,c_step,ac) = interp_3d(ic,a3d_xyt)
          END IF
       END IF
    END DO

  END SUBROUTINE get_aux_vals


  !###########################################################################
  ! get_nbor_xyzt()                                                       ####
  ! Reads in a 4d grid to be used to find var given g* elements.          ####
  !###########################################################################

  FUNCTION get_nbor_xyzt(gis,gjs,gks,gls,ncid,var_varid)
    IMPLICIT NONE

    ! Input Variables
    INTEGER, INTENT(IN), DIMENSION(4) :: gjs, gks, gls
    INTEGER, INTENT(IN), DIMENSION(2) :: gis
    INTEGER, INTENT(IN) ::  ncid, var_varid

    ! Local Variables
    INTEGER, PARAMETER :: NDIMS = 4
    INTEGER :: start(NDIMS), COUNT(NDIMS)
    REAL, DIMENSION(4,4,4,2) :: get_nbor_xyzt


    start = (/gls(1),gks(1),gjs(1),gis(1)/)  ! Reverse order for NetCDF
    count = (/4,4,4,2/)                      ! Number in each dimension.

    ! Read in neighborhood from NetCDF
    CALL check( nf90_get_var(ncid, var_varid, get_nbor_xyzt, start, count ) )

  END FUNCTION get_nbor_xyzt


  !###########################################################################
  ! get_nbor_xyt()                                                        ####
  ! Reads in a 3d grid to be used to find var given g* elements.          ####
  !###########################################################################

  FUNCTION get_nbor_xyt(gis,gks,gls,ncid,var_varid)
    IMPLICIT NONE

    ! Input Variables
    INTEGER, INTENT(IN), DIMENSION(4) :: gks,gls
    INTEGER, INTENT(IN), DIMENSION(2) :: gis
    INTEGER, INTENT(IN) ::  ncid, var_varid

    ! Local Variables
    INTEGER, PARAMETER :: NDIMS = 3
    INTEGER :: start(NDIMS), COUNT(NDIMS)
    REAL, DIMENSION(4,4,2) :: get_nbor_xyt


    start = (/gls(1),gks(1),gis(1)/)   ! Reverse order for NetCDF
    count = (/4,4,2/)                  ! Number in each dimension.

    ! Read in neighborhood from NetCDF
    CALL check( nf90_get_var(ncid, var_varid, get_nbor_xyt, start, count ) )

  END FUNCTION get_nbor_xyt


  !###########################################################################
  ! interp_4d()                                                           ####
  ! Given the g*s arrays and the g* arrays stored in memory the value     ####
  ! at the fractional elements g* are determined.                         ####
  !###########################################################################

  FUNCTION interp_4d(ic, var_xyzt)
    IMPLICIT NONE

    ! Input Variables
    TYPE(interp_coord), INTENT(in) :: ic                 ! Interpolation coordinates
    REAL, INTENT(IN),  DIMENSION(4,4,4,2) :: var_xyzt    ! Neighborhood of values

    ! Local Variables
    REAL :: interp_4d, t1, t2
    INTEGER :: i
    REAL, DIMENSION(4) :: prof1, prof2, psub


    ! Generate two profiles at t1 and t2 using bicubic interpolation
    DO i=1,4
      prof1(i) = bicubic_interpolate_xy(var_xyzt(1:4,1:4,i,1), ic%iwl, ic%iwk)
      prof2(i) = bicubic_interpolate_xy(var_xyzt(1:4,1:4,i,2), ic%iwl, ic%iwk)
    END DO

    ! Determine the t1 and t2 for respective profiles using Newton polynomials
    t1 = neville_interpolate_p(ic%ploc, prof1, ic%p)
    t2 = neville_interpolate_p(ic%ploc, prof2, ic%p)

    ! Linear interpolation for time
    interp_4d = t1*(1.0-ic%iwi) + t2*ic%iwi

  END FUNCTION interp_4d


  !###########################################################################
  ! interp_3d()                                                           ####
  ! Given the g*s arrays and the g* arrays stored in memory the value     ####
  ! at the fractional elements g* are determined.                         ####
  !###########################################################################

  FUNCTION interp_3d(ic, var_xyt)
    IMPLICIT NONE

    ! Input Variables
    TYPE(interp_coord), INTENT(in) :: ic                 ! Interpolation coordinates
    REAL, INTENT(IN),  DIMENSION(4,4,2) :: var_xyt       ! Neighborhood of values

    ! Local Variables
    REAL :: interp_3d, t1, t2


    ! Generate values at t1 and t2 using bicubic interpolation
    t1 = bicubic_interpolate_xy(var_xyt(1:4,1:4,1), ic%iwl, ic%iwk)
    t2 = bicubic_interpolate_xy(var_xyt(1:4,1:4,2), ic%iwl, ic%iwk)

    ! Linear interpolation for time
    interp_3d = t1*(1.0-ic%iwi) + t2*ic%iwi

  END FUNCTION interp_3d


  !###########################################################################
  ! cubic_interpolate()                                                   ####
  ! Perform cubic interpolation to find value at fractional element       ####
  !###########################################################################

  ! In the cubic interpolation functions p is an array of values (e.g. wind)
  ! w, x, y, and z are fractional elements (/-1:2/) with a center point defined
  ! as 0.5.

  FUNCTION cubic_interpolate(p, w)
    IMPLICIT NONE

    ! Input Variables
    REAL, INTENT(IN) :: w
    REAL, INTENT(IN), DIMENSION(4) :: p

    ! Local Variables
    REAL :: cubic_interpolate


    cubic_interpolate = p(2) + 0.5*w*(p(3) - p(1) + w*(2.0*p(1) - 5.0*p(2) + 4.0*p(3) - p(4) + w*(3.0*(p(2) - p(3)) + p(4) - p(1))))

  END FUNCTION cubic_interpolate


  !###########################################################################
  ! bicubic_interpolate_xy()                                              ####
  ! Perform cubic interpolation to find value at fractional element       ####
  !###########################################################################

  FUNCTION bicubic_interpolate_xy(p, w, x)
    IMPLICIT NONE

    ! Input Variables
    REAL, INTENT(IN) :: w, x
    REAL, INTENT(IN), DIMENSION(4,4) :: p

    ! Local Variables
    REAL, DIMENSION(4) :: p_temp
    REAL :: bicubic_interpolate_xy


    p_temp(1) = cubic_interpolate(p(1,:), x)
    p_temp(2) = cubic_interpolate(p(2,:), x)
    p_temp(3) = cubic_interpolate(p(3,:), x)
    p_temp(4) = cubic_interpolate(p(4,:), x)
    bicubic_interpolate_xy = cubic_interpolate(p_temp, w)

  END FUNCTION bicubic_interpolate_xy


  !###########################################################################
  ! neville_interpolate_p()                                               ####
  ! Perform 3rd order Newton polynomial interpolation using               ####
  ! Neville's algortihin                                                  ####
  !###########################################################################

  FUNCTION neville_interpolate_p(psub, prof, inp)
    IMPLICIT NONE

    ! Input Variables
    REAL, INTENT(IN) :: inp
    REAL, INTENT(IN), DIMENSION(4) :: psub, prof

    ! Local Variables
    REAL, DIMENSION(4) :: f
    REAL :: neville_interpolate_p
    INTEGER :: n, i, j


    n = 4   ! Number of values in array

    ! Solve for the unique polynomial of order n-1

    f = prof
   
    DO j=2,n
       DO i=n,j+1,-1
          f(i) = ( (inp-psub(i-j))*f(i) - (inp-psub(i))*f(i-1) ) / ( psub(i)-psub(i-j) )
       END DO
    END DO

    neville_interpolate_p = f(n)

  END FUNCTION neville_interpolate_p



END MODULE mod_calctraj
