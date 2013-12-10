  subroutine check(status)
  use netcdf
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check


    PROGRAM read_netcdf
    USE netcdf
    IMPLICIT NONE
!    real lev(:)
    character (len = *), parameter :: FILE_NAME = "/free4/abrammer/WRFOUT/10d50km/wrfout_d01_2010-07-05_00:00:00"
    character (len=*), parameter :: FMT1 = "(7F15.5)", FMT0="(7A15)"
    character (len = 25)  lat_name, lon_name, lev_name, tim_name, var_name
    integer, parameter :: text_out = 20
    real, dimension(:), allocatable :: lat, lon, lev, time
    integer ncid, varId, dimId, ndim, dimlen, uId, inds(4), ninds(4),i, numAtts, it
    integer, dimension(nf90_max_var_dims) :: dimIds
    real ti(3), li(3), loni(3),lati(3), var(4,4,4), pro1(4), val
	real start_time, start_lon, start_lat, start_lev, end_time
	real t2, t1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! External Functions
    real bicubic_interp, neville_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Lets get Rolling.

    type wind
    	real :: val1, val2, val
    	real, dimension(:), allocatable :: lat
    	real, dimension(:), allocatable :: lon
    	real, dimension(:), allocatable :: lev
    	character (len = 25):: units, name
    	integer :: id, y_stag
    end type wind
    type (wind) u
    type (wind) v
    type (wind) w

    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

    var_name = "XTIME"
    call check( nf90_inq_varid(ncid, var_name, varId ) )
    call check( nf90_inquire_variable(ncid, varId, ndims = ndim, dimids = dimIds, natts = numAtts) )
    call check( nf90_inquire_variable(ncid, varId, dimids = dimIds(:ndim) ) )
    call check( nf90_inquire_dimension(ncid, dimids(1), len=dimlen) )
    allocate ( time(dimlen) )
    call check (nf90_get_var(ncid, varId, time) )

!    ***********************************
!    define dimensions

	u%name = "U"
	call define_coords(u)
	print*, size(u%lon), size(u%lat),u%y_stag
	print*, maxval(u%lon), maxval(u%lat)
	print*, minval(u%lon), minval(u%lat)

	v%name = "V"
	call define_coords(v)
	print*, size(v%lon), size(v%lat),v%y_stag
	print*, maxval(v%lon), maxval(v%lat)
	print*, minval(v%lon), minval(v%lat)

	w%name = "W"
	call define_coords(w)
	print*, size(w%lon), size(w%lat), w%y_stag
	print*, maxval(w%lon), maxval(w%lat)
	print*, minval(w%lon), minval(w%lat)


!    ***********************************
!   starting point

    start_time = 170.
    end_time = 190.
    start_lon = -65.3
    start_lat = 15.2
    start_lev = 15431.0


    call integrate_trajectory(1, start_time, end_time, start_lev, start_lat, start_lon)

	contains

	REAL function torad(val)
	REAL, PARAMETER :: pi = 3.14159265358979  ! pi
	real val
    torad = val* pi/180.
    end function torad

    REAL function todeg(val)
	REAL, PARAMETER :: pi = 3.14159265358979  ! pi
	real val
    todeg = val* 180./pi
    end function todeg


	real function move_parcel(u,v,w,lev_in, lat_in, lon_in)
	real u,v,w
	REAL, PARAMETER :: pi = 3.14159265358979  ! pi
    REAL, PARAMETER :: radius = 6371220.0     ! Earth radius (m)
    dimension move_parcel(3)
    REAL :: latr0, lonr0, latr1, lonr1, lat_in, lon_in       ! Initial and guess lat/lon (radians)
    REAL :: dir, rdist, adj , newlat, newlon, time_step, newlev, lev_in          ! Direction, radial distance, backward/forward switch

    time_step = 30
    adj = 1.

	rdist = SQRT( u**2.0+v**2.0 )*(adj*time_step/radius)      ! is adj needed?
    dir = pi+ATAN2( adj*-1.0*u, adj*-1.0*v )

    latr0 = torad(lat_in)
    lonr0 = torad(lon_in)

    ! Calculate new lat and lon using Haversine's forumula
    latr1 = ASIN( SIN(latr0)*COS(rdist)+COS(latr0)*SIN(rdist)*COS(dir) )
    lonr1 = lonr0+ATAN2( SIN(dir)*SIN(rdist)*COS(latr0), COS(rdist)-SIN(latr0)*SIN(latr1) )

    newlev = lev_in + (adj*time_step)*w

    ! Convert to degrees
    newlat = todeg(latr1)
    newlon = todeg(lonr1)
    move_parcel = (/newlev, newlat, newlon/)
    end function move_parcel


    subroutine integrate_trajectory(tr, start_time, end_time, start_lev, start_lat, start_lon)
    real start_time, end_time, start_lev, start_lat, start_lon
    integer tr

    open(unit=text_out, file="traj1.txt", status="replace", action="write")
    write(text_out,FMT0) "TIME","LEV","LAT","LON","U","V","W"

	call get5dval(ncid, u,v,w, start_time, start_lev, start_lat, start_lon)
    write(text_out,FMT1) start_time, start_lev, start_lat, start_lon, u%val1, v%val1, w%val1

    call cpu_time ( t1 )
    do while (start_time .lt. end_time)
       call petterson(u,v,w, start_time, start_lev, start_lat, start_lon)
       write(text_out,FMT1) start_time, start_lev, start_lat, start_lon, u%val1, v%val1, w%val1
    end do
    call cpu_time ( t2 )
    write ( *, * ) 'Elapsed CPU time = ', t2 - t1
    end subroutine integrate_trajectory


    subroutine petterson(u,v,w,in_time,lev_in, lat_in, lon_in)
    type (wind) :: u,v,w, u0, v0, w0
    real locs(3), newlocs(3), lev_in, lat_in, lon_in, in_time
    integer it
! move parcel based on initial vector
! Get winds at new location.
! move parcel based on avg of initial and new vector
! iterate with updated new vector until convergent.
    u0 = u
    v0 = v
    w0 = v

    newlocs = move_parcel(u%val1,v%val1,w%val1, lev_in, lat_in, lon_in)
    in_time = in_time + 0.5
    call get5dval(ncid, u,v,w, in_time, newlocs(1), newlocs(2), newlocs(3) )  ! get initial guess

    do it=1, 3
    locs  = newlocs
    newlocs = move_parcel(0.5*(u0%val1+u%val1),0.5*(v0%val1+v%val1),0.5*(w0%val1+w%val1), lev_in, lat_in, lon_in)
    call get5dval(ncid, u,v,w, in_time, newlocs(1), newlocs(2), newlocs(3) )
    if(maxval(locs - newlocs).eq.0)then ! check for convergence.
        exit
    end if
    end do
    lev_in = newlocs(1)
    lat_in = newlocs(2)
    lon_in = newlocs(3)
    end subroutine petterson



	subroutine define_coords(vari)
	type (wind) vari
	integer varId, tId
	character (len =25) dimname, newname
	real, dimension(:), allocatable :: lonvar, latvar
    call check( nf90_inq_varid(ncid, vari%name, vari%id ) )
    call check( nf90_inquire_variable(ncid, vari%id, ndims = ndim, dimids = dimIds, natts = numAtts) )
    call check( nf90_inquire_variable(ncid, vari%id, dimids = dimIds(:ndim) ) )
    do i=1,ndim
		call check( nf90_inquire_dimension(ncid, dimids(i),name = dimname, len=dimlen) )
		if(dimname == "west_east")then
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(ncid, "XLONG", tId) )
			call check (nf90_get_var(ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		end if
		if(dimname == "west_east_stag")then
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(ncid, "XLONG_U", tId ) )
			call check (nf90_get_var(ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		end if
		if(dimname == "south_north")then
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(ncid, "XLAT", tId ) )
			call check (nf90_get_var(ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
			vari%lat = latvar
		end if
		if(dimname == "south_north_stag")then
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(ncid, "XLAT_V", tId ) )
			call check (nf90_get_var(ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
			vari%lat = latvar
		end if
		if(dimname == "bottom_top")then
		 allocate (vari%lev(dimlen))
		vari%y_stag = 0
		end if
		if(dimname == "bottom_top_stag")then
		vari%y_stag = 1
		allocate (vari%lev(dimlen))
		end if
    end do
	return
	end subroutine define_coords


    real function coord_2_int(x, xx)
    real x(:), xx, i, int
    dimension coord_2_int(3)

    do i=1, size(x)-1
       if(x(2) .gt. x(1)) then
            if( xx >= x(i) .and. xx < x(i+1) )then
            int = 1.0* (xx - x(i)) / (x(i+1) - x(i) )
            exit
            end if
        else
            if( xx <= x(i) .and. xx > x(i+1) )then
            int = 1.0* (xx - x(i)) / (x(i+1) - x(i) )
            exit
            end if
        end if
    end do
    coord_2_int(1) = int
    coord_2_int(2) = i
    coord_2_int(3) = ( x(i+1) - x(i) )
    end function coord_2_int


	subroutine get_levels(ncid, u,v,w, in_time, in_lat, in_lon)
	real, parameter :: g=9.81
	type (wind) :: u,v,w
	real loni(3), lati(3),  in_lat, in_lon, in_time
	real, dimension(:,:,:), allocatable :: define_levels, pert_levels
	real, dimension(:), allocatable :: profile, ustag_lev
	integer inds(4), ncid, l, i
	integer varId, tId, coords(2), ndim
	character (len =25) dimname, newname
	real bicubic_interp, neville_interp

	 call check( nf90_inq_varid(ncid, "PHB", tId ) )
	 call check( nf90_inquire_variable(ncid, tId, ndims = ndim, dimids = dimIds) )
	 do i=1,ndim
		 call check( nf90_inquire_dimension(ncid, dimids(i),name = dimname, len=dimlen) )
		 if(dimname == "bottom_top_stag")then
			allocate (define_levels(4,4,dimlen))
			allocate (pert_levels(4,4,dimlen))
			allocate (profile(dimlen))
			define_levels = (define_levels + pert_levels) / g
			loni =  coord_2_int(w%lon, in_lon)
   			lati =  coord_2_int(w%lat, in_lat)
   			inds = int( (/loni(2)-1, lati(2)-1, 1, 2 /) )
			call check (nf90_get_var(ncid, tId, define_levels, inds, (/ 4,4,dimlen,1/) ) )
			call check( nf90_inq_varid(ncid, "PH", tId ) )
			call check (nf90_get_var(ncid, tId, pert_levels, inds, (/ 4,4,dimlen,1/) ) )
		   do l=1,dimlen
		   profile(l) = bicubic_interp(define_levels(:,:,l), lati(1), loni(1) )
		   end do
		   w%lev = profile
			exit
		end if
	 end do
	 allocate (ustag_lev(dimlen-1))
	 ustag_lev =  (profile(2:dimlen) + profile(1:dimlen-1) )/ 2.
	 u%lev = ustag_lev
	 v%lev = ustag_lev
    end subroutine get_levels

	subroutine get5dval(ncid, u,v,w, in_time, lev, lat, lon)
	type (wind) u,v,w
	real in_time, lev, lat, lon
	integer ncid
	real linear_interp
	ti =  coord_2_int(time, in_time)
	call get_levels(ncid, u,v,w, time(ti(2)), lat, lon)
	u%val1 = get4dval(ncid, u, time(ti(2)), lev, lat, lon)
	v%val1 = get4dval(ncid, v, time(ti(2)),lev, lat, lon)
	w%val1 = get4dval(ncid, w, time(ti(2)), lev, lat, lon)
	if(ti(1).ne.0)then
	  call get_levels(ncid, u,v,w, time(ti(2)+1), lat, lon)
	  u%val2 = get4dval(ncid, u, time(ti(2)+1),lev, lat, lon)
	  v%val2 = get4dval(ncid, v, time(ti(2)+1),lev, lat, lon)
	  w%val2 = get4dval(ncid, w, time(ti(2)+1),lev, lat, lon)
	  u%val1 = linear_interp( u%val1, u%val2, ti(1) )
	  v%val1 = linear_interp( v%val1, v%val2, ti(1) )
	  w%val1 = linear_interp( w%val1, w%val2, ti(1) )
	end if
	end subroutine get5dval


	real function get4dval(ncid, reqvar, in_time, in_lev, in_lat, in_lon)
	real ti(3), levi(3), loni(3), lati(3),  in_time, in_lat, in_lon, in_lev
	real pro1(4), var(4,4,4)
	integer inds(4), ninds(4), ncid
	type (wind) reqvar
	real bicubic_interp, neville_interp
	ninds = (/4,4,4,1/)

    ti =  coord_2_int(time, in_time)
    levi =  coord_2_int(reqvar%lev, in_lev)
    loni =  coord_2_int(reqvar%lon, in_lon)
    lati =  coord_2_int(reqvar%lat, in_lat)
    inds = int( (/loni(2)-1, lati(2)-1, levi(2)-1, ti(2) /) )

    call check (nf90_get_var(ncid, reqvar%id, var, inds, ninds ) )

	!  read in bottom left - bottom right - top left -  topright
	do i=1,4
	pro1(i) = bicubic_interp(var(:,:,i), lati(1), loni(1) )
	end do
	val= neville_interp(reqvar%lev(levi(2)-1:levi(2)+2), pro1, start_lev, size(pro1) )


	if(ti(1).ne.0)then
		print*, "not zero"
	end if

    get4dval = val
    end function get4dval


	end PROGRAM
