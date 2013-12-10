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
    character (len = *), parameter :: VAR_FILE = "/glade/scratch/abrammer/runs/nadine_2012/run3/wind_d03_2012-09-06_00:00:00.nc"
    character (len = *), parameter :: META_FILE = "/glade/scratch/abrammer/runs/nadine_2012/run3/wrfout_d03_2012-09-06_00:00:00.nc"
    character (len = *), parameter :: namefile = "namelist.traj"
    character (len=*), parameter :: FMT1 = "(7F15.5)", FMT0="(7A15)"
    character (len = 25)  lat_name, lon_name, lev_name, tim_name, var_name
    real, dimension(:), allocatable :: lat, lon, lev, time

    integer, parameter :: text_out = 20
    integer  meta_ncid, ncid, varId, dimId, ndim, dimlen, uId, inds(4), ninds(4),i, numAtts, it,x,t, step
    integer, dimension(nf90_max_var_dims) :: dimIds

    real ti(3), li(3), loni(3),lati(3), var(4,4,4), pro1(4), val
	real start_time, start_lon, start_lat, start_lev, end_time
	real t2, t1
	namelist /time_opt/start_time,end_time,step
	namelist /traj_opt/start_lev,start_lat,start_lon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! External Functions
    real bicubic_interp, neville_interp
    integer istat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Lets get Rolling.
    type wind
    	real :: val1, val2, val
    	real, dimension(:), allocatable :: lat
    	real, dimension(:), allocatable :: lon
    	real, dimension(:), allocatable :: lev
    	real, dimension(:,:,:,:), allocatable :: grid
    	character (len = 25):: units, name
    	integer :: id, y_stag
    end type wind
    type (wind) u
    type (wind) v
    type (wind) w

    istat = 1
    call check( nf90_open(VAR_FILE, NF90_NOWRITE, ncid) )
    call check( nf90_open(META_FILE, NF90_NOWRITE, meta_ncid) )

    !var_name = "XTIME"
    !call check( nf90_inq_varid(ncid, var_name, varId ) )
    !call check( nf90_inquire_variable(ncid, varId, ndims = ndim, dimids = dimIds, natts = numAtts) )
    !call check( nf90_inquire_variable(ncid, varId, dimids = dimIds(:ndim) ) )
    !call check( nf90_inquire_dimension(ncid, dimids(1), len=dimlen) )
    !allocate ( time(dimlen) )
    !call check (nf90_get_var(ncid, varId, time) )

     allocate(time(12))
     do t=1,12
        time(t) = (t*5)-5
        end do

!    ***********************************
!    define dimensions

	u%name = "U"
	call define_coords(u)
	print*, "U Variable"
	print*, "Longitude  (", size(u%lon),"): ",minval(u%lon),":",maxval(u%lon)
	print*, "Latitude  (", size(u%lat),"): ",minval(u%lat),":",maxval(u%lat)

	v%name = "V"
	call define_coords(v)
	print*, "V Variable"
	print*, "Longitude  (", size(v%lon),"): ",minval(v%lon),":",maxval(v%lon)
	print*, "Latitude  (", size(v%lat),"): ",minval(v%lat),":",maxval(v%lat)

	w%name = "W"
	call define_coords(w)
	print*, "W Variable"
	print*, "Longitude  (", size(w%lon),"): ",minval(w%lon),":",maxval(w%lon)
	print*, "Latitude  (", size(w%lat),"): ",minval(w%lat),":",maxval(w%lat)




!    ***********************************
!   starting point
        open(1,file=namefile)
        read(1,time_opt)
        read(1,traj_opt)


    call cpu_time ( t1 )
    call read_time()
    call cpu_time ( t2 )
    write ( *, * ) 'Elapsed CPU time = ', t2 - t1
    print*, "*********************"

    call integrate_trajectory(1, start_time, end_time, start_lev, start_lat, start_lon)




	contains

	subroutine read_time()
    use netcdf
    integer tId
    real vara(2226, 1125, 51)
	call check( nf90_inq_varid(ncid, "PH", tId ) )
	call check (nf90_get_var(ncid, tId, vara, (/1,1,1,1/), (/ 2226,1125,51,1/) ) )
    end subroutine


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
			call check( nf90_inq_varid(meta_ncid, "XLONG", tId) )
			call check (nf90_get_var(meta_ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		end if
		if(dimname == "west_east_stag")then
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLONG_U", tId ) )
			call check (nf90_get_var(meta_ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		end if
		if(dimname == "south_north")then
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLAT", tId ) )
			call check (nf90_get_var(meta_ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
			vari%lat = latvar
		end if
		if(dimname == "south_north_stag")then
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLAT_V", tId ) )
			call check (nf90_get_var(meta_ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
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

    subroutine integrate_trajectory(tr, start_time, end_time, start_lev, start_lat, start_lon)
    real start_time, end_time, start_lev, start_lat, start_lon
    integer tr
    character (len=3) :: stri
    write(stri,'(I3)') tr
    open(unit=text_out, file="traj"//trim(ADJUSTL(stri))//".txt", status="replace", action="write")
    write(text_out,FMT0) "TIME","LEV","LAT","LON","U","V","W"

	call get5dval(ncid, u,v,w, start_time, start_lev, start_lat, start_lon)
    write(text_out,FMT1) start_time, start_lev, start_lat, start_lon, u%val1, v%val1, w%val1

    do while (start_time .lt. end_time)
    call cpu_time ( t1 )
       call petterson(u,v,w, start_time, start_lev, start_lat, start_lon)
    call cpu_time ( t2 )
    write ( *, * ) 'Elapsed CPU time = ', t2 - t1
       if(.not.istat)then
		return
	end if
    if( mod(start_time, 5.).eq.0 )then
	write(text_out,FMT1) start_time, start_lev, start_lat, start_lon, u%val1, v%val1, w%val1
    end if
    end do
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
    w0 = w

    newlocs = move_parcel(u%val1,v%val1,w%val1, lev_in, lat_in, lon_in)
    in_time = in_time + step/60.
    call get5dval(ncid,u,v,w, in_time, newlocs(1), newlocs(2), newlocs(3) )  ! get initial guess

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


    real function move_parcel(u,v,w,lev_in, lat_in, lon_in)
    real u,v,w
    REAL, PARAMETER :: pi = 3.14159265358979  ! pi
    REAL, PARAMETER :: radius = 6371220.0     ! Earth radius (m)
    dimension move_parcel(3)
    REAL :: latr0, lonr0, latr1, lonr1, lat_in, lon_in       ! Initial and guess lat/lon (radians)
    REAL :: dir, rdist, adj , newlat, newlon, time_step, newlev, lev_in          ! Direction, radial distance, backward/forward switch

    time_step = step
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


	subroutine get5dval(ncid, u,v,w, in_time, lev, lat, lon)
	type (wind) u,v,w
	real in_time, lev, lat, lon, ti(3)
	integer ncid
	real linear_interp
	call coord_2_int(time, in_time, ti)

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
	 call check( nf90_inq_varid(meta_ncid, "PHB", tId ) )
	 call check( nf90_inquire_variable(meta_ncid, tId, ndims = ndim, dimids = dimIds) )
	 do i=1,ndim
		 call check( nf90_inquire_dimension(meta_ncid, dimids(i),name = dimname, len=dimlen) )
		 if(dimname == "bottom_top_stag")then
			allocate (define_levels(4,4,dimlen))
			allocate (pert_levels(4,4,dimlen))
			allocate (profile(dimlen))
			call coord_2_int(w%lon, in_lon,loni)
   			call coord_2_int(w%lat, in_lat, lati)
   			call coord_2_int(time, in_time, ti)

   			if(.not.istat)then
				print*, istat, "Error in get_levels"
				print*, loni, lati
				return
			end if

			inds = int( (/loni(2)-1, lati(2)-1, 1, 1 /) )  ! meta time ind should always be one. Constant field so doesn't matter
			call check (nf90_get_var(meta_ncid, tId, define_levels, inds, (/ 4,4,dimlen,1/) ) )

			inds = int( (/loni(2)-1, lati(2)-1, 1, ti(2) /) )
			call check( nf90_inq_varid(ncid, "PH", tId ) )
			call check (nf90_get_var(ncid, tId, pert_levels, inds, (/ 4,4,dimlen,1/) ) )

		   define_levels = (define_levels + pert_levels) / g
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


	real function get4dval(ncid, reqvar, in_time, in_lev, in_lat, in_lon)
	real ti(3), levi(3), loni(3), lati(3),  in_time, in_lat, in_lon, in_lev
	real pro1(4), var(4,4,4)
	integer inds(4), ninds(4), ncid
	type (wind) reqvar
	real bicubic_interp, neville_interp
	ninds = (/4,4,4,1/)

    call coord_2_int(time, in_time, ti)
    call coord_2_int(reqvar%lev, in_lev,levi)
    call coord_2_int(reqvar%lon, in_lon,lati)
    call coord_2_int(reqvar%lat, in_lat,loni)
    if(.not.istat)then
	return
    end if
    inds = int( (/loni(2)-1, lati(2)-1, levi(2)-1, ti(2) /) )

    call check (nf90_get_var(ncid, reqvar%id, var, inds, ninds ) )

	!  read in bottom left - bottom right - top left -  topright
	do i=1,4
	pro1(i) = bicubic_interp(var(:,:,i), lati(1), loni(1) )
	end do
	val= neville_interp(reqvar%lev(levi(2)-1:levi(2)+2), pro1, start_lev, size(pro1) )


	if(ti(1).ne.0)then
!   would also be a istat time.
	end if

    get4dval = val
    end function get4dval



    subroutine coord_2_int(x, xx, coords)
    real x(:), xx, int, coords(3)
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
!    if(.and. ((i .le. 2) .or. (i .ge. size(x)-3)) )then
!       istat = istat+1
! add check to inds in get levels maybe instead.
! can then check if dimension concurrently.
!   end if
    coords(1) = int
    coords(2) = i
    coords(3) = ( x(i+1) - x(i) )
    end subroutine coord_2_int






	end PROGRAM
