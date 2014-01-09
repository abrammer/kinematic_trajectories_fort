    PROGRAM wrf_trajectories
    USE netcdf
    IMPLICIT NONE
    character (len = 999), dimension(:), allocatable :: VAR_FILE
    
    character (len = 499) :: windfiles(5), metafile
    character (len =999) ::  META_FILE 
    character (len = *), parameter :: namefile = "namelist.traj",  FMT1 = "(7F15.5)", FMT0="(7A15)"
    character (len = 25)  lat_name, lon_name, lev_name, tim_name, var_name
    real, dimension(:), allocatable :: lat, lon, lev, time, file_times
    real ti(3), li(3), loni(3),lati(3), var(4,4,4), pro1(4), val
    real start_time, start_lon, start_lat, start_lev, end_time
    real t2, t1,  filetimes(99), outmin, dx,dy,dl
    integer, parameter :: text_out = 20
    integer  meta_ncid, varId, dimId, ndim, dimlen, uId, inds(4), ninds(4),i, numAtts, it,x,t, step, no_of_parcels
    integer nfiles,f,  ntime,mintime, maxtime, timedt,ny,nx,nl,sz,sy,sl
    integer, dimension(:), allocatable :: ncid
    integer, dimension(nf90_max_var_dims) :: dimIds
        
    namelist /fileio/metafile,nfiles, windfiles,filetimes
	namelist /time_opt/start_time,end_time,step,outmin, mintime,maxtime,timedt
	namelist /traj_opt/no_of_parcels,start_lev,start_lat,start_lon, dy, dx, dl, ny,nx,nl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! External Functions
    real bicubic_interp, neville_interp
    integer istat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Lets get Rolling.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Define types for 4d grids, containing their own metadata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type wind
    	real :: val1, val2, val
    	real, dimension(:), allocatable :: lat
    	real, dimension(:), allocatable :: lon
    	real, dimension(:), allocatable :: lev
    	real, dimension(:,:,:,:), allocatable :: grid
    	real, dimension(:,:,:), allocatable :: base
    	character (len = 25):: units, name
    	integer :: id, y_stag, ti
    end type wind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Define types each separate trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type parcel
    	real :: u1,v1,w1,u2,v2,w2, lev, lat, lon, time
    	logical :: bound = .true.
    end type parcel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Initialise 4d grids and trajectory variables. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (wind) :: u,v,w,ph
    type (parcel), dimension (:), allocatable :: traj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Read the namelist file. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(1,file=namefile)
	read(1,fileio)  
	read(1,time_opt)
	read(1,traj_opt)
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Run through some initial file io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	META_FILE =TRIM(ADJUSTL(METAFILE))
    meta_ncid = open_file(META_FILE)

	allocate(file_times(nfiles))
	file_times = filetimes(:nfiles)

	allocate (VAR_FILE(nfiles))
	allocate(ncid(nfiles))
	do f=1,nfiles
	  VAR_FILE(f) = adjustl(windfiles(f))
	  ncid(f) = open_file(var_file(f))
	end do



    !var_name = "XTIME"
    !call check( nf90_inq_varid(ncid, var_name, varId ) )
    !call check( nf90_inquire_variable(ncid, varId, ndims = ndim, dimids = dimIds, natts = numAtts) )
    !call check( nf90_inquire_variable(ncid, varId, dimids = dimIds(:ndim) ) )
    !call check( nf90_inquire_dimension(ncid, dimids(1), len=dimlen) )
    !allocate ( time(dimlen) )
    !call check (nf90_get_var(ncid, varId, time) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Run through coordinate definition info. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Messy time coordinate definitions. 
     ntime = ( (maxtime-mintime)/timedt )+ 1
     print*, ntime
     allocate(time(ntime))
     do t=0, ntime-1
        time(t+1) = mintime+ (t*timedt)
     end do
        print*, time(ntime)
        
        
!    ***********************************
!    define dimensions
print*, "*****************************"
	u%name = "U"
	call define_coords(u)

	v%name = "V"
	call define_coords(v)

	w%name = "W"
	call define_coords(w)

	ph%name = "PH"
	call define_coords(ph)


!    ***********************************
!   starting point

!    allocate( traj(no_of_parcels) )
!    do t=1,no_of_parcels
!       traj(t)%lat = start_lat
!       traj(t)%lev = start_lev+(t*400)
!       traj(t)%lon = start_lon
!    end do
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  Initialise trajectories. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	print*, no_of_parcels, nx*ny*nl
    allocate( traj(no_of_parcels) )
    t=1
    do sx=1, nx
	  do sy=1, ny
		do sl=1, nl
		   traj(t)%lat = start_lat+(sy*dy)
		   traj(t)%lev = start_lev+(sl*dl)
		   traj(t)%lon = start_lon+(sx*dx)
		   t=t+1
		end do
	  end do
    end do


    print*, "Starting", no_of_parcels," trajectories"

!!    grab first time segment of data
    call coord_2_int(time, start_time, ti)
    call grab_grid(u, int(ti(2)))
    call grab_grid(v, int(ti(2)))
    call grab_grid(w, int(ti(2)))
    call grab_grid(ph, int(ti(2)))


    print*, "*********************"
    call integrate_trajectory(traj, start_time, end_time)


!***********
!***********
!***********
! All functions  below 
!***********
!***********
!***********

	contains

    integer function open_file(fname)
        character (len=*) fname
        call check( nf90_open(fname, NF90_NOWRITE, open_file) )
        end function


    subroutine grab_grid(reqvar, ti)
    type(wind) reqvar
    integer ti, fti, inc,i,lncid
    real fi(3), pfi(3)
        do i=0,1
         call coord_2_int(file_times, time(ti+i), fi)
         call coord_2_int(time, file_times(fi(2)),  pfi)
         fti = (ti+i)-(pfi(2)-1)

         lncid = ncid(fi(2))
         call check( nf90_inq_varid(lncid, reqvar%name, reqvar%id ) )
         call check( nf90_get_var(lncid, reqvar%id, reqvar%grid(:,:,:,1+i), (/1,1,1,fti/), (/ size(reqvar%lon),size(reqvar%lat),size(reqvar%lev),1/) ) )
        end do
        reqvar%ti = ti
    end subroutine


	subroutine define_coords(vari)
	type (wind) vari
	integer varId, tId, dimids(6)
	character (len =25) dimname, newname
	real, dimension(:), allocatable :: lonvar, latvar

    call check( nf90_inq_varid(ncid(1), vari%name, vari%id ) )
    call check( nf90_inquire_variable(ncid(1), vari%id, ndims = ndim, dimids = dimids, natts = numAtts) )
    do i=1,ndim
		call check( nf90_inquire_dimension(ncid(1), dimIds(i),name = dimname, len=dimlen) )
		select case (dimname)
		case( "west_east")
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLONG", tId) )
			call check (nf90_get_var(meta_ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		case( "west_east_stag")
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLONG_U", tId ) )
			call check (nf90_get_var(meta_ncid, tId, lonvar, (/1,1,1,1/), (/dimlen,1,1,0/) ) )
			vari%lon = lonvar
		case( "south_north")
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLAT", tId ) )
			call check (nf90_get_var(meta_ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
			vari%lat = latvar
		case( "south_north_stag")
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(meta_ncid, "XLAT_V", tId ) )
			call check (nf90_get_var(meta_ncid, tId,latvar, (/1,1,1,1/), (/ 1,dimlen,1,0/) ) )
			vari%lat = latvar
		case( "bottom_top")
            allocate (vari%lev(dimlen))
            vari%y_stag = 0
		case( "bottom_top_stag")
            vari%y_stag = 1
            allocate (vari%lev(dimlen))
		end select
    end do
    
    allocate(vari%grid( size(vari%lon), size(vari%lat), size(vari%lev),2 ) )
    if(vari%name =="PH")then
        allocate(vari%base( size(vari%lon), size(vari%lat), size(vari%lev) ) )
	    call check( nf90_inq_varid(meta_ncid, "PHB", tId ) )
		call check (nf90_get_var(meta_ncid, tId, vari%base, (/1,1,1,1/), (/ size(vari%lon), size(vari%lat), size(vari%lev),1/) ) )
    else
       allocate(vari%base(1,1,1)  )  !! just so we don't have unallocated grids? Not sure whether this matters
    end if
    
    print*, vari%name
	print*, "Longitude  (", size(vari%lon),"): ",minval(vari%lon),":",maxval(vari%lon)
	print*, "Latitude   (", size(vari%lat),"): ",minval(vari%lat),":",maxval(vari%lat)
	return
	end subroutine define_coords


    subroutine integrate_trajectory(traj, stime,end_time)
    type (parcel):: traj(:)
    real time, end_time, stime,minstep
    integer t, ob
    character (len=3) :: stri

	minstep = step/60D0  ! have an issue with rounding errors at the moment. 

    do t=1,no_of_parcels
	   call get5dval(u,v,w, stime , traj(t)%lev, traj(t)%lat, traj(t)%lon, traj(t) )
	   if(.not.traj(t)%bound)then
	   		ob = ob+1
		 cycle
	   end if
	   write(stri,'(I3)') t
	   open(unit=(100+t), file="traj"//trim(ADJUSTL(stri))//".txt", status="replace", action="write")
	   write((100+t),FMT0) "TIME","LEV","LAT","LON","U","V","W"
	   write((100+t),FMT1)  stime , traj(t)%lev, traj(t)%lat, traj(t)%lon, traj(t)%u1, traj(t)%v1, traj(t)%w1
    end do

	if(ob)then
	print*, ob,"Trajectories are not starting due to being out of bounds"
	end if

   open(unit=6, carriagecontrol='fortran')
    do time=stime, end_time, minstep
	   do t=1,no_of_parcels
			if(traj(t)%bound)then
				 call petterssen(u,v,w,  time ,traj(t))
				 if( mod( time, outmin).lt.minstep .and. traj(t)%bound)then
					 write((100+t),FMT1)  time , traj(t)%lev, traj(t)%lat, traj(t)%lon, traj(t)%u1,traj(t)%v1,traj(t)%w1
				 end if
			 end if
		 end do
   call progress( (time-stime)/(end_time-stime)*100. )
   end do

   end subroutine integrate_trajectory



    subroutine petterssen(u,v,w,in_time,traj)
    real, parameter :: converg = 0.000001   !!! check for convergence
    type (wind) :: u,v,w
    type (parcel) :: traj
    real locs(3), newlocs(3), lev_in, lat_in, lon_in, in_time
    real u0,v0,w0
    integer it
	! Save initial winds.
    u0 = traj%u1
    v0 = traj%v1
    w0 = traj%w1

	! Make an initial movement with initial winds.
    newlocs = move_parcel(traj%u1,traj%v1,traj%w1, traj%lev, traj%lat, traj%lon)
    call get5dval(u,v,w, in_time, newlocs(1), newlocs(2), newlocs(3), traj )  ! get initial guess
	! Use first guess location, and make new movement using average of winds from X0 and X1, iterate until they are the same
    do it=1, 3		! This rarely take 3 iterations.
	   locs  = newlocs
	   newlocs = move_parcel(0.5*(u0+traj%u1),0.5*(v0+traj%v1),0.5*(w0+traj%w1), traj%lev, traj%lat, traj%lon)
	   call get5dval(u,v,w, in_time, newlocs(1), newlocs(2), newlocs(3), traj )
!	   	if(.not.traj%bound)then
!		   return
!	   end if

	   if(maxval(abs(locs - newlocs)).le.converg)then 	! check for convergence.
		  exit
	   end if
	   if(it.eq.3)then
		   print*, maxval(abs(locs - newlocs)), "solution not converged : Likely Due to CFL Issue. Try reducing the timestep"
	   end if
    end do
    traj%lev = newlocs(1)
    traj%lat = newlocs(2)
    traj%lon = newlocs(3)
    end subroutine petterssen


    real function move_parcel(Vu,Vv,Vw,lev_in, lat_in, lon_in)
    dimension move_parcel(3)
    real Vu,Vv,Vw
    REAL, PARAMETER :: pi = 3.14159265358979  ! pi
    REAL, PARAMETER :: radius = 6371220.0     ! Earth radius (m)
    REAL :: latr0, lonr0, latr1, lonr1, lat_in, lon_in       ! Initial and guess lat/lon (radians)
    REAL :: dir, rdist, adj , newlat, newlon, time_step, newlev, lev_in          ! Direction, radial distance, backward/forward switch

    time_step = step
    adj = -1.

	rdist = SQRT( Vu**2.0+Vv**2.0 )*(adj*time_step/radius)      ! is adj needed?
    dir = pi+ATAN2( adj*-1.0*Vu, adj*-1.0*Vv )

    latr0 = torad(lat_in)
    lonr0 = torad(lon_in)

    ! Calculate new lat and lon using Haversine's forumula
    latr1 = ASIN( SIN(latr0)*COS(rdist)+COS(latr0)*SIN(rdist)*COS(dir) )
    lonr1 = lonr0+ATAN2( SIN(dir)*SIN(rdist)*COS(latr0), COS(rdist)-SIN(latr0)*SIN(latr1) )

    newlev = lev_in + (adj*time_step)*Vw

    newlat = todeg(latr1)  ! convert to degs
    newlon = todeg(lonr1)
    move_parcel = (/newlev, newlat, newlon/)
    end function move_parcel


	REAL function torad(val)
	REAL, PARAMETER :: pi = 3.14159265358979
	real val
    torad = val* pi/180.
    end function torad

    REAL function todeg(val)
	REAL, PARAMETER :: pi = 3.14159265358979
	real val
    todeg = val* 180./pi
    end function todeg


	subroutine get5dval(u,v,w, in_time, lev, lat, lon, traj)
	type (wind) u,v,w
    type (parcel) traj
	real in_time, lev, lat, lon, ti(3), t1, t2
	real linear_interp

	call coord_2_int(time, in_time, ti)
	if(.not.get_levels(u,v,w, time(ti(2)), lat, lon))then
	  traj%bound = .false.
	  return
	end if
	traj%u1 = get4dval(u, time(ti(2)), lev, lat, lon, traj%bound)
	traj%v1 = get4dval(v, time(ti(2)),lev, lat, lon, traj%bound)
	traj%w1 = get4dval(w, time(ti(2)), lev, lat, lon, traj%bound)
	if(ti(1).ne.0)then  ! if we're exactly on a data time, may as well take it. 
	  if(.not.get_levels(u,v,w, time(ti(2)+1), lat, lon))then
		traj%bound = .false.
		return
	  end if
	  traj%u2 = get4dval(u, time(ti(2)+1),lev, lat, lon, traj%bound)
	  traj%v2 = get4dval(v, time(ti(2)+1),lev, lat, lon, traj%bound)
	  traj%w2 = get4dval(w, time(ti(2)+1),lev, lat, lon, traj%bound)
	  traj%u1 = linear_interp( traj%u1, traj%u2, ti(1) )
	  traj%v1 = linear_interp( traj%v1, traj%v2, ti(1) )
	  traj%w1 = linear_interp( traj%w1, traj%w2, ti(1) )
	end if
	end subroutine get5dval



	logical function get_levels(u,v,w, in_time, in_lat, in_lon)
	type (wind) u,v,w
	real, parameter :: g=9.81
	real loni(3), lati(3),  in_lat, in_lon, in_time, t2, t1
	real :: define_levels(4,4,51)
	real :: profile(51), ustag_lev(50)
	integer inds(2), l, i, tt
	integer varId, tId, coords(2), ndim
	character (len =25) dimname, newname
	real bicubic_interp, neville_interp
	get_levels = .true.
	
     call coord_2_int(time, in_time, ti)
     call coord_2_int(ph%lat, in_lat, lati)
     call coord_2_int(ph%lon, in_lon,loni)

     inds = int( (/loni(2)-1, lati(2)-1/) )
	  if(any(inds.le.0) .or. (inds(1).ge.size(ph%lon)-2) .or. (inds(2).ge.size(ph%lat)-2) )then
		get_levels = .false.
		return
	  end if

     tt = ( ti(2) - ph%ti) + 1
     if(tt.gt.2)then
         call grab_grid(ph, ph%ti+1 )
		 tt = ( ti(2) - ph%ti) + 1
     end if

     define_levels = ( ph%base(inds(1):inds(1)+3, inds(2):inds(2)+3, :) + ph%grid(inds(1):inds(1)+3, inds(2):inds(2)+3, :, tt)) / g
     do l=1,size(profile)
         profile(l) = bicubic_interp(define_levels(:,:,l), lati(1), loni(1) )
     end do

     w%lev = profile
	 ustag_lev =  (profile(2:dimlen) + profile(1:dimlen-1) )/ 2.
	 u%lev = ustag_lev
	 v%lev = ustag_lev
    end function get_levels


	real function get4dval( reqvar, in_time, in_lev, in_lat, in_lon, bound)
	real ti(3), levi(3), loni(3), lati(3),in_time,in_lev, in_lat, in_lon
	real pro1(4), var(4,4,4), t1, t2
	integer inds(4), ninds(4), tt
	logical bound
	type (wind) reqvar
	real bicubic_interp, neville_interp

    call coord_2_int(time, in_time, ti)
    call coord_2_int(reqvar%lev, in_lev,levi)
    call coord_2_int(reqvar%lon, in_lon,lati)
    call coord_2_int(reqvar%lat, in_lat,loni)

    inds = int( (/loni(2)-1, lati(2)-1, levi(2)-1, ti(2) /) )
	  if(any(inds.le.0) .or. (inds(1).ge.size(reqvar%lon)-2) .or. (inds(2).ge.size(reqvar%lat)-2).or. (inds(3).ge.size(reqvar%lev)-2) )then
		bound = .false.
		return
	  end if
    
    tt = ( ti(2) - reqvar%ti) + 1
    if(tt.gt.2)then
        call grab_grid(reqvar, reqvar%ti+1 )
	    tt = ( ti(2) - reqvar%ti) + 1
    end if

    var = reqvar%grid(inds(1):inds(1)+3, inds(2):inds(2)+3, inds(3):inds(3)+3, tt)
	!  read in bottom left - bottom right - top left -  topright

	do i=1,4
	   pro1(i) = bicubic_interp(var(:,:,i), lati(1), loni(1) )
	end do
	val= neville_interp(reqvar%lev(inds(3):inds(3)+3), pro1, in_lev, size(pro1) )

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
    coords(1) = int
    coords(2) = i
    coords(3) = ( x(i+1) - x(i) )
    end subroutine coord_2_int






	end PROGRAM

     subroutine check(status)
     use netcdf
       integer, intent ( in) :: status
       if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
       end if
     end subroutine check


subroutine progress(j)
  implicit none
  integer k
  real j
  character(len=57)::bar="???% |"
  write(unit=bar(1:3),fmt="(i3)") int(j)
  do k=1, int(j/2)
    bar(6+k:6+k)="*"
  enddo
! print the progress bar.
  write(6,fmt="(a1,a1,a57)") '+',char(13), bar
  return
end subroutine progress
