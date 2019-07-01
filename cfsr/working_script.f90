    PROGRAM cfsr_trajectories
    USE netcdf
!    IMPLICIT NONE   
    character (len = 999), dimension(:), allocatable :: VAR_FILE
    character (len = 499) :: windfiles(5), metafile
    character (len = 499) :: filename, filedir, filesuf    
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
    	real, dimension(:), allocatable :: time
    	real, dimension(:,:,:,:), allocatable :: grid
    	character (len = 25):: units, name
    	integer :: id, ti, fid
    end type wind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Define types each separate trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type parcel
    	real :: u1,v1,w1,u2,v2,w2, lev, lat, lon, time
    	logical :: bound = .true.
    end type parcel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! End of Definitions lets initialise stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (wind) :: u,v,w,p
    type (parcel), dimension (:), allocatable :: traj



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialise all grids 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    filedir = "/cfsr/data/2012/"
    filesuf = ".2012.0p5.anl.nc"
    
    u%name = "u"
    v%name = "v"
    w%name = "w"
    p%name = "psfc"

    call open_file(u)
    call open_file(v)
    call open_file(w)
    call open_file(p)

    
    call define_coords(u)
    call define_coords(v)
    call define_coords(w)
    call define_coords(p)

    start_time = 1858944
    
    call coord_2_int(u%time, start_time, ti)
    call grab_grid(u, int(ti(2)))
    call grab_grid(v, int(ti(2)))
    call grab_grid(w, int(ti(2)))
    call grab_grid(p, int(ti(2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Read the namelist file. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(1,file=namefile)
	read(1,fileio)  
	read(1,time_opt)
	read(1,traj_opt)
	



    	contains

     subroutine open_file(vari)
    	type (wind) vari
        character (len=99) fname
        fname =  trim(filedir)//trim(vari%name)//trim(filesuf)
        print*, trim(fname)
        call check( nf90_open(fname, NF90_NOWRITE, vari%fid) )
        return
        end subroutine


	subroutine define_coords(vari)
	type (wind) vari
	integer varId, tId, dimids(6), dimlen
	character (len =25) dimname, newname
	real, dimension(:), allocatable :: lonvar, latvar, levvar, timevar

    call check( nf90_inq_varid(vari%fid, vari%name, vari%id ) )
    call check( nf90_inquire_variable(vari%fid, vari%id, ndims = ndim, dimids = dimids, natts = numAtts) )
    do i=1,ndim
		call check( nf90_inquire_dimension(vari%fid, dimIds(i),name = dimname, len=dimlen) )
		select case (dimname)
		case( "lon")
			allocate ( vari%lon(dimlen) )
			allocate ( lonvar(dimlen) )
			call check( nf90_inq_varid(vari%fid, "lon", tId) )
			call check (nf90_get_var(vari%fid, tId, lonvar, (/1/), (/dimlen/) ) )
			vari%lon = lonvar
		case( "lat")
			allocate ( vari%lat(dimlen) )
			allocate ( latvar(dimlen) )
			call check( nf90_inq_varid(vari%fid, "lat", tId ) )
			call check (nf90_get_var(vari%fid, tId,latvar, (/1/), (/ dimlen/) ) )
			vari%lat = latvar
		case( "lev")
            allocate (vari%lev(dimlen))
            allocate ( levvar(dimlen) )
            call check( nf90_inq_varid(vari%fid, "lev", tId ) )
            call check (nf90_get_var(vari%fid, tId,levvar, (/1/), (/ dimlen/) ) )
            vari%lev = levvar
		case( "time")
            allocate (vari%time(dimlen))
            allocate ( timevar(dimlen) )
            call check( nf90_inq_varid(vari%fid, "time", tId ) )
            call check (nf90_get_var(vari%fid, tId,timevar, (/1/), (/ dimlen/) ) )
            vari%time = timevar
		end select
    end do
    
   if(vari%name.eq."psfc")then
      allocate (vari%lev(1))  
      vari%lev = 1  
   end if

    allocate(vari%grid( size(vari%lon), size(vari%lat), size(vari%lev),2 ) )
    
	print*, vari%name
	print*, "Longitude  (", size(vari%lon),"): ",minval(vari%lon),":",maxval(vari%lon)
	print*, "Latitude   (", size(vari%lat),"): ",minval(vari%lat),":",maxval(vari%lat)
	print*, "Levels   (", size(vari%lev),"): ",minval(vari%lev),":",maxval(vari%lev)
	print*, "Time   (", size(vari%time),"): ",minval(vari%time),":",maxval(vari%time)


	return
	end subroutine define_coords



    subroutine grab_grid(reqvar, ti)
    type(wind) reqvar
    integer ti, fti, inc,i,lncid
    real fi(3), pfi(3)
        do i=0,1
         call coord_2_int(reqvar%time,reqvar%time(ti+i),  pfi)
         fti = pfi(2)
         call check( nf90_inq_varid(reqvar%fid, reqvar%name, reqvar%id ) )
         call check( nf90_get_var(reqvar%fid, reqvar%id, reqvar%grid(:,:,:,1+i), (/1,1,1,fti/), (/ size(reqvar%lon),size(reqvar%lat),size(reqvar%lev),1/) ) )
        end do
        reqvar%ti = ti
!     	print*, reqvar%name,  minval(reqvar%grid),":",maxval(reqvar%grid)
    end subroutine


    subroutine coord_2_int(x, xx, coords)
    real x(:), xx, int, coords(3)
    if(xx.lt.minval(x) .or. xx.gt.maxval(x))then
    print*, "out of bounds"
    end if
	 do i=1, size(x)-1
		if(x(2) .gt. x(1)) then
			 if( xx >= x(i) .and. xx < x(i+1) )then
			 int = 1.0* (xx - x(i)) / (x(i+1) - x(i) )
			 exit
			 end if
		 else
			 if( xx <= x(i) .and. xx > x(i+1) )then
			 int = 1.0* (xx - x(i)) / (x(i+1) - x(i) )
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
