;###############################################################
; make_grid_fil.gfs_0p5.ncl
; Matthew Janiga (janiga@atmos.albany.edu
; Last Update: Oct. 24, 2011
;
; Description: Build a grid file for input to the trajectory model
; using 0.5 degree analyses from the GFS.
;###############################################################


;;---- NCL Libraries ----;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"


;;---------------------------------------------------------------------------;;
;; Build the wind file

undef("create_wind_n_aux_file")
procedure create_wind_n_aux_file(start_date:double, end_date:double, opt:logical)
begin

  ;----------------------------------------------------
  ; Program settings

  ; These could be commented out and specified as input

  opt@u_name = "UGRD_P0_L100_GLL0"
  opt@v_name = "VGRD_P0_L100_GLL0"
  opt@w_name = "VVEL_P0_L100_GLL0"
  opt@sp_name = "PRES_P0_L1_GLL0"

  opt@t_name = "TMP_P0_L100_GLL0"
  opt@rh_name = "RH_P0_L100_GLL0"
  opt@cape_name = "CAPE_P0_L1_GLL0"

  opt@p_name = "lv_ISBL0"
  opt@lat_name = "lat_0"
  opt@lon_name = "lon_0"

  opt@lev_max = 100000.0
  opt@lev_min = 10000.0

  ; Based on the opt@in_dir determine files in the directory

  file_list = opt@in_dir+systemfunc("/bin/ls -1 "+opt@in_dir)+".grib"
  nfile = dimsizes(file_list)

  ;----------------------------------------------------
  ; Initialize the arrays that will be written to disk
  
  ; Metadata

  data_p = addfile(file_list(0), "r") 

  latu = (/data_p->$opt@lat_name(::-1)$/)
  lat = latu(::-1)
  lat@units = "degrees_north"
  lat@long_name = "latitude"
  nlat = dimsizes(lat)

  lon = fspan(-180.0, 179.5, 720)
  lon@units = "degrees_east"
  lon@long_name = "longitude"
  nlon = dimsizes(lon)

  levf = (/data_p->$opt@p_name(::-1)$/)
  levf = levf(::-1)
  lev = levf(0:20)
  lev@units = "Pa"
  lev@long_name = "pressure"
  nlev = dimsizes(lev)


  time = start_date+fspan(0, end_date-start_date, toint((end_date-start_date)/(6.0*3600.0))+1)
  time@units = opt@time_units
  ntime = dimsizes(time)


  ; Wind field

  u = new((/ntime, nlev, nlat, nlon/), "float")
  u@units = "m s^-1"
  u@long_name = "zonal_wind"

  v = new((/ntime, nlev, nlat, nlon/), "float")
  v@units = "m s^-1"
  v@long_name = "meridional_wind"

  w = new((/ntime, nlev, nlat, nlon/), "float")
  w@units = "Pa s^-1"
  w@long_name = "pressure_vertical_velocity"

  sp = new((/ntime, nlat, nlon/), "float")
  sp@units = "Pa"
  sp@long_name = "surface_pressure"
  
  u!0 = "time"
  u!1 = "lev"
  u!2 = "lat"
  u!3 = "lon"

  u&time = time
  u&lev = lev
  u&lat = lat
  u&lon = lon

  copy_VarMeta(u, v)
  copy_VarMeta(u, w)

  sp!0 = "time"
  sp!1 = "lat"
  sp!2 = "lon"

  sp&time = time
  sp&lat = lat
  sp&lon = lon

  ; Auxiliary variables

  th = new((/ntime, nlev, nlat, nlon/), "float")
  th@units = "K"
  th@long_name = "potential_temperature"

  q = new((/ntime, nlev, nlat, nlon/), "float")
  q@units = "g kg^-1"
  q@long_name = "specific_humidity"

  cape = new((/ntime, nlat, nlon/), "float")
  cape@units = "J kg^-1"
  cape@long_name = "convective_available_potential_energy"
  
  copy_VarMeta(u, th)
  copy_VarMeta(u, q)
  copy_VarMeta(sp, cape)

  ;----------------------------------------------------
  ; Fill the grid file arrays

  ; NOTE: Flip to center on Prime Meridian, have lat increase, and have pressure increase

  do i=0,nfile-1

    print("Reading GRIB file "+i)
    data_p = addfile(file_list(i), "r") 

    ; Required variables

    u(i,:,:,:) = (/ lonFlip( data_p->$opt@u_name$({lev}, ::-1, :) ) /)
    v(i,:,:,:) = (/ lonFlip( data_p->$opt@v_name$({lev}, ::-1, :) ) /)
    w(i,:,:,:) = (/ lonFlip( data_p->$opt@w_name$({lev}, ::-1, :) ) /)

    sp(i,:,:) = (/ lonFlip( data_p->$opt@sp_name$(::-1, :) ) /)


    ; Auxiliary variables
    
    p = conform_dims(dimsizes(u(i,:,:,:)), lev, 0)

    t = lonFlip( data_p->$opt@t_name$({lev}, ::-1, :) )
    th(i,:,:,:) = (/ t*(100000.0/p)^(287.0/1004.0) /)

    rh = lonFlip( data_p->$opt@t_name$({lev}, ::-1, :) )
    q(i,:,:,:) = (/ mixhum_ptrh(p, t, rh, -2) /)

    cape(i,:,:) = (/ lonFlip( data_p->$opt@cape_name$(::-1, :) ) /)

  end do


  ;----------------------------------------------------
  ; Write to disk wind file

  ; File and Variable Attributes

  fAtt = True

  fAtt@creation_date  = systemfunc ("date")        
  fAtt@source         = "/nmc2/grib/"
  fAtt@comments       = "Matthew Janiga (matthew.janiga@gmail.com)" 

  ; File pointer and attributes

  system("/bin/rm -f "+opt@wind_file)  
  out_file_p = addfile(opt@wind_file, "c" )
 
  fileattdef( out_file_p, fAtt )            ; Set file attributes   

  ; Write arrays to record

  out_file_p->u = u
  out_file_p->v = v
  out_file_p->w = w
  out_file_p->sp = sp

  print("wind_file netCDF file initialized on "+systemfunc("date"))


  ;----------------------------------------------------
  ; Write to disk aux file

  ; File and Variable Attributes

  fAtt = True

  fAtt@creation_date  = systemfunc ("date")        
  fAtt@source         = "/nmc2/grib/"
  fAtt@comments       = "Matthew Janiga (matthew.janiga@gmail.com)" 

  ; File pointer and attributes

  system("/bin/rm -f "+opt@aux_file)  
  out_file_p = addfile(opt@aux_file, "c" )
 
  fileattdef( out_file_p, fAtt )            ; Set file attributes   

  ; Write arrays to record

  out_file_p->th = th
  out_file_p->q = q
  out_file_p->cape = cape

  print("wind_file netCDF file initialized on "+systemfunc("date"))


end


;;---------------------------------------------------------------------------;;
;; Create the trajectory file

undef("create_traj_file")
procedure create_traj_file(opt:logical)
begin

  ; Specify domain and resolution of grid

  nlat = floattoint( (opt@grid_lat_max-opt@grid_lat_min)/opt@grid_res + 1 )
  nlon = floattoint( (opt@grid_lon_max-opt@grid_lon_min)/opt@grid_res + 1 )


  ; Create a matrix of trajectories

  ntraj = nlat*nlon
  traj = ispan(1,ntraj,1)
  traj@long_name = "trajectory_index"

  lat_box = new((/nlat,nlon/), "float")
  do ln=0,nlon-1
    lat_box(:,ln) = fspan(opt@grid_lat_min,opt@grid_lat_max,nlat)
  end do

  lon_box = new((/nlat,nlon/), "float")
  do la=0,nlat-1
    lon_box(la,:) = fspan(opt@grid_lon_min,opt@grid_lon_max, nlon)
  end do


 ; Apply metadata

  lat = ndtooned(lat_box)
  lat!0 = "traj"
  lat&traj = traj
  lat@units = "degrees_north"
  lat@long_name = "latitude"

  lon = ndtooned(lon_box)
  lon!0 = "traj"
  lon&traj = traj
  lon@units = "degrees_east"
  lon@long_name = "longitude"

  time = conform_dims(dimsizes(lat),opt@init_time,-1)
  time@units = opt@time_units
  time!0 = "traj"
  time&traj = traj

  lev = conform_dims(dimsizes(lat),opt@init_pres,-1)
  lev@units = "Pa"
  lev!0 = "traj"
  lev&traj = traj


  ;----------------------------------------------------
  ; Write to disk

 ; File and variable attributes

  fAtt = True

  fAtt@creation_date  = systemfunc ("date")        
  fAtt@comments       = "Matthew Janiga (matthew.janiga@gmail.com)" 
  fAtt@description    = "Grid of initial points for simtraj" 

  ; File pointer and attributes

  system("/bin/rm -f "+opt@traj_file)  
  out_fil_p = addfile(opt@traj_file, "c" )
 
  fileattdef( out_fil_p, fAtt )            ; Set file attributes   

  ; Write coordinates to record

  out_fil_p->time = time
  out_fil_p->lev = lev
  out_fil_p->lat = lat
  out_fil_p->lon = lon

  print("traj_file netCDF file initialized on "+systemfunc("date"))

end


;;---------------------------------------------------------------------------;;
;; Create the namelist file

undef("create_nml_file")
procedure create_nml_file(opt:logical)
begin

  nml_s = new((/11/), "string")

  nml_s(0) = "&PARAMS"
  nml_s(1) = "gridfile= '"+opt@wind_file+"'"
  nml_s(2) = "auxfile = '"+opt@aux_file+"'"
  nml_s(3) = "trajfile = '"+opt@traj_file+"'"
  nml_s(4) = "outfile = '"+opt@output_file+"'"
  nml_s(5) = "naux = 3"
  nml_s(6) = "aux_dim_t = 4,4,3"
  nml_s(7) = "aux_name_t = 'th','q','cape'"
  nml_s(8) = "nstep = 73"
  nml_s(9) = "time_step = -3600"
  nml_s(10) = "/"

  ; Write namelist to disk

  asciiwrite(opt@nml_file, nml_s)

  print("namelist file initialized on "+systemfunc("date"))

end


;;---------------------------------------------------------------------------;;
;; Call the executable

undef("run_simtraj")
procedure run_simtraj(opt:logical)
begin

;  print("execution started on "+systemfunc("date"))

  sys_cmd = opt@exc+" "+opt@nml_file

 ; print(sys_cmd)
  system(sys_cmd)

;  print("execution completed on "+systemfunc("date"))

end


;;---------------------------------------------------------------------------;;
;; Main()

begin

  ;---------------------
  ; Program settings

  out_dir = "/ct13/janiga/programs/simtrajectory/test_data/runs/gfs_0p5/20111030_300.lint_polyp/"

  opt = True

  ; Input grib files

  opt@in_dir = "/ct13/janiga/programs/simtrajectory/test_data/gfs_0p5/grib/"

  ; Files used by the trajectory program

  opt@wind_file = "/ct13/janiga/programs/simtrajectory/test_data/gfs_0p5/nc/wind_grid.nc"
  opt@aux_file = "/ct13/janiga/programs/simtrajectory/test_data/gfs_0p5/nc/aux_grid.nc"

  opt@traj_file = out_dir+ "traj_input.nc"
  opt@nml_file = out_dir + "exc_input.warm_sector.back.60min.3pet.nml"

  ; Time

  opt@time_units = "seconds since 2010-10-26 00:00:00"
  start_date = ut_inv_calendar(2010,10,26,00,0,0,opt@time_units, 0)
  end_date = ut_inv_calendar(2010,10,31,18,0,0,opt@time_units, 0)

  ; Initial conditions of the trajectories

  opt@init_time = ut_inv_calendar(2010,10,30,0,0,0,"seconds since 2010-10-26 00:00:00", 0)

  opt@init_pres = 30000.0

  opt@grid_lat_min = 10.0
  opt@grid_lat_max = 50.0
  opt@grid_lon_min = -70.0
  opt@grid_lon_max = -10.0
  opt@grid_res = 1.0

  ; Output files

  opt@output_file = out_dir +"warm_sector.back.60min.3pet.nc"

  ; Executable

  opt@exc = "/ct13/janiga/programs/simtrajectory/exc/simtraj_lint_polyp"


  ;-----------------------------------------------
  ; Create input files for the trajectory program

;  create_wind_n_aux_file(start_date, end_date, opt)
  create_traj_file(opt)
  create_nml_file(opt)


  ;-----------------------------------------------
  ; Run the program
  
  print("execution started on "+systemfunc("date"))

  run_simtraj(opt)

  print("execution completed on "+systemfunc("date"))

  ;-----------------------------------------------
  ; Plot the results

  
end