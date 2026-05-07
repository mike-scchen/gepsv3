load "/nwpr/gfs/xb157/include/dmskey_io.ncl"
load "/nwpr/gfs/xb157/include/mjo_diagnostic_coda.ncl"

  setfileoption("nc", "Format", "LargeFile")
  case=""
  root_path = "/nwpr/gfs/xb157/GPU/GEPSv3/jobs/18091800_000.ufs"+case+"/TCo383"
  case_name = "gfs"
  outfile = case+"_atmo.nc"

  init_yr = "2018"
  init_mm = "09"
  init_dd = "18"
  init_hr = "00"
  init_mn = "00"
  
  tau_ini = 0
  tau_gap = 6
  tau_end = 0
  dms_list = (/"500000","850100","850200","850210","200200","SSL010"/)
  var_name = (/"HGT500", "850t","U850","V850","U200","V200"/)
  var_unit = (/"m", "C","m/s","m/s","m/s","m/s"/)
  var_scrp = (/"HGT500",\
               "T850",\
               "U850",\
               "V850",\
               "U200",\
               "V200"/)
  
  ggdef = "GI0GH1191936"

  init_date = init_yr + init_mm + init_dd
  init_datetime = init_date + init_hr + init_mn
  in_path = root_path ; + "/" + init_date + "/" + case_name

  TCo_file = "/nwpr/gfs/xb157/diag/TCo_grid.nc"
  TCo = addfile(TCo_file, "r")
  lat = TCo->lat
  lon = TCo->lon

  ia = ind(lat .ge. -90.0d0  .and. lat .le. 90.0d0)
  io = ind(lon .ge. 0.0d0 .and. lon .le. 360.0d0)
;  ia = ind(lat .ge. 15.0d0  .and. lat .le. 35.0d0)
;  io = ind(lon .ge. 110.0d0 .and. lon .le. 140.0d0)

  latp = lat(ia)
  lonp = lon(io)

  print(latp)

;  lat = asciiread("../grid_dims/t319_lat.txt" , -1, "double")
;  lon = asciiread("../grid_dims/t319_lon.txt" , -1, "double")
  forecast_tau = ispan(tau_ini, tau_end, tau_gap)
  delete(lat@_FillValue)
  delete(lon@_FillValue)

  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  ntime = dimsizes(forecast_tau)
  ndims = (/ntime, nlat, nlon/)
  ndmskey = dimsizes(dms_list)

  nlatp = dimsizes(latp)
  nlonp = dimsizes(lonp)
  ndimsp = (/ntime, nlatp, nlonp/)
  dimName = dmskey2nc_create2D(outfile, init_datetime, forecast_tau, latp, lonp, ndimsp)

  do n = 0, ndmskey-1
    print(dms_list(n))
    dmskey_ggdef = dms_list(n) + ggdef
    time2D_value = read_dmskey_time2D(in_path, init_datetime, dmskey_ggdef, forecast_tau, ndims)
    time2D_value@short_name = var_name(n)
    time2D_value@long_name = var_scrp(n)
    time2D_value@unit = var_unit(n)
    time2D_value@dmskey = dms_list(n)

   ; avg_zonal = mean_zonal(time2D_value, lat, (/-10.0,10.0/), ndims)
   ; avg_zonal@short_name = "avg_"+var_name(n)
   ; avg_zonal@long_name = "zonal mean of " + var_scrp(n)
   ; avg_zonal@unit = var_unit(n)
   ; avg_zonal@dmskey = dms_list(n)

   ; if( n .eq. 0)
   ;   avg_zonal = avg_zonal/2.5e6*86400.0
   ; end if
    ierr = dmskey2nc_putvar2D(outfile, time2D_value(:,ia,io), dimName)
   ; ierr = dmskey2nc_putvar2D(outfile, avg_zonal, (/dimName(0),dimName(2)/))
    delete(time2D_value)
   ; delete(avg_zonal)
    print(n)
  end do 

