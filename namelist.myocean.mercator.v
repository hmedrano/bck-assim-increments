!! ------------------- 
!! Namelist for SOSIE 
!! ------------------- 
!!
!!
!! *********************** 
!! Input characteristics : 
!! *********************** 
!! 
&ninput
ivect     = 0    ! this is not a vector interpolation
lregin    = T    ! input grid is regular (lon and lat are "1D")
cf_in     = 'INPUTFILE'
cv_in     = 'v'
cv_t_in   = 'time_counter'
jt1       = IDX    ! we want to interpolate each time record
jt2       = IDX    !           //
jplev     = 0
cf_x_in   = 'INPUTFILE'
cv_lon_in = 'longitude'
cv_lat_in = 'latitude'
cf_lsm_in = 'missing_value'   ! we use 'missing_value' of input field to determine
cv_lsm_in = ''                ! the land-sea-mask
ldrown    = T            ! we want to propagate sea values onto the land-sea mask
ewper     = 0          ! input field does have east-west periodicity with 0 overlapping point
vmax      =  99999
vmin      = -30000
/
!!
!! ***********************************
!!  IF 3D INTERPOLATION ( jplev = 0 )
!! ***********************************
!!
&n3d
cf_z_in  = 'INPUTFILE'
cv_z_in  = 'depth'
cf_z_out = 'MASKFILE'
cv_z_out = 'deptht'
cv_z_out_name = 'deptht'
/
!!
!!
!! ***************************** 
!! Output Grid characteristics : 
!! ***************************** 
!! 
!!
&noutput
lregout    = F
cf_x_out   = 'MASKFILE'
cv_lon_out = 'nav_lon'
cv_lat_out = 'nav_lat'
cf_lsm_out = 'MASKFILE'
cv_lsm_out = 'vmask'
lmout      = T
rmaskvalue = 0
lct        = F      ! we use time from input file
t0         = 0
t_stp      = 0
/
!! 
!! 
!! 
!! 
!! ******************************* 
!! Netcdf output characteristics : 
!! ******************************* 
!! 
&nnetcdf
cmethod  = 'akima'
cv_l_out = 'nav_lon'
cv_p_out = 'nav_lat'
cv_t_out = 'time_counter'
cv_out   = 'vomecrty'
cu_out   = 'm/s'
cu_t     = 'hours since 1950-01-01'
cd_out   = '.'
csource  = 'mercator'
ctarget  = 'nemomesh'
cextra   = 'interpol'
/ 
!!
