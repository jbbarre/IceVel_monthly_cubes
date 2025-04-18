pro merge_unw_ncdf_test,i,j 

path='/mnt/pennell-z0/eric/ANTARCTICA/'

geo0=load_geo_param()

result=unw2vel_feather_sp_fortran(i,j,geo0,/verbose,/debug)


end
