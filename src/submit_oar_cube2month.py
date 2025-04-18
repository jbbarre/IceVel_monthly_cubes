
from cube_class import cube_class
from glob import glob
from cube2monthly_map_using_linear_regression import cube2monthly_map_using_linear_regression, check_vy_for_s2
from cube_calib_using_low_speed import cube_calib_using_low_speed
import os, sys
import time 

def load_cubes(file):
    c=cube_class()
    print('Start loading cubes')

    rootdir = '/summer/ice_speed/surface_flow_velocity/'

    files = glob(rootdir+'ANTARCTICA/SENTINEL2/2018/10d/MOSAIC/cubes/'+file)
    c.load(files[0])
    
    print('Loading RADARSAT-2....')
    for yr in range(2016,2022):
        if glob(rootdir+'ANTARCTICA/RADARSAT-2_'+str(yr)+'/MOSAIC/cubes/'+file):
            l = glob(rootdir+'ANTARCTICA/RADARSAT-2_'+str(yr)+'/MOSAIC/cubes/'+file)
            print(l[0])
            c.load(l[0],concat=True)
    
    print('Loading LANDSAT....')
    for yr in range(2013,2022):
        for cy in range(16,416,16):
            if glob(rootdir+'ANTARCTICA/LANDSAT/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                for l in glob(rootdir+'ANTARCTICA/LANDSAT/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                    print(l)
                    c.load(l,concat=True)    
    
    print('Loading S1....')
    for yr in range(2014,2022):
        for cy in range(6,30,6):
            if glob(rootdir+'ANTARCTICA/SENTINEL1/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                for l in glob(rootdir+'ANTARCTICA/SENTINEL1/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                    print(l)
                    c.load(l,concat=True)   

    print('Loading S2....')
    for yr in range(2016,2022):
        for cy in range(5,405,5):
            if glob(rootdir+'ANTARCTICA/SENTINEL2/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                for l in glob(rootdir+'ANTARCTICA/SENTINEL2/'+str(yr)+'/'+str(cy)+'d/MOSAIC/cubes/'+file):
                    print(l)
                    c.load(l,concat=True)   


    c.describe()
    c.sieve_empty_layers(sieve_n=1000)
    c.describe() 

    return c               
#status, vxo, vyo, stdxo, stdyo, nno, medx0, medy0, stdxo0, stdyo0, errxo, erryo, nno0 = c.cube2vel()

def submit_oar_cube2monthly(i,j):
    
    print('Starting the process')
    start_time= time.time()
    nam = 'x{0:05d}_y{1:05d}'.format(int(i),int(j))

    if len(glob('cube_monthly_'+nam+'*.nc')) != len(range(2015,2022)):

        file = 'c_'+nam+'_post0450_yr*nc'

        c=load_cubes(file)

        print('check vy for S2')
        c=check_vy_for_s2(c)
        
        print ('Start calibration...')
        c2 = cube_calib_using_low_speed(c)
        c=None
        for i in range(2015,2022):
            if not os.path.exists('cube_monthly_'+nam+'_yr'+str(i)+'.nc'):
                ds = cube2monthly_map_using_linear_regression(c2,i,i,deltamax=36) #deltamax = cycles in days

                comp = dict(zlib=True, complevel=5)
                encoding = {var: comp for var in ds.data_vars}

                ds.to_netcdf('cube_monthly_'+nam+'_yr'+str(i)+'.nc',encoding=encoding)
                ds.close()
    print(f'process cube in {time.time()-start_time} s.')
if __name__ == "__main__":

    #submit_oar_cube2monthly(sys.argv[1],sys.argv[2])
    submit_oar_cube2monthly('02695','07105')
