
#import matplotlib.pyplot as plt
from tracemalloc import stop
import numpy as np
import datetime
from progress_bar import progress_bar
from scipy.optimize import curve_fit
import xarray as xr
import time

def linear_speed_func(x, a, b, e):
    '''
    Return vx and vy values for a given speed magnitude and direction
    assume that flow direction remains the same
    x : datetime values
    a : slope (acceleration/decc trend)
    b : value at the origin (x=0)
    e : theta = flow direction
    '''
    vx = (a*x + b )*np.cos(e)
    vy = (a*x + b )*np.sin(e)
    
    f = np.array( [vx,vy] )
    return f.ravel()

def gaussian_func(x, amplitude, mean, stddev):
    '''
     gaussian function
    '''
    return amplitude * np.exp(-((x - mean)**2 / 2 / stddev**2))

def check_vy_for_s2(cube):
    
    sources = cube.source_()
    data2=[]
    
    for i, d in enumerate(cube.data):
        if 'SENTINEL1' in sources[i]:
            data2.append(d)
    
    cube2=cube.copy_structure()
    cube2.data=data2
    data2=None
    cube2.nz = len(cube2.data)
    status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = cube2.cube2vel()
    cube2=None
    
    for i, d in enumerate(cube.data):
        if 'SENTINEL-2' in sources[i] or 'SENTINEL2' in sources[i] :
            if np.nansum(np.abs(d.vy-vyo)) > np.nansum(np.abs(-d.vy-vyo)):
                    d.vy=-d.vy
    sources=None
    
    return cube

def linear_fit_speed(days,vxy,err=None, display=False,einit=None,filtering=True,nsigma=7):
    '''
    days : datetime list
    vxy : [[vx],[vy]] where vx and vy are lists containing speed in x and y direction
    param err : [[errx],[erry]]
    param display : for debuging
    param einit : initial values can be provided
    param filtering : False = just linear fit

    output
    popt2 : slope, origin, direction 
    pcov2 : Mask for outliers
    '''
    if display:
        import matplotlib.pyplot as plt
        figs, axs = plt.subplots(5,figsize=(8,15))

    nanMask = vxy[0,:].mask
    try:
        binit = np.median(np.sqrt(vxy[1,~nanMask]**2+vxy[0,~nanMask]**2))
        if einit is None:
            einit = np.median(np.arctan2(vxy[1,~nanMask],vxy[0,~nanMask]))
            if einit < 0: einit = einit + 2.0*np.pi
            
        if filtering:
            
            if err is None:
                popt, pcov = curve_fit(linear_speed_func, days[~nanMask], vxy[:,~nanMask].ravel(),p0=(0.,binit,einit),bounds=((-np.inf,0,0) , (np.inf, np.inf,2.0*np.pi) ) )
            else:
                popt, pcov = curve_fit(linear_speed_func, days[~nanMask], vxy[:,~nanMask].ravel(),sigma=err[:,~nanMask].ravel(),p0=(0.,binit,einit),bounds=((-np.inf,0,0) , (np.inf, np.inf,2.0*np.pi) ) )
                
            f = linear_speed_func(days, *popt).reshape(2, len(days))

            range_vx = (np.min(vxy[0,~nanMask]-f[0,~nanMask]),np.max(vxy[0,~nanMask]-f[0,~nanMask])) #(-popt[1]/2.,popt[1]/2.)
            range_vy = (np.min(vxy[1,~nanMask]-f[1,~nanMask]),np.max(vxy[1,~nanMask]-f[1,~nanMask]))

            nx, bins = np.histogram(vxy[0,~nanMask]-f[0,~nanMask], 500, range=range_vx)
            x_vx = bins[:-1]+(bins[1]-bins[0])/2. # bins in x values
            popt_gaussian_vx, _ = curve_fit(gaussian_func, x_vx, nx, p0=(100,0.,50), bounds=( (0,-np.inf,0), np.inf) ) #bounds=([-np.inf, -np.inf, 0], np.inf)
            stddevx = popt_gaussian_vx[2]
            meanx = popt_gaussian_vx[1]

            ny, bins = np.histogram(vxy[1,~nanMask]-f[1,~nanMask], 500, range=range_vy)
            x_vy = bins[:-1]+(bins[1]-bins[0])/2.
            popt_gaussian_vy, _ = curve_fit(gaussian_func, x_vy, ny, p0=(100,0.,50),bounds=( (0,-np.inf,0), np.inf))
            stddevy = popt_gaussian_vy[2]
            meany = popt_gaussian_vy[1]

            if display:
                #print(range_vx)
                axs[0].plot(days[~nanMask]/365.25,np.sqrt(vxy[1,~nanMask]**2+vxy[0,~nanMask]**2),'.')
                axs[0].plot(days[~nanMask]/365.25,np.sqrt(f[0,~nanMask]**2+f[1,~nanMask]**2))
                axs[0].set_xlabel('days since January, 1st 2017')
                axs[0].set_ylabel('v (m/yr)')

                n, bins, patches = axs[1].hist((vxy[0,~nanMask]-f[0,~nanMask]), 500, range=range_vx)
                axs[1].plot(x_vx, gaussian_func(x_vx, *popt_gaussian_vx))
                axs[1].set_title('Distribution linear fit minus data for vx')
                axs[1].set_xlim( range_vx)

                n, bins, patches = axs[2].hist((vxy[1,~nanMask]-f[1,~nanMask]), 500, range=range_vy)
                axs[2].plot(x_vy, gaussian_func(x_vy, *popt_gaussian_vy))
                axs[2].set_title('Distribution linear fit minus data for vy')
                axs[2].set_xlabel('v (m/day)')
                axs[2].set_xlim( range_vy)

                axs[3].plot(days[~nanMask], vxy[0,~nanMask],'.')
                axs[3].set_ylabel('vx (m/year)')
                #axs[3].set_xlabel('days since '+str(time_ref))

                axs[4].plot(days[~nanMask], vxy[1,~nanMask],'.')
                axs[4].set_ylabel('vy (m/year)')
                #axs[4].set_xlabel('days since '+str(time_ref))
            #print(popt)

            nanMask = nanMask | (np.abs(vxy[0,:]-f[0,:]-meanx) > nsigma*stddevx) | (np.abs(vxy[1,:]-f[1,:]-meany) > nsigma*stddevy)
        
        else:
            popt = (0,binit,einit)

        if err is None:
            popt2, pcov2 = curve_fit(linear_speed_func, days[~nanMask], vxy[:,~nanMask].ravel(),p0=popt,bounds=((-np.inf,0,0) , (np.inf, np.inf,2.0*np.pi) ) )
        else:
            popt2, pcov2 = curve_fit(linear_speed_func, days[~nanMask], vxy[:,~nanMask].ravel(),sigma=err[:,~nanMask].ravel(),p0=popt,bounds=((-np.inf,0,0) , (np.inf, np.inf,2.0*np.pi) ) )

        if display:
            axs[0].plot(days[~nanMask]/365.25,np.sqrt(vxy[1,~nanMask]**2+vxy[0,~nanMask]**2),'.')
            axs[0].plot(days[~nanMask]/365.25,np.sqrt(f[0,~nanMask]**2+f[1,~nanMask]**2))
            
            f = linear_speed_func(days, *popt).reshape(2, len(days))
            
            n, bins, patches = axs[1].hist((vxy[0,~nanMask]-f[0,~nanMask]), 500, range=range_vx)
            axs[1].plot(x_vx, gaussian_func(x_vx, *popt_gaussian_vx))
           
            n, bins, patches = axs[2].hist((vxy[1,~nanMask]-f[1,~nanMask]), 500, range=range_vy)
            axs[2].plot(x_vy, gaussian_func(x_vy, *popt_gaussian_vy))

            axs[3].plot(days[~nanMask], vxy[0,~nanMask],'.',color='black')
            axs[3].plot(days[~nanMask], f[0,~nanMask])

            axs[4].plot(days[~nanMask], vxy[1,~nanMask],'.',color='black')
            axs[4].plot(days[~nanMask], f[1,~nanMask])
            
    except:
        popt2 = [np.nan,np.nan,np.nan]
        pcov2 = np.zeros( (3,3) )
    
    return popt2, pcov2, nanMask

def cube2monthly_map_using_linear_regression(cube,year_start,year_end,display=False,fortnightly=False,deltamax=45,deltamin=5,verbose=False):

    '''
        cube from cube_class
        year_start 
        year_end
        display
        fortnightly : 12 points per year if False and 24 points per year if True. point are taken on the 15th of each month if False. 1th and 15th if False.  
        
    '''

    if fortnightly:
        dayofmonth=[1,15]
        n=24
    else:
        dayofmonth=[15]
        n=12

    #creating output arrays
    #v_yearly=np.zeros( (len(range(year_start,year_end+1)),250,250) )
    #direction_yearly=np.zeros( (len(range(year_start,year_end+1)),250,250) )
    v_monthly=np.zeros( (len(range(year_start,year_end+1))*n,250,250) )
    direction_monthly=np.zeros( (len(range(year_start,year_end+1))*n,250,250) )
    #count_yearly = np.zeros( (len(range(year_start,year_end+1)),250,250) )

    #err_v_yearly=np.zeros( (len(range(year_start,year_end+1)),250,250) )
    #err_direction_yearly=np.zeros( (len(range(year_start,year_end+1)),250,250) )
    err_v_monthly=np.zeros( (len(range(year_start,year_end+1))*n,250,250) )
    err_direction_monthly=np.zeros( (len(range(year_start,year_end+1))*n,250,250) )
    count_monthly = np.zeros( (len(range(year_start,year_end+1))*n,250,250) )

    dates_monthly=[]

    for yr in range(year_start,year_end+1):
        
        # select data +/- 1.5 years = fit over 3 years.
        time_ref = datetime.date(yr,7,1)
        d1 = time_ref - datetime.timedelta(days=547)
        d2 = time_ref + datetime.timedelta(days=547) 
        c_year = cube.pick_time_(d1,d2)

        date1 = c_year.date1_(datetime=True)
        date2 = c_year.date2_(datetime=True)
        
        intervals = np.asarray(date2) - np.asarray(date1)
        dates= np.asarray(date1) + (np.asarray(date2)-np.asarray(date1))/2.
        date1=None
        date2=None

        delta = [ d - time_ref for d in dates ] 
        days = np.asarray([ d.days for d in delta])
        delta=None

        #errxy = np.concatenate( [ c_year.errx_(),c_year.erry_() ] ).reshape(2,c_year.nz)
        
        if c_year.nz > 20:
            
            vx = c_year.vx_()
            vy = c_year.vy_()

            errx = c_year.errx_()
            erry = c_year.erry_()
            
            nz = c_year.nz
            c_year = None

            errxy = np.ma.concatenate( [errx, erry ] ).reshape(2,nz)
            errx=None
            erry=None
            vx.mask = vx.mask | vy.mask | (vx > 6000) | (vx < -6000) | (vy > 6000) | (vy < -6000) # IF USED FOR RAPID GLACIERS, THIS SHOULD BE ADJUSTED (GOOD FOR ANTARCTICA)
            vy.mask = vx.mask

            for i in range(250):
                for j in range(250):
            
                    if verbose: progress_bar(j/250*100.,msg='{0}/250 row + year {1}  '.format(i,yr))

                    vxy = np.ma.concatenate( [vx[:,j,i], vy[:,j,i] ] ).reshape(2,nz)
                    #count_yearly[yr-year_start,j,i]=vxy[0,:].count()
    
                    if vxy[0,:].count() > 10:
                        # first fit for filtering

                        popt, _ , nanMask = linear_fit_speed(days, vxy,err=errxy,display=display) #
                        
                        ## remove temporal baseline longer than 60 days
                        if popt[1] < 200:
                            deltam=0
                        else:
                            deltam=deltamin
                        
                        nanMask = nanMask | (intervals > datetime.timedelta(deltamax)) | (intervals <= datetime.timedelta(deltam))
                        nanMask = nanMask | (days < -200) | (days > 200) # keep only 400 days (1100 initially)
                        vxy = np.ma.vstack( (vxy[0,~nanMask], vxy[1,~nanMask]) )
                        errxy1 = np.ma.vstack( (errxy[0,~nanMask], errxy[1,~nanMask]) )
                        dates_masked = dates[~nanMask]
                        nanMask=None
                        
                        #print('popt:',popt)  
                        for mo in range(12):
                            for ida, da in enumerate(dayofmonth):    
                                time_ref2 = datetime.date(yr,mo+1,da)
                                
                                # +/- 1.5 months or 3 months total
                                d1 = time_ref2 - datetime.timedelta(days=deltamax)
                                d2 = time_ref2 + datetime.timedelta(days=deltamax) 
                                
                                dates_masked2 = dates_masked[(dates_masked > d1) & (dates_masked < d2)]
                                delta2 = [ d - time_ref2 for d in dates_masked2 ] 
                                days2 = np.asarray([ d.days for d in delta2])
                                
                                # check if more than 2 values and if they are values before or after central date
                                #print('len(days2):',len(days2)) 
                                if len(days2) > 4 and np.min(days2) <=0 and np.max(days2) > 0:
                                    
                                    vxy2 = vxy[:,(dates_masked > d1) & (dates_masked < d2)]
                                    errxy2 = np.maximum(errxy1[:,(dates_masked > d1) & (dates_masked < d2)],popt[1]*10./100.)
                                    #errxy2 = np.clip(np.zeros(vxy2.shape)+popt[1]*5./100.,10,100)

                                    #errxy2[0,:] = errxy2[0,:]* (np.clip(np.abs(days2)-5,5,40)/5.) 
                                    #errxy2[1,:] = errxy2[1,:]* (np.clip(np.abs(days2)-5,5,40)/5.)

                                    count_monthly[(yr-year_start)*n+mo*len(dayofmonth)+ida,j,i]=len(days2)
                                
                                    try:
                                        # err has been removed all measurement are treated equally !
                                        popt2, pcov2 = curve_fit(linear_speed_func, days2, vxy2.ravel(),sigma=errxy2.ravel(),absolute_sigma=True,p0=popt,bounds=((-np.inf,0,0) , (np.inf, np.inf,2.0*np.pi) ), maxfev=50 )
                                    except:
                                        popt2 = [np.nan,np.nan,np.nan]
                                        pcov2 = np.zeros( (3,3) )
                                    
                                    #print('popt2:',popt2)    
                                    v_monthly[(yr-year_start)*n+mo*len(dayofmonth)+ida,j,i] = popt2[1]
                                    direction_monthly[(yr-year_start)*n+mo*len(dayofmonth)+ida,j,i] = popt2[2]
                                    err_v_monthly[(yr-year_start)*n+mo*len(dayofmonth)+ida,j,i] = pcov2[1,1]
                                    err_direction_monthly[(yr-year_start)*12+mo*len(dayofmonth)+ida,j,i] = pcov2[2,2]

    #time, y and x axis
    xmap = np.arange(cube.geo.xmin,cube.geo.xmin+cube.geo.npix*cube.geo.xposting,cube.geo.xposting)
    ymap = np.arange(cube.geo.ymax,cube.geo.ymax+cube.geo.nrec*cube.geo.yposting,cube.geo.yposting)
    dates_monthly = [datetime.datetime(yr,mo+1,da) for yr in range(year_start,year_end+1) for mo in range(12) for da in dayofmonth]

    if verbose:
        print('len(xmap),len(ymap),len(dates_monthly)',len(xmap),len(ymap),len(dates_monthly))

        print('v_monthly',v_monthly.shape)
        print('direction_monthly',direction_monthly.shape)
        print('err_v_monthly',err_v_monthly.shape)
        print('err_direction_monthly',err_direction_monthly.shape)
        print('count_monthly',count_monthly.shape)

    ds_monthly = None
    ds_monthly = xr.Dataset(
        {"v": (("time","y", "x"), v_monthly),
         "direction": (("time","y", "x"), direction_monthly),
         "err_v": (("time","y", "x"), np.sqrt(err_v_monthly)),
         "err_direction": (("time","y", "x"), np.sqrt(err_direction_monthly)),
         "cnt" : (("time","y", "x"), count_monthly)},
        coords={
            "x": xmap,
            "y": ymap,
            "time": dates_monthly,
        },
    )
    ds_monthly.v.attrs['units'] = "meter/year"
    ds_monthly.v.attrs['long_name']='Surface ice speed (Magnitude)'
    ds_monthly.v.attrs['standard_name']="land_ice_velocity"
    ds_monthly.v.attrs['description'] = 'Surface ice speed obtained after linear regression on Sentinel-1, Landsat-8, Sentinel-2 and Radarsat-2 observations'

    ds_monthly.direction.attrs['units'] = "radian"
    ds_monthly.direction.attrs['long_name']='Direction of the ice flow'
    ds_monthly.direction.attrs['standard_name']="land_ice_direction_velocity"
    ds_monthly.direction.attrs['description'] = 'Flow direction represented as the angle between the flow vector and x-axis. vx = v*cos(direction) and vy = v*sin(direction)'

    ds_monthly.cnt.attrs['units'] = "None"
    ds_monthly.cnt.attrs['long_name']='Number of velocity measurements used'
    ds_monthly.cnt.attrs['standard_name']="land_ice_cnt_velocity"
    ds_monthly.cnt.attrs['description'] = 'Number of measurements used for the regression'

    ds_monthly.err_v.attrs['units'] = "meter/year"
    ds_monthly.err_v.attrs['long_name']='Estimated error (variance) on the surface ice speed (Magnitude)'
    ds_monthly.err_v.attrs['standard_name']="land_ice_error_velocity"
    ds_monthly.err_v.attrs['description'] = 'The error on surface ice speed is given as the variance obtained from the results of the regression.'

    ds_monthly.err_direction.attrs['units'] = "radian"
    ds_monthly.err_direction.attrs['long_name']='Estimated error (variance) on the flow direction'
    ds_monthly.err_direction.attrs['standard_name']="land_ice_error_direction_velocity"
    ds_monthly.err_direction.attrs['description'] = 'The error on surface ice speed is given as the variance obtained from the results of the regression.'    

    ds_monthly.attrs['Conventions'] = "CF-1.7"
    ds_monthly.attrs['Title'] = 'Monthly Time series of surface ice velocity '
    ds_monthly.attrs['Author'] = "Jeremie Mouginot"
    ds_monthly.attrs['version'] = "11mar2022"
    ds_monthly.attrs['Projection'] = 'Polar Stereographic South (71S, 0W)'
    ds_monthly.attrs['proj4'] = "+init=epsg:3413"
    ds_monthly.attrs['datasets'] = 'L8+S1+S2+R2'

    return ds_monthly #, ds_yearly            
