
import pandas as pd
import numpy as np
from bilinear import bilinear

def speed2flowline_with_dates(vx,vy,date1,date2,dt,geo,verbose=False,display=False):

    '''
    vx : numpy array containing displacement in m/yr along x-axis (put NaN for no data)
    vy : numpy array containing displacement in m.yr along y-axis (put NaN for no data)
    date1 : starting date (datetime) ; date1 is older than date2
    date2 : ending date (datetime)
    dt : temporal spacing in days
    geo = geo_param() from fparam.py giving the extent/spacing of vx and vy
    '''

    if vx.shape != vy.shape:
        print('ERROR: vx and vy have different shape.')
        return None, None

    dates = pd.date_range(start=date1,end=date2,freq=str(int(dt))+'D') # not exactly at the position on date2 as date2-date1 is not necessarly a factor of dt

    if verbose:
        #geo.describe()
        print(dates)
        
    n = len(dates)
    positionx = np.zeros( (n+1,geo.nrec,geo.npix) )
    positiony = np.zeros( (n+1,geo.nrec,geo.npix) )

    positionx[0,:,:], positiony[0,:,:] = np.meshgrid(np.arange(geo.npix),np.arange(geo.nrec))
    positionx[0,np.isnan(vx)]=np.nan
    positiony[0,np.isnan(vx)]=np.nan


    for it, t in enumerate(dates[:-1]):
        if verbose:
            print(t)

        # compute the new position based on speed.
        wy, wx = np.where(~np.isnan(positiony[it,:,:]))
        positiony[it+1,wy,wx] = positiony[it,wy,wx]+ bilinear(vy,positionx[it,wy,wx],positiony[it,wy,wx])/365.25*dt/geo.yposting
        positionx[it+1,wy,wx] = positionx[it,wy,wx]+ bilinear(vx,positionx[it,wy,wx],positiony[it,wy,wx])/365.25*dt/geo.xposting

        # check that new position are not outside the map :
        wy, wx = np.where( (positionx[it+1,:,:] > geo.npix) | (positionx[it+1,:,:] < 0) | (positiony[it+1,:,:] > geo.nrec) | (positiony[it+1,:,:] < 0) )

        positionx[it+1,wy,wx] = np.nan
        positiony[it+1,wy,wx] = np.nan

    # compute the last position based on speed.
    dt_final = date2 - dates[-1]
    wy, wx = np.where(~np.isnan(positiony[it,:,:]))
    positiony[-1,wy,wx] = positiony[-2,wy,wx]+ bilinear(vy,positionx[-2,wy,wx],positiony[-2,wy,wx])/365.25*float(dt_final.days)/geo.yposting
    positionx[-1,wy,wx] = positionx[-2,wy,wx]+ bilinear(vx,positionx[-2,wy,wx],positiony[-2,wy,wx])/365.25*float(dt_final.days)/geo.xposting

    # check that last position are not outside the map :
    wy, wx = np.where( (positionx[-1,:,:] > geo.npix) | (positionx[-1,:,:] < 0) | (positiony[-1,:,:] > geo.nrec) | (positiony[-1,:,:] < 0) )

    positionx[-1,wy,wx] = np.nan
    positiony[-1,wy,wx] = np.nan

    return positionx, positiony #