from cube_class import cube_class
import numpy as np
#from numba import jit # numba 0.44 MINIMUM !
import time
from mjd2date import mjd2date, date2mjd # /ST_RELEASE/UTILITIES/PYTHON/mjd2date.py
from datetime import datetime, timedelta
import pandas as pd # pandas 0.25.3 MINIMUM !
import statsmodels.api as sm

''' Supplementary processing functions for cube_class '''

# ====== = ====== DATA FILTRATION ====== = ======
def filter_space(self, space_window=3):
    '''Median space filter with mooving window of given size. Velocity magnitude is filtered and vx-vy are restored.'''
    from scipy import ndimage

    for z in range(self.nz):
        vx = self.data[z].vx
        vy = self.data[z].vy
        # vv = np.sqrt(vx ** 2 + vy ** 2)
        # vvf = np.where(vv.mask == True, np.nan, ndimage.median_filter(vv, size=space_window))
        # vxf = vx * (vvf / vv)
        # vyf = vy * (vvf / vv)
        vxf = np.where(vx.mask == True, np.nan, ndimage.median_filter(vx, size=space_window))
        vyf = np.where(vy.mask == True, np.nan, ndimage.median_filter(vy, size=space_window))
        self.data[z].vx = np.ma.array(data=vxf, mask=[vxf == 0 | self.data[z].vx.mask])
        self.data[z].vy = np.ma.array(data=vyf, mask=[vyf == 0 | self.data[z].vy.mask])
    print('     Filtering is done')
    return self

def filter_time_stat(self, vvar, threshold=1, replace=True, extended_out=False):
    '''
    >> ! << NOT RECOMENDED FOR REGIONS  with strong seasonal dynamic
    VX AND VY ARE FILTERED INDEPENDENTLY
    Simple statistical filter: throw out the vx/vy there abs(time_median - value) > max(pixel_std, stdmax)
    :param vvar: filtering criteria, fraction of the pixel's mean velocity that is acceptable
                    as +- variation from the v[z] value if it is > than STD*threshold (e.g. 1STD, 3STD).
                       Not more than 0.5 could be recomended (~ 0.2).
    :param threshold: STD*-threshold- (e.g. 1STD, 3STD)
    :param replace: True = replace the value in the source cube, False = np.array outputs
    :param extended_out: True = include statistics matrix to outputs
    :return:
            filtered velocity: vxf, vyf
            median and STD maps: medvx, medvy, stdvx, stdvy ( -> extended_out==True)
    '''

    vx = self.vx_().filled(np.nan)
    vy = self.vy_().filled(np.nan)
    stdmaxVX = np.abs(np.nanmedian(vx, 0) * vvar)
    stdmaxVY = np.abs(np.nanmedian(vy, 0) * vvar)
    #  initialize output values {{{
    medvx = np.zeros((self.ny, self.nx), dtype=np.float32)
    medvy = np.zeros((self.ny, self.nx), dtype=np.float32)
    stdvx = np.zeros((self.ny, self.nx), dtype=np.float32)
    stdvy = np.zeros((self.ny, self.nx), dtype=np.float32)
    vxf = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    vyf = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    # RETURN: SOURCE DATA MEDIANS AND STD - 2D.
    for j in range(self.ny):
        for i in range(self.nx):
            medvx[j, i], medvy[j, i], stdvx[j, i], stdvy[j, i] = nb_stat(vx[:, j, i], vy[:, j, i])
    stdvx = np.clip(stdvx, 0.1, stdmaxVX)
    stdvy = np.clip(stdvy, 0.1, stdmaxVY)
    # RETURN: FILTERED SPEEDS - 3D
    for j in range(self.ny):
        for i in range(self.nx):
            vxf[:, j, i], vyf[:, j, i] = nb_filtering(vx[:, j, i], vy[:, j, i], medvx[j, i], medvy[j, i],
                                                      stdvx[j, i], stdvy[j, i], threshold)
    if replace:
        for z in range(self.nz):
            self.data[z].vx[:, :] = vxf[z, :, :]
            self.data[z].vy[:, :] = vyf[z, :, :]
        self.filename=self.filename.replace('.nc', '_TimeStat.nc')
    else:
        if extended_out:
            return vxf, vyf, medvx, medvy, stdvx, stdvy
        else:
            return vxf, vyf  # }}}

def filter_time_RollM(self, RollWindowV=15, RollWindowA=60, RollMinN=3, window=1, verbose=False):
    ''' Apply Rolling MEDIAN filter to velocity magnitude & direction through the pixels (time-dimension), restore vx-vy.
    Output error is local rolling STD of the magnitude.
    /!\ Following the standard rolling behaviour, edges layers where rolling window can't do the processing will be filled with NaN
    :param RollWindowV: the rolling window size in days for the velocity magnitude, int
    :param RollWindowA: the rolling window size in days for the direction (30 days min recommended), int
    :param RollMinN: the minimum quantity of valid values in a window to do calculations, int
    :param window: a space window to do calculation (1 = one pixel, 3 = mean in a 3*3 space box, etc), int '''

    vxf = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    vyf = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    vstdf = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)

    # get data from cube to python variables - decrease extreamly the processing time
    vxg = self.vx_()
    vyg = self.vy_()
    errx = self.errx_()
    erry = self.erry_()
    # make the variables related to the datetime type
    date1 = np.array([d for d in self.date1_()])
    date2 = np.array([d for d in self.date2_()])
    date_c = date2 - ((date2 - date1) / 2)

    # count Rolling Median per lines
    for j in range(self.ny):
        if verbose: print('      line # {0}/{1}'.format(j, self.ny))
        start_time = time.time()
        vxf[:, j, :], vyf[:, j, :], vstdf[:, j, :] = RollMed_line(self, j, vxg, vyg, errx, erry, date_c,
                                                                   RollWindowV, RollWindowA, RollMinN, window)
        t = time.time() - start_time
        if j == 0:
            print('     ---> Estimated filtering time at local machine (not luke-OAR mode): ... ')
            print("            {0} min ({1} sec per line) ".format(t * self.ny / 60, t))

    # prepare data
    df = get_point_data(0, 0, vxg, vyg, errx, erry, date_c, window)
    text = '{} days Rolling Mean. Error is the Rolling STD. '.format(RollWindowV)
    source=[text + df.index.astype(str)[z] for z in range(self.nz)]
    date1 = [date2mjd(df.index[z] - timedelta(RollWindowV / 2)) for z in range(self.nz)]
    date2 = [date2mjd(df.index[z] + timedelta(RollWindowV / 2)) for z in range(self.nz)]
    sensor = ['' for n in range(self.nz)]
    # empty the source cube
    self.insert(vx=[], vy=[], errx=[], erry=[], date1=[], date2=[], source=[], sensor=[], err2D=True, replace=True)
    # fill the cube with filtered data
    self.insert(vx=vxf, vy=vyf, errx=vstdf, erry=vstdf,
                date1=date1, date2=date2, source=source, sensor=sensor, err2D=True)
    self.filename = self.filename.replace('.nc', '_RollM{0}.nc'.format(RollWindowV))
    if verbose: print('  New cube name: ', self.filename)
    return self

def filter_time_lowess_jer(cube, dates_out, lowess_window=20, verbose=False, i=None, j=None, degree=1, it=3, ci=0.95,
                           family='symmetric'):  # {{{

    '''
        dates_out = mjd = modified Julian dates
        ci = Confidence Interval
        it = Iteration
        degree = 1 or 2 (linear or quadratic)
        lowess_window = Number of sample to take into account

        return vx_f, vy_f, vx_std_f, vy_std_f, vx_ci_f, vy_ci_f

    '''
    from skmisc.loess import loess

    date1 = np.array(cube.date1_())
    date2 = np.array(cube.date2_())

    dates = (date1 + date2) / 2.0
    del date1
    del date2

    vx_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vy_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vx_std_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vy_std_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vx_ci_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vy_ci_f = np.zeros((len(dates_out), cube.ny, cube.nx))

    wx = (1. / np.array(cube.errx_())) ** 2
    wy = (1. / np.array(cube.erry_())) ** 2

    if j is None and i is None:
        # processing entire cube
        jj = range(cube.ny)
        ii = range(cube.nx)
    else:
        jj = [j]
        ii = [i]
        # processing entire cube
    for j in jj:
        start_time = time.time()

        for i in ii:
            vxj = cube.vx_(j=j, i=i)
            vyj = cube.vy_(j=j, i=i)

            if vxj.count() == 0:  # no valid data at all
                vx_f[:, j, i] = np.nan
                vy_f[:, j, i] = np.nan

            else:

                frac = np.clip(lowess_window / float(vxj.count()), 0.01, 0.9)
                w, = np.where(~vxj.mask)
                # loess = lowess_bell_shape_kern( dates[w], vxj[w], tau=frac*5)

                lx = loess(dates[w], vxj[w], weights=wx[w], span=frac, degree=degree, family=family, iterations=it)
                lx.fit()

                w2, = np.where((dates_out > np.min(dates[w])) & (dates_out < np.max(dates[w])))

                predvx = lx.predict(dates_out[w2], stderror=True)

                vx_f[w2, j, i] = predvx.values  # np.interp(dates_out,dates2[index_sort],lx[index_sort])
                vx_std_f[w2, j, i] = predvx.stderr
                vx_ci_f[w2, j, i] = (predvx.confidence(1-ci).upper - predvx.confidence(1-ci).lower)/2. # assume symmetric confidence interval

                del lx, predvx

                ly = loess(dates[w], vyj[w], weights=wy[w], span=frac, degree=degree, family=family, iterations=it)
                ly.fit()

                predvy = ly.predict(dates_out[w2], stderror=True)
                vy_f[w2, j, i] = predvy.values
                vy_std_f[w2, j, i] = predvy.stderr
                vy_ci_f[w2, j, i] = (predvy.confidence(1-ci).upper - predvy.confidence(1-ci).lower)/2.
                del ly, predvy
                del w
        t = time.time() - start_time
        if verbose: print('#line, time:', j, t)

    return vx_f, vy_f, vx_std_f, vy_std_f, vx_ci_f, vy_ci_f # }}}

def filter_time_RollBox(cube, time_interval, time_step, boxsize, func='mean', weighted=False, verbose=False):
    '''
    time_interval = ['01-01-2015', '31-12-2019']
    time_step = 'W' # make culculation at earch 7th day, interval for box time shift and output dates
    boxsize = 21 # days; total size of the data collecting time box
    func = 'mean' or 'median'
     '''

    if time_step in ['SM', 'M']: time_step = time_step + 'S'  # steps related to the start of an interval
    dates_out = pd.date_range(start=datetime.strptime(time_interval[0], "%d-%m-%Y"),
                              end=datetime.strptime(time_interval[1], "%d-%m-%Y"), freq=time_step)
    dates_out = dates_out.insert(len(dates_out), pd.datetime.strptime(time_interval[1], "%d-%m-%Y"))

    c_m = cube.copy_structure()
    c_m.filename = c_m.filename.replace('.nc', '_TS{0}_box{1}.nc'.format(time_step, boxsize))

    for i, date_c in enumerate(dates_out):
        start = date_c - timedelta(boxsize / 2)
        end = date_c + timedelta(boxsize / 2)
        if verbose: print(start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"))

        source = 'Mean of period {0} to {1}, time_step={2} box_size={3}'.format(
            start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"), time_step, boxsize)

        # extract cube elements using pick_time2
        cube2 = cube.pick_time2(start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"))
        if cube2.nz != 0:
            if func=='mean':
                if not weighted or cube2.nz==1:
                            c_m.insert(np.nanmean(cube2.vx_(), axis=0), np.nanmean(cube2.vy_(), axis=0),
                            np.nanmean(cube2.errx_(), axis=0), np.nanmean(cube2.erry_(), axis=0),
                            date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=False)
                else:
                    status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = cube2.cube2vel_py()
                    c_m.insert(vxo, vyo, np.nanmean(cube2.errx_(), axis=0), np.nanmean(cube2.erry_(), axis=0),
                                date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=True)

            else:
                c_m.insert(np.nanmedian(cube2.vx_(), axis=0), np.nanmedian(cube2.vy_(), axis=0),
                           np.nanmedian(cube2.errx_(), axis=0), np.nanmedian(cube2.erry_(), axis=0),
                           date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=False)
        else:
            continue
    return c_m

def filter_time_lowess2(cube, dates_out, lowess_window=20, verbose=False, i=None, j=None, lowess='statsmodels',
                        it=3):  # {{{
    '''
    -> sm.nonparametric.lowess, Vx-Vy case
    -> cycle-in-cycle realisation

    Apply LOWESS with re-interpolation to the given dates.
    :return  vx_f, vy_f : np.arrays[len(dates_out), cube.ny, cube.nx]

    :param date_out : date_out as Julian date
    :param i,j : do not processe entire cube, take a given pixel
    '''

    import statsmodels.api as sm

    date1 = np.array(cube.date1_())
    date2 = np.array(cube.date2_())

    dates = (date1 + date2) / 2.0
    del date1
    del date2

    vx_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vy_f = np.zeros((len(dates_out), cube.ny, cube.nx))

    if j is None and i is None:
        # processing entire cube
        jj = range(cube.ny)
        ii = range(cube.nx)
    else:
        jj = [j]
        ii = [i]
        # processing entire cube
    for j in jj:
        start_time = time.time()

        for i in ii:
            vxj = cube.vx_(j=j, i=i)
            vyj = cube.vy_(j=j, i=i)

            if vxj.count() == 0:  # no valid data at all
                vx_f[:, j, i] = np.nan
                vy_f[:, j, i] = np.nan

            else:

                frac = np.clip(lowess_window / float(vxj.count()), 0.01, 0.9)

                w, = np.where(~vxj.mask)

                print(frac, it, float(vxj.count()) * 0.01 )
                loess = sm.nonparametric.lowess(vxj[w], dates[w], frac=frac, it=it,
                                                    delta=0.1, return_sorted=False)

                print(loess)
                index_sort = np.argsort(dates[w])
                dates2 = dates[w] #central dates that are not empty

                vx_f[:, j, i] = np.interp(dates_out, dates2[index_sort], loess[index_sort])

                del loess
                
                loess = sm.nonparametric.lowess(vyj[w], dates[w], frac=frac, it=it,
                                                    delta=0.1, return_sorted=False) #delta=float(vxj.count()) * 0.01

                vy_f[:, j, i] = np.interp(dates_out, dates2[index_sort], loess[index_sort])
                del loess
                del w

        t = time.time() - start_time
        if verbose: print('#line, time:', j, t)

    return vx_f, vy_f  # }}}

def filter_time_lowess3(cube, lowess_window=20, verbose=False, i=None, j=None, it=3, save2cube=''):  # {{{
    '''
    -> sm.nonparametric.lowess, Vx-Vy case
    -> cycle-in-cycle realisation
    '''

    date1 = np.array(cube.date1_())
    date2 = np.array(cube.date2_())

    dates = (date1 + date2) / 2.0
    del date1
    del date2

    vx_f = np.zeros((len(dates),cube.ny,cube.nx))
    vy_f = np.zeros((len(dates),cube.ny,cube.nx))

    if j is None and i is None:
        #processing entire cube
        jj = range(cube.ny)
        ii = range(cube.nx)
    else:
        jj = [j]
        ii = [i]
        # processing entire cube
    for j in jj:
        start_time = time.time()
        
        for i in ii:
            vxj = cube.vx_(j=j,i=i)
            vyj = cube.vy_(j=j,i=i)

            if vxj.count() == 0:  # no valid data at all
                vx_f[:,j,i] = np.nan
                vy_f[:,j,i] = np.nan
            else:
                frac = lowess_window / vxj.count()
                if frac > 0.9: frac = 0.9  # less data in the pixel than -lowess_window- -> frac > 1
                vxj = np.where(vxj.mask == True, np.nan, vxj.data)  # as lowess ignore mask
                vyj = np.where(vyj.mask == True, np.nan, vyj.data)
                # vx_f[:, j, i] = sm.nonparametric.lowess(vxj[w], dates[w], frac=frac, it=it,
                #                               return_sorted=False)
                vx_f[:, j, i] = sm.nonparametric.lowess(vxj, dates, frac=frac, it=it, return_sorted=False)
                while np.sum(~np.isnan(vx_f[:, j, i])) == 0:  # lowess instability bug fixing
                    frac = frac + 0.005
                    vx_f[:, j, i]  = sm.nonparametric.lowess(vxj, dates, frac=frac, it=it,return_sorted=False)

                vy_f[:, j, i] = sm.nonparametric.lowess(vyj, dates, frac=frac, it=it, return_sorted=False)
        t = time.time() - start_time
        if verbose: print('#line, time:',j,t)

    if save2cube!='' :
        print(' RETURN: CUBE VARIABLE')
        c2 = cube.copy_structure()
        text = 'LOWESS VX-VY fitting, day '
        source = [text + str(mjd2date(dates[i])) + '. ' + save2cube for i in range(len(dates))]
        c2.insert(vx_f, vy_f, cube.errx_(), cube.erry_(), dates, dates,
                  source, sensor=['' for i in range(len(dates))], err2D=False)
        c2.filename = c2.filename.replace('.nc', '_LWS.nc')
        print(      'New cube filename: ', c2.filename)
        return c2
    else:
        if verbose: print(' RETURN: 2 NUMPY ARRAYS')
        return vx_f, vy_f # }}}  # }}} # }}}

def filter_time_cylowess(cube, dates_out, lowess_window=20, verbose=False,i=None, j=None, it=3, delta=15.0):
    
    '''

    cylowess.pyx is Cython code for a faster version of the lowess function
    in statsmodels.nonparametric.lowess.

    https://pypi.org/project/cylowess/
    https://github.com/livingsocial/cylowess

    dates_out : 
    '''
    
    import cylowess as cyl
    
    date1 = np.array(cube.date1_()) # modified julian dates
    date2 = np.array(cube.date2_()) # modified julian dates

    dates0 = (date1 + date2) / 2.0 # modified julian dates
    dates = dates0 # modified julian dates
    #dates = np.concatenate( (dates0,dates0,dates0) )
    del date1
    del date2

    date_final = date2mjd(dates_out) # datetinme -> modified julian dates
    
    vx_f = np.zeros((len(dates_out), cube.ny, cube.nx))
    vy_f = np.zeros((len(dates_out), cube.ny, cube.nx))

    #vx0 = cube.vx_()
    #vy0 = cube.vy_()

    if j is None and i is None:
        # processing entire cube
        jj = range(1,cube.ny-1)
        ii = range(1,cube.nx-1)
    else:
        jj = [j]
        ii = [i]
        # processing entire cube
        
    for j in jj:
        start_time = time.time()

        for i in ii:
            vxj = cube.vx_(j=j, i=i) # np.ma.concatenate( (vx0[:,j,i],(vx0[:,j,i+1]+vx0[:,j,i-1])/2., (vx0[:,j+1,i]+vx0[:,j-1,i])/2.) ) #c
            vyj = cube.vy_(j=j, i=i)
            #vxj = vx0[:,j,i]
            #vyj = np.ma.concatenate( (vy0[:,j,i],(vy0[:,j,i+1]+vx0[:,j,i-1])/2., (vy0[:,j+1,i]+vy0[:,j-1,i])/2.) )
            #vyj = vy0[:,j,i]#cube.vy_(j=j, i=i)
            #print(vxj.shape, dates.shape)
            if vxj.count() == 0:  # no valid data at all
                vx_f[:, j, i] = np.nan
                vy_f[:, j, i] = np.nan

            else:
                
                frac = np.clip(lowess_window / float(vxj.count()), 0.001, 0.9)

                print(frac, lowess_window / float(vxj.count()))
                w, = np.where(~vxj.mask)
                
                d = dates[w]
                #delta = (np.max(d) - np.min(d))*0.1
                xloess, yloess = cyl.lowess(vxj[w], d, frac=frac, it=it,delta=delta)
                
                #index_sort = np.argsort(xloess)
                vx_f[:, j, i] = np.interp(date_final, xloess, yloess)
                del yloess
                
                xloess,yloess = cyl.lowess(vyj[w], d, frac=frac, it=it, delta=delta)
                vy_f[:, j, i] = np.interp(date_final, xloess, yloess)
                del yloess
                del w

        t = time.time() - start_time
        if verbose: print('#line, time:', j, t)

    return vx_f, vy_f  # }}}

def filter_time_lowess(cube, lowess_window=20, save2cube='', drop_duplicates = True, mode='replace_orig',  verbose=False):

    '''
    -> sm.nonparametric.lowess, Vx-Vy case
    -> function-in-cycle realisation (copy of the filter_time_lowess but with Vx-Vy)

    Apply LOWESS fitting
    Does NOT re-estimate errors, keeps the source-data values
    :param lowess_window: proxi of number of surrounding points that are taken into accounting for the regression line calculation
    :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
    :param drop_duplicates: before do LOWESS, keep only 1 value per day if many layers are presented (median is taken)
    '''

    if mode in ['replace_orig']:
        pass
    else:
        raise ValueError('-mode- parameter not correctly recognised. See the function description to make a choise.')

    print('     Okey, Google: who can optimize this function? (parallelisation could be implemented)')
    print('     Not you? Well, go get a coffee and waite. It will take some time.')

    vx_f = np.zeros((cube.nz, cube.ny, cube.nx))
    vy_f = np.zeros((cube.nz, cube.ny, cube.nx))

    dates = np.array(cube.date1_()) + np.array(cube.offset_()) / 2
    if drop_duplicates: cube = reset_duplicates(cube, verbose)
    vx_o = cube.vx_()
    vy_o = cube.vy_()

    for j in range(cube.ny):
        start_time = time.time()
        vxj = vx_o[:, j, :]
        vyj = vy_o[:, j, :]
        vxjf, vyjf = lowess_line_vxvy(vxj, vyj, dates, lowess_window=lowess_window)
        for z in range(cube.nz):
            vx_f[z, j, :] = vxjf[z, :]
            vy_f[z, j, :] = vyjf[z, :]

        t = time.time() - start_time
        if j == 0:
            print('     ---> Estimated filtering time at local machine (not luke-OAR mode): ... ')
            print("            {0} min ({1} sec per line) ".format(t * cube.ny / 60, t))
        if verbose: print('      line # {0}/{1}'.format(j, cube.ny))

    if save2cube != '':
        print(' RETURN: CUBE VARIABLE')
        c2 = cube.copy_structure()
        text = 'LOWESS fitting, day '
        source = [text + str(mjd2date(dates[i])) + '. ' + save2cube for i in range(len(dates))]
        c2.insert(vx_f, vy_f, cube.errx_(), cube.erry_(), dates, dates,
                  source, sensor=['' for i in range(len(dates))], err2D=False)
        c2.filename = c2.filename.replace('.nc', '_LWS.nc')
        print('New cube filename: ', c2.filename)
        return c2
    else:
        print(' RETURN: 2 NUMPY ARRAYS')
        return vx_f, vy_f  # }}}  # }}}

def filter_time_spline(cube, dates_out, smooth=0.05, save2cube='', verbose=False, i=None, j=None):  # {{{
    ''' Apply the CUBIC spline fitting
        Does NOT re-estimate errors, fill it with 0
        cube_out.date1 = cube_out.date2 = dates_out
        :param dates_out : date_out as Julian date (np.array) or pandas interval text ('D', 'W', 'SM', etc)
        :param smooth: between 0 and 1 (values between 0.001 and 0.0001 seem to work when weghts are not set),
                        try 0.05 as reference
        :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
        :param i, j: if you want apply the fitting only at one-pixel, xy-yf np.array as output

        need : pip install csaps
    '''
    #import csaps
    from csaps import csaps
    
    date1 = np.array(cube.date1_())
    date2 = np.array(cube.date2_())

    dates_in = (date1 + date2) / 2.0
    wx = (1./np.array(cube.errx_()))**2
    wy = (1./np.array(cube.erry_()))**2

    #if type(dates_out)!=np.ndarray:
    #    try:
    #        add_name = dates_out
    #        dates_out = pd.date_range(start=mjd2date(dates_in.min()), end=mjd2date(dates_in.max()), freq=dates_out, closed='left')
    #        dates_out = np.array([date2mjd(d) for d in dates_out])
    #        dates_out = np.insert(dates_out, [0], dates_in.min())
    #        print ('    len(dates_out)=', len(dates_out))
    #    except: raise ValueError('dates_out variable is not correctly recognised!')

    # make dates uniq
    for d in dates_in:
        if np.sum(dates_in == d) > 1:
            ind, = np.where(dates_in == d)
            dates_in[ind] = dates_in[ind]+ np.arange(len(ind))*0.001

    del date1, date2

    vx_f = np.zeros((len(dates_out),cube.ny,cube.nx))
    vy_f = np.zeros((len(dates_out),cube.ny,cube.nx))

    if j is None and i is None:
        #processing entire cube
        jj = range(cube.ny)
        ii = range(cube.nx)
        flag = True
    else:
        jj = [j]
        ii = [i]
        flag = False

    for j in jj:
        start_time = time.time()
        
        for i in ii:
            vxj = cube.vx_(j=j,i=i)
            vyj = cube.vy_(j=j,i=i)
         
            if vxj.count() == 0:  # no valid data at all
                print(i,j,'no valid data ..')
                vx_f[:,j,i] = np.nan
                vy_f[:,j,i] = np.nan
            
            else:
                
                w, = np.where(~vxj.mask)
                vxj2 = vxj[w]
                dates2 = dates_in[w]
                wx2 = wx[w]
                index_sort = np.argsort(dates2)
                
                # spline parameters
                #0: The smoothing spline is the least-squares straight line fit to the data, 1: The natural cubic spline interpolant
                #sp = csaps.UnivariateCubicSmoothingSpline(dates2[index_sort], vxj2[index_sort], wx2[index_sort], smooth=smooth)
                #vx_f[:,j, i] = sp(dates_out)
                vx_f[:,j, i] = csaps(dates2[index_sort], vxj2[index_sort], dates_out, weights=wx2[index_sort] ,smooth=smooth)
                del vxj2, wx2#, sp

                vyj2 = vyj[w]
                wy2 = wy[w]
                #sp = csaps.UnivariateCubicSmoothingSpline(dates2[index_sort], vyj2[index_sort], wy2[index_sort], smooth=smooth)
                #vy_f[:,j, i] = sp(dates_out)
                vy_f[:,j, i] = csaps(dates2[index_sort], vyj2[index_sort], dates_out,weights=wy2[index_sort], smooth=smooth)
              
                
             
                del vyj2,w,wy2#,sp

        t = time.time() - start_time
        if verbose: print('   line #{0}, time:{1}'.format(j,t))

    ind, = np.where(dates_out > max(dates_in))
    for k in ind:
        vx_f[k,:,:]=0
        vy_f[k,:,:]=0

    if save2cube!='' and flag:
        if verbose: print(' RETURN: CUBE VARIABLE')
        c2 = cube.copy_structure()
        text = 'Cubic spline fitting with -{0}- step, day '.format(add_name)
        source = [text + str(mjd2date(dates_out[i])) + '. ' + save2cube for i in range(len(dates_out))]
        err = np.zeros((len(dates_out))) # errors are filled with 0 TODO TODO TODO
        c2.insert(vx_f, vy_f, err, err, date1=dates_out, date2=dates_out,
                  source=source, sensor=['' for i in range(len(dates_out))], err2D=False)
        try: c2.filename = c2.filename.replace('.nc', '_Spl_{0}.nc'.format(add_name))
        except: c2.filename = c2.filename.replace('.nc', '_Spl.nc')
        print(      'New cube filename: ', c2.filename)
        return c2
    elif save2cube!='' and not flag:
        if verbose: print(' -save2cube- and -i-j- given at the same time. RETURN: 2 NUMPY ARRAYS')
        return vx_f, vy_f # }}}
    else:
        if verbose: print(' RETURN: 2 NUMPY ARRAYS')
        return vx_f, vy_f # }}}

def filter_time_lowess1(cube, lowess_window=20, save2cube='', mode='replace_orig',  verbose=False):  # {{{
    '''
    -> sm.nonparametric.lowess, V-a case
    -> function-in-cycle realisation

    Apply LOWESS fitting to velocity magnitude & direction through the pixels (time-dimension), restore vx-vy.
    Does NOT re-estimate errors, keeps the source-data values
    :param lowess_window: proxi of number of surrounding points that are taken into accounting for the regression line calculation
    :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
    :param mode:
                 'replace_orig' -> replace pixels value with lowess estimation, only for non-masked source pixels
                 [does not yet implemented]'replace_all' -> replace pixels value with lowess estimation for non-masked source pixels, linearly interpolate by time masked pixels
    '''

    if mode in ['replace_orig']:
        pass
    else:
        raise ValueError('-mode- parameter not correctly recognised. See the function description to make a choise.')

    print('     Okey, Google: who can optimize this function? (parallelisation could be implemented)')
    print('     Not you? Well, go get a coffee and waite. It will take some time.')
    
    dates = np.array(cube.date1_()) + np.array(cube.offset_()) / 2
    vx_f = np.zeros((cube.nz,cube.ny,cube.nx))
    vy_f = np.zeros((cube.nz,cube.ny,cube.nx))

    for j in range(cube.ny):
        start_time = time.time()
        vxj = cube.vx_()[:, j, :]
        vyj = cube.vy_()[:, j, :]
        vxjf, vyjf = lowess_line(vxj, vyj, dates, lowess_window=lowess_window)
        for z in range(cube.nz):
            vx_f[z, j, :] = vxjf[z, :]
            vy_f[z, j, :] = vyjf[z, :]

        t = time.time() - start_time
        if j == 0:
            print('     ---> Estimated filtering time at local machine (not luke-OAR mode): ... ')
            print("            {0} min ({1} sec per line) ".format(t * cube.ny / 60, t))
        if verbose: print('      line # {0}/{1}'.format(j, cube.ny))

    if save2cube!='' :
        print(' RETURN: CUBE VARIABLE')
        c2 = cube.copy_structure()
        text = 'LOWESS fitting, day '
        source = [text + str(mjd2date(dates[i])) + '. ' + save2cube for i in range(len(dates))]
        c2.insert(vx_f, vy_f, cube.errx_(), cube.erry_(), dates, dates,
                  source, sensor=['' for i in range(len(dates))], err2D=False)
        c2.filename = c2.filename.replace('.nc', '_LWS.nc')
        print(      'New cube filename: ', c2.filename)
        return c2
    else:
        print(' RETURN: 2 NUMPY ARRAYS')
        return vx_f, vy_f # }}}  # }}}

# ====== = ====== GLOBAL MEANS FOR DATASET ====== = ======

def cube2vel(self, plot_preview=False):  # {{{
    ''':return:
           status
           filtered weighted mean maps: vxo, vyo
           filtered weighted standard deviation maps: stdxo stdyo
           filtered count map: nno
           median maps: medxo medyo
           standard deviation maps: stdxo0 stdyo0
           err maps : errxo erryo
           count map : nno0 '''
    from scipy import ndimage
    
    try:
        #import cube2vel  # FORTRAN core
        from cube2vel_ifort import cube2vel_feather_sp1
    except:
        raise ModuleNotFoundError(
            "You don't have compiled FORTRAN cube2vel CORE CODE. Manage the problem or use the -cube2vel_py-")
    
    print('self.nz',self.nz)
    # NO MAPS - RETURN 0 {{{
    if self.nz == 0:
        return False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # }}}

    start_time = time.time()
    if self.nz == 1:  # ONLY 1 maps, no filtering/averaging possible {{{
        errx = np.zeros((self.ny, self.nx)) + self.data[0].errx
        erry = np.zeros((self.ny, self.ny)) + self.data[0].erry
        return True, \
               np.transpose(self.data[0].vx.filled(fill_value=0)), \
               np.transpose(self.data[0].vy.filled(fill_value=0)), \
               np.zeros((self.ny, self.nx)), \
               np.zeros((self.ny, self.nx)), \
               np.transpose(np.int16(self.data[0].vx.mask)), \
               np.transpose(self.data[0].vx.filled(fill_value=0)), \
               np.transpose(self.data[0].vy.filled(fill_value=0)), \
               np.zeros((self.ny, self.nx)), \
               np.zeros((self.ny, self.nx)), \
               np.transpose(errx.filled(fill_value=0)), \
               np.transpose(erry.filled(fill_value=0)), \
               np.transpose(np.int16(self.data[0].vx.mask)) 
               
        # True, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 # }}}

    vx = self.vx_()

    start_time_w = time.time()
    weightx = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    weighty = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
    feathering = np.float32(
        (vx.mask == 0))  # TODO - the feathering creation AFTER filtering of vx as new invalid pixels are appearing
    feather_size = 10.
    for z in range(self.nz):
        # weightmap = feather * weight, where weight=1/err**2 and feather is the distance to the edge of good value
        weightx[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / self.data[z].errx ** 2
        weighty[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / self.data[z].erry ** 2
    del feathering  # free memory
    print("--- %s seconds for WETHGT job ---" % (time.time() - start_time_w))

    vy = self.vy_()

    # manage the invalid pixels in Fortran style - as 0
    vx = np.where(vx.mask == True, 0, vx)
    vy = np.where(vy.mask == True, 0, vy)

    # switch axis z,y,x (PYTHON) -> x,y,z (FORTRAN) {{{
    errx = self.errx_()
    erry = self.erry_()
    # }}}
    # initialize output values {{{
    medxo = np.zeros((self.nx, self.ny), dtype=np.float32)
    medyo = np.zeros((self.nx, self.ny), dtype=np.float32)
    stdxo0 = np.zeros((self.nx, self.ny), dtype=np.float32)
    stdyo0 = np.zeros((self.nx, self.ny), dtype=np.float32)
    stdxo = np.zeros((self.nx, self.ny), dtype=np.float32)
    stdyo = np.zeros((self.nx, self.ny), dtype=np.float32)
    errxo = np.zeros((self.nx, self.ny), dtype=np.float32)
    erryo = np.zeros((self.nx, self.ny), dtype=np.float32)
    vxo = np.zeros((self.nx, self.ny), dtype=np.float32)
    vyo = np.zeros((self.nx, self.ny), dtype=np.float32)
    nno0 = np.zeros((self.nx, self.ny), dtype=np.int16)
    nno = np.zeros((self.nx, self.ny), dtype=np.int16)
    # }}}
    start_time_c = time.time()

    medxo, medyo, stdxo0, stdyo0, stdxo, stdyo, errxo, erryo, vxo, vyo, nno0, nno = \
        cube2vel_feather_sp1(vx.T, vy.T, errx, erry, weightx.T, weighty.T)
    
    #print(vx.shape,vy.shape,weightx.shape,weighty.shape)
    #print(vx.dtype,vy.dtype,errx.dtype,erry.dtype,weightx.dtype,weighty.dtype)
    del vx,vy,errx,erry,weightx,weighty

    # FORTRAN=PYTHON
    # velox=medxo,veloy=medyo,stdx=stdx0,stdy=stdx0,stdx2=stdxo,stdy2=stdyo,errfx=errxo,errfy=erryo,w_vx=vxo,w_vy=vyo,nn=nno0,nn2=nno
    print("--- %s seconds for FORTRAN CORE job ---" % (time.time() - start_time_c))

    # manage the invalid pixels in Fortran style - as 0
    medxo = np.where(medxo == 0, np.nan, medxo)
    medyo = np.where(medyo == 0, np.nan, medyo)
    stdxo0 = np.where(stdxo0 == 0, np.nan, stdxo0)
    stdyo0 = np.where(stdyo0 == 0, np.nan, stdyo0)
    stdxo = np.where(stdxo == 0, np.nan, stdxo)
    stdyo = np.where(stdyo == 0, np.nan, stdyo)
    errxo = np.where(errxo == 0, np.nan, errxo)
    erryo = np.where(erryo == 0, np.nan, erryo)
    vxo = np.where(vxo == 0, np.nan, vxo)
    vyo = np.where(vyo == 0, np.nan, vyo)
    nno0 = np.where(nno0 == 0, np.nan, nno0)
    nno = np.where(nno == 0, np.nan, nno)
    print("->>> %s seconds for full job ---" % (time.time() - start_time))

    if plot_preview:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].pcolor(np.sqrt(np.rot90(vxo) ** 2 + np.rot90(vyo) ** 2), cmap='coolwarm', vmin=0, vmax=500);
        ax[0, 0].set(title='Filtered mean V')
        ax[0, 1].pcolor(np.sqrt(np.rot90(stdxo) ** 2 + np.rot90(stdyo) ** 2), cmap='coolwarm', vmin=0, vmax=50);
        ax[0, 1].set(title='STD')
        ax[1, 0].pcolor(np.sqrt(np.rot90(errxo) ** 2 + np.rot90(erryo) ** 2), cmap='coolwarm', vmin=0, vmax=50);
        ax[1, 0].set(title='Error')
        ax[1, 1].pcolor(np.rot90(nno), cmap='coolwarm', vmin=0, vmax=100);
        ax[1, 1].set(title='Layers N')
        fig.tight_layout()
        plt.show(block=False)

    return True, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0# }}}

def cube2vel_2(cube, nsigm=0.25, verbose=False):  # {{{
    ''':return:
           status
           filtered weighted mean maps: vxo, vyo
           filtered weighted standard deviation maps: stdxo stdyo
           filtered count map: nno
           median maps: medxo medyo
           standard deviation maps: stdxo0 stdyo0
           err maps : errxo erryo
           count map : nno0 '''
    
    from scipy import ndimage
    
    try:
        import cube2vel_sp2_ifort
    except:
        raise ModuleNotFoundError(
            "You don't have compiled FORTRAN cube2vel CORE CODE. Manage the problem or use the -cube2vel_py-")
    
    print('cube.nz',cube.nz, cube.ny, cube.nx)
    # NO MAPS - RETURN 0 {{{
    if cube.nz == 0:
        return False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # }}}

    if cube.nz == 1:  # ONLY 1 maps, no filtering/averaging possible {{{
        errx = np.zeros((cube.ny, cube.nx)) + cube.data[0].errx
        erry = np.zeros((cube.ny, cube.ny)) + cube.data[0].erry
        return True, \
               np.transpose(cube.data[0].vx.filled(fill_value=0)), \
               np.transpose(cube.data[0].vy.filled(fill_value=0)), \
               np.zeros((cube.ny, cube.nx)), \
               np.zeros((cube.ny, cube.nx)), \
               np.transpose(np.int16(cube.data[0].vx.mask)), \
               np.transpose(cube.data[0].vx.filled(fill_value=0)), \
               np.transpose(cube.data[0].vy.filled(fill_value=0)), \
               np.zeros((cube.ny, cube.nx)), \
               np.zeros((cube.ny, cube.nx)), \
               np.transpose(errx.filled(fill_value=0)), \
               np.transpose(erry.filled(fill_value=0)), \
               np.transpose(np.int16(cube.data[0].vx.mask))  # }}}

    vx = np.float32(cube.vx_())

    weightx = np.zeros((cube.nz, cube.ny, cube.nx), dtype=np.float32)
    weighty = np.zeros((cube.nz, cube.ny, cube.nx), dtype=np.float32)
    feathering = np.float32(
        (vx.mask == 0))  # TODO - the feathering creation AFTER filtering of vx as new invalid pixels are appearing
    feather_size = 10.
    for z in range(cube.nz):
        if verbose: print(z,'feathering ..')
        # weightmap = feather * weight, where weight=1/err**2 and feather is the distance to the edge of good value
        weightx[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / cube.data[z].errx ** 2
        weighty[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / cube.data[z].erry ** 2
    del feathering  # free memory
    if verbose: print('feathering done ..')


    vy = np.float32(cube.vy_())

    # manage the invalid pixels in Fortran style - as 0
    vx = np.where(vx.mask == True, 0, vx)
    vy = np.where(vy.mask == True, 0, vy)
    if verbose: print('set invalid pixels to 0 ..')
        
    # switch axis z,y,x (PYTHON) -> x,y,z (FORTRAN) {{{
    errx = cube.errx_()
    erry = cube.erry_()
    # }}}
    # initialize output values {{{
    medxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    medyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdxo0 = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdyo0 = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    errxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    erryo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    vxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    vyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    nno0 = np.zeros((cube.nx, cube.ny), dtype=np.int32)
    nno = np.zeros((cube.nx, cube.ny), dtype=np.int32)
    if verbose: print('Output matrices created ..')
        
    if verbose:
        print('vx.dtype',vx.dtype)
        print('vy.dtype',vy.dtype)
        print('weightx.dtype',weightx.dtype)
        print('weighty.dtype',weighty.dtype)
        print('errx.dtype',errx.dtype)
        print('erry.dtype',erry.dtype)
        print('nsigm.dtype',np.float32(nsigm).dtype)
        
    # }}}
    start_time_c = time.time()

    #from create_mask import create_mask
    #create_mask('/bettik/jmougino/MASK_GLIMS/World-glacier-outline-simplified.shp',cube.geo, \
    #        output='mask__x'+'{:05d}'.format(cube.i)+'_y'+'{:05d}'.format(cube.j)+'.dat')
    #
    #mask = np.int16(np.fromfile('mask__x'+'{:05d}'.format(cube.i)+'_y'+'{:05d}'.format(cube.j)+'.dat', \
    #             dtype=np.uint8).byteswap().reshape(cube.ny, cube.nx))
    nsigmf = np.float32(nsigm)
    if verbose: print('Launching ..')
    #vxo, vyo = cube2vel_sp2_ifort.cube2vel_feather_sp1(vx.T, vy.T, errx, erry, weightx.T, weighty.T, nsigmf)
     
    medxo, medyo, stdxo0, stdyo0, stdxo, stdyo, errxo, erryo, vxo, vyo, nno0, nno = \
        cube2vel_sp2_ifort.cube2vel_feather_sp1(vx.T, vy.T, errx, erry, weightx.T, weighty.T,np.float32(nsigm))
    
    if verbose: print('Finished ..')
    #print(vx.shape,vy.shape,weightx.shape,weighty.shape)
    #print(vx.dtype,vy.dtype,errx.dtype,erry.dtype,weightx.dtype,weighty.dtype)
    del vx,vy,errx,erry,weightx,weighty
    
    if verbose: print('Not needed variables deleted..')
      
    #print(medxo.dtype,
    # FORTRAN=PYTHON
    # velox=medxo,veloy=medyo,stdx=stdx0,stdy=stdx0,stdx2=stdxo,stdy2=stdyo,errfx=errxo,errfy=erryo,w_vx=vxo,w_vy=vyo,nn=nno0,nn2=nno
    #print("--- %s seconds for FORTRAN cube2vel ---" % (time.time() - start_time_c))

    # manage the invalid pixels in Fortran style - as 0
    medxo = np.where(medxo == 0, np.nan, medxo)
    medyo = np.where(medyo == 0, np.nan, medyo)
    stdxo0 = np.where(stdxo0 == 0, np.nan, stdxo0)
    stdyo0 = np.where(stdyo0 == 0, np.nan, stdyo0)
    stdxo = np.where(stdxo == 0, np.nan, stdxo)
    stdyo = np.where(stdyo == 0, np.nan, stdyo)
    errxo = np.where(errxo == 0, np.nan, errxo)
    erryo = np.where(erryo == 0, np.nan, erryo)
    vxo = np.where(vxo == 0, np.nan, vxo)
    vyo = np.where(vyo == 0, np.nan, vyo)
    nno0 = np.where(nno0 == 0, np.nan, nno0)
    nno = np.where(nno == 0, np.nan, nno)

    if verbose: print('Replace zeros by NaN ..')
    return True, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0# }}}

def cube2vel_3(cube, nsigm=0.25, verbose=False,shp_mask_file='/bettik/jmougino/MASK_GLIMS/World-glacier-outline-simplified.shp'):  # {{{
    ''':return:
           status
           filtered weighted mean maps: vxo, vyo
           filtered weighted standard deviation maps: stdxo stdyo
           filtered count map: nno
           median maps: medxo medyo
           standard deviation maps: stdxo0 stdyo0
           err maps : errxo erryo
           count map : nno0 '''
    
    from scipy import ndimage
    
    try:
        import cube2vel_sp3_ifort
    except:
        raise ModuleNotFoundError(
            "You don't have compiled FORTRAN cube2vel CORE CODE. Manage the problem or use the -cube2vel_py-")
    
    print('cube.nz',cube.nz, cube.ny, cube.nx)
    # NO MAPS - RETURN 0 {{{
    if cube.nz == 0:
        return False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # }}}

    if cube.nz == 1:  # ONLY 1 maps, no filtering/averaging possible {{{
        errx = np.zeros((cube.ny, cube.nx)) + cube.data[0].errx
        erry = np.zeros((cube.ny, cube.ny)) + cube.data[0].erry
        return True, \
               np.transpose(cube.data[0].vx.filled(fill_value=0)), \
               np.transpose(cube.data[0].vy.filled(fill_value=0)), \
               np.zeros((cube.ny, cube.nx)), \
               np.zeros((cube.ny, cube.nx)), \
               np.transpose(np.int16(cube.data[0].vx.mask)), \
               np.transpose(cube.data[0].vx.filled(fill_value=0)), \
               np.transpose(cube.data[0].vy.filled(fill_value=0)), \
               np.zeros((cube.ny, cube.nx)), \
               np.zeros((cube.ny, cube.nx)), \
               np.transpose(errx.filled(fill_value=0)), \
               np.transpose(erry.filled(fill_value=0)), \
               np.transpose(np.int16(cube.data[0].vx.mask))  # }}}

    vx = np.float32(cube.vx_())

    weightx = np.zeros((cube.nz, cube.ny, cube.nx), dtype=np.float32)
    weighty = np.zeros((cube.nz, cube.ny, cube.nx), dtype=np.float32)
    feathering = np.float32(
        (vx.mask == 0))  # TODO - the feathering creation AFTER filtering of vx as new invalid pixels are appearing
    feather_size = 10.
    for z in range(cube.nz):
        if verbose: print(z,'feathering ..')
        # weightmap = feather * weight, where weight=1/err**2 and feather is the distance to the edge of good value
        weightx[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / cube.data[z].errx #** 2
        weighty[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
                                   feather_size) / feather_size / cube.data[z].erry #** 2
    del feathering  # free memory
    if verbose: print('feathering done ..')


    vy = np.float32(cube.vy_())

    # manage the invalid pixels in Fortran style - as 0
    vx = np.where(vx.mask == True, 0, vx)
    vy = np.where(vy.mask == True, 0, vy)
    if verbose: print('set invalid pixels to 0 ..')
        
    # switch axis z,y,x (PYTHON) -> x,y,z (FORTRAN) {{{
    errx = cube.errx_()
    erry = cube.erry_()
    # }}}
    # initialize output values {{{
    medxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    medyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdxo0 = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdyo0 = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    stdyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    errxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    erryo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    vxo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    vyo = np.zeros((cube.nx, cube.ny), dtype=np.float32)
    nno0 = np.zeros((cube.nx, cube.ny), dtype=np.int16)
    nno = np.zeros((cube.nx, cube.ny), dtype=np.int16)
    if verbose: print('Output matrices created ..')
    # }}}
    start_time_c = time.time()

    from create_mask import create_mask
    create_mask(shp_mask_file,cube.geo, \
            output='mask__x'+'{:05d}'.format(cube.i)+'_y'+'{:05d}'.format(cube.j)+'.dat')
    mask = np.int16(np.fromfile('mask__x'+'{:05d}'.format(cube.i)+'_y'+'{:05d}'.format(cube.j)+'.dat', \
                 dtype=np.uint8).byteswap().reshape(cube.ny, cube.nx))
    
    lapse = np.float32(np.asarray(cube.offset_()))
    
    if verbose: print('Launching ..')
    medxo, medyo, stdxo0, stdyo0, stdxo, stdyo, errxo, erryo, vxo, vyo, nno0, nno = \
        cube2vel_sp3_ifort.cube2vel_feather_sp1(vx.T, vy.T, errx, erry, weightx.T, weighty.T,mask.T,lapse, np.float32(nsigm))
    if verbose: print('Finished ..')
    del vx,vy,errx,erry,weightx,weighty
    
    #print(medxo.dtype,
    # FORTRAN=PYTHON
    # velox=medxo,veloy=medyo,stdx=stdx0,stdy=stdx0,stdx2=stdxo,stdy2=stdyo,errfx=errxo,errfy=erryo,w_vx=vxo,w_vy=vyo,nn=nno0,nn2=nno
    #print("--- %s seconds for FORTRAN cube2vel ---" % (time.time() - start_time_c))

    # manage the invalid pixels in Fortran style - as 0
    medxo = np.where(medxo == 0, np.nan, medxo)
    medyo = np.where(medyo == 0, np.nan, medyo)
    stdxo0 = np.where(stdxo0 == 0, np.nan, stdxo0)
    stdyo0 = np.where(stdyo0 == 0, np.nan, stdyo0)
    stdxo = np.where(stdxo == 0, np.nan, stdxo)
    stdyo = np.where(stdyo == 0, np.nan, stdyo)
    errxo = np.where(errxo == 0, np.nan, errxo)
    erryo = np.where(erryo == 0, np.nan, erryo)
    vxo = np.where(vxo == 0, np.nan, vxo)
    vyo = np.where(vyo == 0, np.nan, vyo)
    nno0 = np.where(nno0 == 0, np.nan, nno0)
    nno = np.where(nno == 0, np.nan, nno)

    return True, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0# }}}

# def cube2vel_py(self, time_filter=False, plot_preview=False, save2cube='', verbose=False):
#     '''
#     Calculate weighted mean velocity and some statistics.
#     :param time_filter: apply time filter filter_time_stat with given number as -stat_stdmax- threshold (filter similar to Fortran), int
#                         if False -> don't apply filtering
#     :param plot_preview: plot the mean velocity, weighted std and error, layers N, bool
#     :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
#     :return:
#            status
#            filtered weighted mean maps: vxo, vyo
#            filtered weighted standard deviation maps: stdxo stdyo
#            filtered count map: nno
#            median maps: medxo medyo
#            standard deviation maps: stdxo0 stdyo0
#            err maps : errxo erryo
#            count map : nno0
#     '''
#     from scipy import ndimage

#     if self.data[0].errx.shape != ():  flag=True
#     else: flag = False

#     if self.nz == 0:
#         print('     WARNING: empty cube')
#         return False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

#     if verbose:
#         print("->>>>>> Use the force of Numba! It's a new cube2vel with the PYTHON core. ")
#         print( "->>>>>> It works quicker for big and huge cubes at a personal laptop (test's time is comparable or better towards Fortran core).")
#         print("->>>>>> Calculations are slightly different from the Fortran code.")
#     start_time = time.time()
#     # {{{
#     errx = np.array(self.errx_())
#     erry = np.array(self.erry_())
#     w_meanvxf = np.zeros((self.ny, self.nx), dtype=np.float32)
#     w_meanvyf = np.zeros((self.ny, self.nx), dtype=np.float32)
#     w_errxf = np.zeros((self.ny, self.nx), dtype=np.float32)
#     w_erryf = np.zeros((self.ny, self.nx), dtype=np.float32)
#     w_stdvxf = np.zeros((self.ny, self.nx), dtype=np.float32)
#     w_stdvyf = np.zeros((self.ny, self.nx), dtype=np.float32)  # }}}

#     if time_filter:
#         vxf, vyf, medvx, medvy, stdvx, stdvy = self.filter_time_stat(stat_stdmax=time_filter, replace=False, extended_out=True)
#     else:
#         vxf = self.vx_().filled(np.nan)
#         vyf = self.vy_().filled(np.nan)
#         medvx = np.round(np.nanmedian(vxf, 0), 2)
#         medvy = np.round(np.nanmedian(vyf, 0), 2)
#         stdvx = np.round(np.nanstd(vxf, 0), 2)
#         stdvy = np.round(np.nanstd(vyf, 0), 2)

#     # RETURN: NUMBER OF VALID PIXELS
#     n = np.count_nonzero(~np.isnan(self.vx_()), axis=0)
#     n = np.where(n == 0, np.nan, n)
#     # RETURN: NUMBER OF VALID PIXELS after filtering
#     nf = np.count_nonzero(~np.isnan(vxf), axis=0)
#     nf = np.where(nf == 0, np.nan, nf)

#     # RETURN : WEIGHTS
#     start_time_w = time.time()
#     weightx = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
#     weighty = np.zeros((self.nz, self.ny, self.nx), dtype=np.float32)
#     feathering = np.float32(~np.isnan(vxf))
#     feather_size = np.float32(5.)
#     for z in range(self.nz):
#         if np.nansum(errx[z]) == 0: # case for post-processed data with no-filled error; equal weight values
#             weightx[z, :, :] = 1
#             weighty[z, :, :] = 1
#         else:
#             # weightmap = feather * weight, where weight=1/err**2 and feather is the distance to the edge of good value
#             weightx[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
#                                        feather_size) / feather_size / errx[z] ** 2
#             weighty[z, :, :] = np.clip(ndimage.morphology.distance_transform_edt(feathering[z, :, :]), 0.,
#                                        feather_size) / feather_size / erry[z] ** 2
#     del feathering  # free memory
#     if  verbose: print("--- %s seconds for WETHGT job ---" % (time.time() - start_time_w))

#     # RETURN : WEIGHTED MEAN - 2D. Use the force of numba magic
#     start_time_w = time.time()
#     for j in range(self.ny):
#         for i in range(self.nx):
#             try:
#                 w_meanvxf[j, i], w_meanvyf[j, i] = nb_wmean(vxf[:, j, i], vyf[:, j, i], weightx[:, j, i],
#                                                             weighty[:, j, i])
#             except:
#                 w_meanvxf[j, i], w_meanvyf[j, i] = np.nan, np.nan  # avoid division per 0 error
#     if verbose: print("--- %s seconds for MEAN_VEL job ---" % (time.time() - start_time_w))
#     # RETURN : WEIGHTED ERROR - 2D. Use the force of numba magic
#     start_time_w = time.time()
#     if flag:
#         terrx = errx * weightx
#         terry = erry * weighty
#     else:
#         terrx = np.einsum('z,zji->zji', errx, weightx)
#         terry = np.einsum('z,zji->zji', erry, weighty)
#     for j in range(self.ny):
#         for i in range(self.nx):
#             try:
#                 w_errxf[j, i], w_erryf[j, i] = nb_werr(terrx[:, j, i], terry[:, j, i], weightx[:, j, i],
#                                                        weighty[:, j, i])
#             except:
#                 w_errxf[j, i], w_erryf[j, i] = np.nan, np.nan  # avoid division per 0 error
#     if verbose: print("--- %s seconds for WEIGHTED_ERR job ---" % (time.time() - start_time_w))
#     # RETURN : WEIGHTED STD - 2D. Use the force of numba magic
#     start_time_w = time.time()
#     for j in range(self.ny):
#         for i in range(self.nx):
#             try:
#                 w_stdvxf[j, i], w_stdvyf[j, i] = nb_wstd(vxf[:, j, i], vyf[:, j, i], weightx[:, j, i],
#                                                          weighty[:, j, i], nf[i, j])
#             except:
#                 w_stdvxf[j, i], w_stdvyf[j, i] = np.nan, np.nan  # avoid division per 0 error
#     if verbose: print("--- %s seconds for WEIGHTED_STD job ---" % (time.time() - start_time_w))

#     if verbose: print("->>>> %s seconds for full job ---" % (time.time() - start_time))
#     # keep old-vertion return names:
#     vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = \
#         w_meanvxf, w_meanvyf, w_stdvxf, w_stdvyf, nf, medvx, medvy, stdvx, stdvy, w_errxf, w_erryf, n
#     del w_meanvxf, w_meanvyf, w_stdvxf, w_stdvyf, nf, medvx, medvy, stdvx, stdvy, w_errxf, w_erryf, n

#     if plot_preview:
#         import matplotlib.pyplot as plt
#         fig, ax = plt.subplots(2, 2)
#         ax[0, 0].pcolor(np.sqrt(np.flip(vxo, 0) ** 2 + np.flip(vyo, 0) ** 2), cmap='coolwarm', vmin=0, vmax=500);
#         ax[0, 0].set(title='Filtered mean V')
#         ax[0, 1].pcolor(np.sqrt(np.flip(stdxo, 0) ** 2 + np.flip(stdyo, 0) ** 2), cmap='coolwarm', vmin=0, vmax=50);
#         ax[0, 1].set(title='STD')
#         ax[1, 0].pcolor(np.sqrt(np.flip(errxo, 0) ** 2 + np.flip(erryo, 0) ** 2), cmap='coolwarm', vmin=0, vmax=50);
#         ax[1, 0].set(title='Error')
#         ax[1, 1].pcolor(np.flip(nno, 0), cmap='coolwarm', vmin=0, vmax=100);
#         ax[1, 1].set(title='Layers N')
#         fig.tight_layout()
#         plt.show(block=False)

#     if save2cube != '':
#         c2 = self.copy_structure()
#         date1 = np.array(self.date1_()).min()
#         date2 = np.array(self.date1_()).max()
#         text = 'cube2vel_py output (weighted mean for velocity), period  '
#         source = text + '{0} to {1}'.format(str(mjd2date(date1)), str(mjd2date(date2))) + '. ' + save2cube
#         c2.insert(vxo, vyo, errxo, erryo, date1, date2, source, sensor='', err2D=True)
#         c2.filename = c2.filename.replace('.nc', '_c2v.')
#         print('     New cube.filename: ',c2.filename)
#         return c2
#     else:
#         return True, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0

# ====== = ====== TIME-BASED MEANS FOR DATASET ====== = ======

def time_regularisation(cube, time_interval, time_step, func='mean', weighted = False, verbose=False):
    ''' Regularize z time step (date1-date2). When working, each velocity layer is attributed to the central date of the date1-date2 interval.
        Does NOT fill data gaps or empty time intervals.
        Does NOT re-estimate errors for weighted=False
    :param cube: cube_class variable
    :param time_interval: start-end dates, ['01-01-2015', '31-12-2019']
    :param time_step: 'D': day, 'W': week, 'SM': 2 weeks, 'M': month, 'Y': year
                    use the number+letter to derive your interval, such as '6M'=half of year
    :param func: function of data aggregation, 'mean' or 'median' or 'std'
    :param weighted: False = simple numpy calculations, True = weighted calculations with cube2vel
    :return: cube_class variables '''

    if float(pd.__version__[:-2]) < 0.25:
        raise EnvironmentError('    >>> Your -pandas- library is too old. Update it! v0.25.3 minimum is required. <<<')

    c_m = cube.copy_structure()
    if cube.data[0].errx.shape == (): err2D=False
    else: err2D=True

    if func=='mean':
        c_m.filename = c_m.filename.replace('.nc', '_{0}.nc'.format(time_step))
        if weighted:
            source_pr = 'Weighted mean by interval'
        else: source_pr = 'Mean by interval'
    elif func=='median':
        c_m.filename = c_m.filename.replace('.nc', '_{0}_Med.nc'.format(time_step))
        source_pr = 'Median by interval'
    elif func=='std':
        c_m.filename = c_m.filename.replace('.nc', '_{0}_STD.nc'.format(time_step))
        source_pr = 'STD by interval'
    else: raise ValueError ('-func- variable are not correct')

    if time_step in ['SM', 'M']: time_step = time_step+'S' # steps related to the start of an interval
    dates_out = pd.date_range(start=datetime.strptime(time_interval[0], "%d-%m-%Y"), end=datetime.strptime(time_interval[1], "%d-%m-%Y"), freq=time_step)
    dates_out = dates_out.insert(len(dates_out), pd.datetime.strptime(time_interval[1],"%d-%m-%Y"))

    for i, start in enumerate(dates_out[0:-1]):
        end = dates_out[i+1] - timedelta(1)
        if i+2==len(dates_out): end=end + timedelta(1) # just make it look nicer, nothing else
        if verbose: print(start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"))

        source = source_pr + ', period {0} to {1}'.format(start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"))

        # extract cube elements using pick_time2
        cube2 = cube.pick_time2(start.strftime("%d-%m-%Y"), end.strftime("%d-%m-%Y"))
        # calculate the average with defined -func-
        c_m = stack2layer(cube2, c_m, func, start, end, source, err2D, weighted, verbose)

    print(' >>> time_reg FINISH')
    return c_m

def average_year(cube, time_step, time_interval=-1, func='median', weighted=False, verbose=False):
    ''' Calculate the typical avereg velocity by given -time_steps- of an entire year
        using the source layers from the given -time_interval-.
        Does NOT fill data gaps or empty time intervals.
        Does NOT reestimate errors.
    :param cube: cube_class variable
    :param time_interval: start-end dates in dd-mm format or "-1" for the entire year, ['21-04', '31-12']
    :param time_step: 'D': day, 'W': week, 'SM': 2 weeks, 'M': month, 'Y': year
                    use the number+letter to derive your interval, such as '6M'=half of year
    :param func: function of data aggregation, 'mean' or 'median' or 'std'
    :param weighted: False = simple numpy calculations, True = weighted calculations with cube2vel
    :return: cube_class variables '''

    if float(pd.__version__[:-2]) < 0.25:
        raise EnvironmentError('    >>> Your -pandas- library is too old. Update it! v0.25.3 minimum is required. <<<')

    if time_interval==-1: time_interval=['01-01', '31-12']

    c_m = cube.copy_structure()
    if cube.data[0].errx.shape == (): err2D=False
    else: err2D=True

    # manage metadata labels
    if func=='mean':
        c_m.filename = c_m.filename.replace('.nc', '_AvY_{0}_Mean.nc'.format(time_step))
        source_pr = 'Mean by interval: vx=mean, vy=mean, error=STD / sqrt(N_valide_pix).'
    elif func=='median':
        c_m.filename = c_m.filename.replace('.nc', '_AvY_{0}.nc'.format(time_step))
        source_pr = 'Median by interval: vx=median, vy=median, error=STD / sqrt(N_valide_pix).'
    elif func=='std':
        c_m.filename = c_m.filename.replace('.nc', '_AvY_{0}_STD.nc'.format(time_step))
        source_pr = 'STD by interval: vx=stdx, vy=stdy, error=N_valide_pix.'
    else: raise ValueError ('-func- variable are not correct')

    # manage the dates in nice format
    if time_step in ['SM', 'M']:
        time_step = time_step+'S' # steps related to the start of an interval
        dates_out = pd.date_range(start=datetime.strptime(time_interval[0], "%d-%m"), end=datetime.strptime(time_interval[1], "%d-%m"), freq=time_step)
    elif time_step not in ['W', 'D']:
        dates_out = pd.date_range(start=datetime.strptime(time_interval[0], "%d-%m"),end=datetime.strptime(time_interval[1], "%d-%m"), freq=time_step)
        dd = []
        for d in dates_out: dd.append(d.replace(day=1)) # pass per list as pandas_DateTimeIndex can't do that
        dates_out = pd.to_datetime(dd); del dd, d
    if not pd.datetime.strptime(time_interval[1],"%d-%m") == dates_out[-1]: # add it at the end to close the interval
        dates_out = dates_out.insert(len(dates_out), pd.datetime.strptime(time_interval[1],"%d-%m"))
    print('  OUTPUT nz = ', len(dates_out)-1)

    for i, start in enumerate(dates_out[0:-1]):
        end = dates_out[i+1] - timedelta(1)
        if i+2==len(dates_out): end=end + timedelta(1) # just make it look nicer, nothing else

        source = source_pr + 'Period {0} to {1}'.format(start.strftime("%d-%m"), end.strftime("%d-%m"))
        # extract cube elements for given -start-end- interval
        cube2 = cube.copy_structure()
        for z in range(cube.nz):
            d_c = mjd2date( (cube.data[z].date1 + cube.data[z].date2)/2 )
            if d_c.month == 2 and d_c.day == 29 : d_c = d_c.replace(day=28) # manage bissextile year
            d_c = d_c.replace(year=1900) # delete year info
            if start <= d_c <= end :
                cube2.insert(cube.data[z].vx, cube.data[z].vy, cube.data[z].errx, cube.data[z].erry,
                             cube.data[z].date1, cube.data[z].date2, cube.data[z].source, cube.data[z].sensor,
                             err2D=err2D, replace=False, verbose=False)
        if verbose: print('  period {0} to {1}, loc_cube.nz={2}'.format(start.strftime("%d-%m"), end.strftime("%d-%m"), cube2.nz))
        # calculate the average with defined -func-
        c_m = stack2layer(cube2, c_m, func, start, end, source, err2D, weighted, verbose=False)

    print(' >>> average_year FINISH')
    return c_m

# ========= = ========== = ========== = ======== = ========== = ===========
#  internal sub-functions  {{{

# @jit(nopython=True, fastmath=True)
# def nb_stat(vx, vy):
#     # COMPUTE MEDIAN
#     medvx = np.round(np.nanmedian(vx), 2)
#     medvy = np.round(np.nanmedian(vy), 2)
#     # COMPUTE STD
#     stdvx = np.round(np.nanstd(vx), 2)
#     stdvy = np.round(np.nanstd(vy), 2)
#     return medvx, medvy, stdvx, stdvy

# @jit(nopython=True, parallel=True, fastmath=True)
# def nb_filtering(vx, vy, medvx, medvy, stdvx, stdvy, threshold=1):
#     #REMOVE VALUES THAT ARE MORE THAN threshold*STD AWAY FOR MEDIAN
#     vxf = np.where(np.absolute(medvx - vx) <= threshold*stdvx, vx, np.nan)
#     vyf = np.where(np.absolute(medvy - vy) <= threshold*stdvy, vy, np.nan)
#     return vxf, vyf

# @jit(nopython=True, parallel=True, fastmath=True)
# def nb_wmean(vxf, vyf, weightx, weighty):
#     # WEIGHTED MEAN VELOCITY
#     w_meanvxf = np.round(np.nansum(vxf * weightx) / np.nansum(weightx), 2)
#     w_meanvyf = np.round(np.nansum(vyf * weighty) / np.nansum(weighty), 2)
#     return w_meanvxf, w_meanvyf

# @jit(nopython=True, fastmath=True)
# def nb_werr(terrx, terry, weightx, weighty):
#     # WEIGHTED ERROR
#     w_errxf = np.round(np.sqrt(np.nansum(terrx ** 2) / (np.nansum(weightx) ** 2)), 2)
#     w_erryf = np.round(np.sqrt(np.nansum(terry ** 2) / (np.nansum(weighty) ** 2)), 2)
#     return w_errxf, w_erryf

# @jit(nopython=True, parallel=True, fastmath=True)
# def nb_wstd(vxf, vyf, weightx, weighty, nf):
#     w_stdvxf = np.round(np.sqrt((np.nansum(weightx * ((vxf - np.nanmean(vxf)) ** 2) / np.nansum(weightx))) * (nf / (nf - 1))),2)
#     w_stdvyf = np.round(np.sqrt((np.nansum(weighty * ((vyf - np.nanmean(vyf)) ** 2) / np.nansum(weighty))) * (nf / (nf - 1))),2)
#     return w_stdvxf, w_stdvyf

# @jit(nopython=True, parallel=True, fastmath=True, error_model="numpy")
# def nb_wstd2(values, weights):
#     w_std = np.round(np.sqrt(np.nansum(weights * ((values - np.nanmean(values)) ** 2)) / np.nansum(weights)), 2)
#     return w_std
# }}}

def lowess_line(vxj, vyj, dates, lowess_window=20): #{{{
    '''Apply LOWESS filter to velocity magnitud throuth the pixel time-dataset, restore vx-vy
    :param lowess_window: proxi of number of surrounding points that are taken into accounting for the regression line calculation
    '''
    vxjf = np.zeros((vxj.shape[0], vxj.shape[1]), dtype=np.float32)
    vyjf = np.zeros((vyj.shape[0], vyj.shape[1]), dtype=np.float32)

    vv = np.sqrt(vxj ** 2 + vyj ** 2)
    vv = np.where(vv.mask == True, np.nan, vv.data)  # as lowess ignore mask
    a = np.where(vyj > 0, np.rad2deg(np.arccos(vxj / vv)), 360 - np.rad2deg(np.arccos(vxj / vv)))

    for i in range(vxj.shape[1]):
        if np.sum(~np.isnan(vv[:, i])) == 0:  # no valid data at all
            vvf = af = np.empty((vxj.shape[1], vxj.shape[1])) * np.nan
        else:
            flag=False
            # Velocity magnitude fitting
            frac = lowess_window / np.sum(~np.isnan(vv[:, i]))
            if frac > 0.9: frac=0.9 # less data in the pixel than -lowess_window- -> frac > 1
            vvf = sm.nonparametric.lowess(vv[:, i], dates, frac=frac, it=3,
                                                delta=np.sum(~np.isnan(vv[:, i])) * 0.01, return_sorted=False)
            while np.sum(~np.isnan(vvf))==0: # lowess instability bug fixing
                frac=frac+0.005
                vvf = sm.nonparametric.lowess(vv[:, i], dates, frac=frac, it=3,
                                              delta=np.sum(~np.isnan(vv[:, i])) * 0.01, return_sorted=False)
            # Direction fitting
            # manage 0-360 zone: if 50% of points are in the 90deg sector around the OX axe, tern the ax to do the fit
            if np.quantile(a[:, i], 0.25) + (360 - np.quantile(a[:, i], 0.75)) < 90:
                flag = True
                a[:, i] = np.where(a[:, i] < 270, a[:, i] + 90, a[:, i] - 270)
            frac = 50 / np.sum(~np.isnan(a[:, i]))
            if frac > 0.9: frac = 0.9  # less data in the pixel than -lowess_window- -> frac > 1
            af = sm.nonparametric.lowess(a[:, i], dates, frac=frac, it=2, return_sorted=False)
            if flag:
                af = np.where(af > 90, af - 90, 270 + af)

        vxjf[:, i] = vvf * np.round(np.cos(np.deg2rad(af)), 2)
        vyjf[:, i] = vvf * np.round(np.sin(np.deg2rad(af)), 2)

    return vxjf, vyjf # }}}

def get_point_data(i, j, vxg, vyg, errx, erry, date_c, window):  # {{{
   '''
   Get data from the cube at [i,j] point as pandas dataframe with averaging in -window-
   SET INDEXING AS DATETIME DATES + KEEP DATES COLUMN AS MJD INTEGER
   :return: pandas Dataframe
   '''

   w = round(window / 2)
   date_as_date = np.array([mjd2date(d) for d in date_c])
   # velocity in W*W window median
   vx = np.nanmedian(np.nanmedian(vxg[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
   vy = np.nanmedian(np.nanmedian(vyg[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
   vv = np.round(np.sqrt(vx ** 2 + vy ** 2), 2)
   a = np.where(vy>0, np.rad2deg(np.arccos(vx/vv)), 360-np.rad2deg(np.arccos(vx/vv))) # direction
   # TODO : decide for error
   try: # 2D error
       errx = np.nanmedian(np.nanmedian(errx[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
       erry = np.nanmedian(np.nanmedian(erry[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
       error = np.round(np.sqrt(errx ** 2 + erry ** 2) / 2, 2)
   except: # 1D error
       error = np.round(np.sqrt(errx ** 2 + erry ** 2) / 2, 2)
   #error = abs(np.round((cube.errx_(i,j) * cube.vx_(i, j) + cube.erry_(i,j) * cube.vy_(i, j)) / vv, 2)) # VV error
   #error = cube.errx_(i,j)
   if not vv.any():
        try: # point is closer to the edge than the window size
            vx = vxg[:, j, i]
            vy = vyg[:, j, i]
            vv = np.round(np.sqrt(vx ** 2 + vy ** 2), 2)
            a = np.where(vy > 0, np.rad2deg(np.arccos(vx / vv)), 360 - np.rad2deg(np.arccos(vx / vv)))
        except:
            #print('    NO VALID VELOCITY DATA at x:y {0}:{1}')
            return -1

   # combine data to Pandas.DataFrame
   d = {'date_c': pd.Series(date_as_date, dtype='datetime64[ns]'),
        'vx': pd.Series(vx),
        'vy': pd.Series(vy),
        'vv': pd.Series(vv),
        'a': pd.Series(a), # direction
        'error': pd.Series(error),
        'date': pd.Series(date_c)}
   df = pd.DataFrame(d)

   df = df.set_index(['date_c']).sort_index()

   return df  # }}}

def RollMed_line(cube, j, vxg, vyg, errx, erry, date_c, RollWindowV, RollWindowA, RollMinN, window):

    vxfj = np.zeros((cube.nz, cube.nx), dtype=np.float32)
    vyfj = np.zeros((cube.nz, cube.nx), dtype=np.float32)
    vstdj = np.zeros((cube.nz, cube.nx), dtype=np.float32)

    for i in range(cube.nx):
        flag = False
        df = get_point_data(i, j, vxg, vyg, errx, erry, date_c, window)

        # manage 0-360 zone: if 50% of points are in the 90deg sector around the OX axe, tern the ax to do the fit
        if df.a.quantile(0.25) + (360 - df.a.quantile(0.75)) < 90:
            flag=True
            df.a=np.where(df.a<270, df.a+90, df.a-270)

        # prefilter - 30-days 3STD
        # df = df[np.abs(df.vv - df.vv.rolling(window=str(30) + 'D', min_periods=3).median()) <=
        #                  3 * df.vv.rolling(window=str(30) + 'D', min_periods=3).std()]
        df.mask(np.abs(df.vv - df.vv.rolling(window=str(30) + 'D', min_periods=3).median()) >
                         3 * df.vv.rolling(window=str(30) + 'D', min_periods=3).std(), np.nan, inplace=True)
        # count rolling median
        dff = df.rolling(window='{0}D'.format(RollWindowV), min_periods=RollMinN).median()
        dff['a'] = df['a'].rolling(window='{0}D'.format(RollWindowA), min_periods=RollMinN).median()
        dfs = pd.DataFrame(df.rolling(window='{0}D'.format(RollWindowV), min_periods=RollMinN).std())
        if flag:
            dff.a=np.where(dff.a>90, dff.a-90, 270 + dff.a)
        dff.index = dff.index - timedelta(RollWindowV / 2)
        dfs.index = dfs.index - timedelta(RollWindowV / 2)

        vxfj[:, i] = dff.vv * np.round(np.cos(np.deg2rad(dff.a)), 2)
        vyfj[:, i] = dff.vv * np.round(np.sin(np.deg2rad(dff.a)), 2)
        vstdj[:, i] = dfs.vv

        vxfj[:, i] = dff.vx
        vyfj[:, i] = dff.vy
        vstdj[:, i] = dfs.vv

    return vxfj, vyfj, vstdj

def lowess_line_vxvy(vxj, vyj, dates, lowess_window=20): #{{{
    '''Apply LOWESS filter to velocity magnitud throuth the pixel time-dataset, restore vx-vy
    :param lowess_window: proxi of number of surrounding points that are taken into accounting for the regression line calculation
    '''
    vxjf = np.zeros((vxj.shape[0], vxj.shape[1]), dtype=np.float32)
    vyjf = np.zeros((vyj.shape[0], vyj.shape[1]), dtype=np.float32)

    vxj = vxj.filled(np.nan) ## np.where(vxj.mask == True, np.nan, vxj.data)  # as lowess ignore mask
    vyj = vyj.filled(np.nan) ## np.where(vyj.mask == True, np.nan, vyj.data)

    for i in range(vxj.shape[1]):
        if np.sum(~np.isnan(vyj[:, i])) == 0:  # no valid data at all
            vxjf[:, i] = vyjf[:, i] = np.empty((vxj.shape[0])) * np.nan #np.empty((vxj.shape[1], vxj.shape[1])) * np.nan
        else:
            # VX fitting
            frac = lowess_window / np.sum(~np.isnan(vxj[:, i]))
            if frac > 0.9: frac=0.9 # less data in the pixel than -lowess_window- -> frac > 1
            vxjf[:, i] = sm.nonparametric.lowess(vxj[:, i], dates, frac=frac, it=3,
                                                delta=np.sum(~np.isnan(vxj[:, i])) * 0.01, return_sorted=False)
            while np.sum(~np.isnan(vxjf[:,i]))==0: # frac instability bug fixing
                frac=frac+0.005
                vxjf[:,i] = sm.nonparametric.lowess(vxj[:, i], dates, frac=frac, it=3,
                                               delta=np.sum(~np.isnan(vxj[:, i])) * 0.01, return_sorted=False)
            # VY fitting
            vyjf[:, i] = sm.nonparametric.lowess(vyj[:, i], dates, frac=frac, it=3,
                                           delta=np.sum(~np.isnan(vyj[:, i])) * 0.01, return_sorted=False)
            while np.sum(~np.isnan(vyjf[:, i])) == 0:
                frac = frac + 0.001
                vyjf[:, i] = sm.nonparametric.lowess(vyj[:, i], dates, frac=frac, it=3,
                                                     delta=np.sum(~np.isnan(vyj[:, i])) * 0.01, return_sorted=False)
            # # fix the wright-side outbursts instability
            # vxjf[:, i] = np.where(abs(vxjf[:, i]) > np.nanquantile(abs(vxjf[:, i]), 0.995), np.nan, vxjf[:, i])
            # vyjf[:, i] = np.where(abs(vyjf[:, i]) > np.nanquantile(abs(vyjf[:, i]), 0.995), np.nan, vyjf[:, i])

    return vxjf, vyjf # }}}

def stack2layer(cube2, c_m, func, start, end, source, err2D, weighted, verbose):
    '''Collaps the stack of layers to one by given function (mean, median, std) and insert to the -c_m- cube.  '''

    # MEAN error = STD / sqrt(pix_N)
    # count the pix_N (number of valid pixels)
    N2d = np.zeros([cube2.ny, cube2.nx])
    for ii in range(cube2.nz):
        N2d = N2d + np.float32(~cube2.data[ii].vx.mask)


    if cube2.nz == 0:
        if verbose: print('     empty date')
        if weighted or err2D:  # make 2D error to synchronize processing
            c_m.insert(np.full((cube2.ny, cube2.nx), 0), np.full((cube2.ny, cube2.nx), 0),
                       np.full((cube2.ny, cube2.nx), 0), np.full((cube2.ny, cube2.nx), 0),
                       date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=True, verbose=verbose)
        else:
            c_m.insert(np.full((cube2.ny, cube2.nx), 0), np.full((cube2.ny, cube2.nx), 0),
                       0, 0,
                       date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=False, verbose=verbose)
    elif cube2.nz == 1:
        if func != 'std':
            c_m.insert(cube2.data[0].vx, cube2.data[0].vy,
                       cube2.data[0].errx, cube2.data[0].erry,
                       date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=err2D, verbose=verbose)
        else:
            c_m.insert(np.full((cube2.ny, cube2.nx), 0), np.full((cube2.ny, cube2.nx), 0),
                       #np.full(cube2.data[0].errx.shape, 0), np.full(cube2.data[0].errx.shape, 0),
                       N2d, N2d,
                       date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=err2D, verbose=verbose)
    else:
        if weighted:
            status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = cube2.cube2vel_py()
            if status:
                if func == 'mean':
                    c_m.insert(vxo, vyo, errxo, erryo,
                               date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=True, verbose=verbose)
                elif func == 'median':
                    c_m.insert(medxo, medyo, errxo, erryo,
                               date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=True, verbose=verbose)
                elif func == 'std':
                    c_m.insert(stdxo, stdyo,
                               N2d, N2d, #errxo, erryo,
                               date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=True, verbose=verbose)
        else:
            if func == 'mean': # mean_error ar described on the top in error
                c_m.insert(np.nanmean(cube2.vx_(), axis=0), np.nanmean(cube2.vy_(), axis=0),
                           np.nanstd(cube2.vx_(), axis=0)/np.sqrt(N2d), np.nanstd(cube2.vy_(), axis=0)/np.sqrt(N2d), #np.nanmean(cube2.errx_(), axis=0), np.nanmean(cube2.erry_(), axis=0),
                           date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=err2D, verbose=verbose)
            elif func == 'median': # mean_error ar described on the top in error
                c_m.insert(np.nanmedian(cube2.vx_(), axis=0), np.nanmedian(cube2.vy_(), axis=0),
                           np.nanstd(cube2.vx_(), axis=0) / np.sqrt(N2d), np.nanstd(cube2.vy_(), axis=0) / np.sqrt(N2d), #np.nanmedian(cube2.errx_(), axis=0), np.nanmedian(cube2.erry_(), axis=0),
                           date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=err2D, verbose=verbose)
            elif func == 'std': # valid_pix_N in error
                c_m.insert(np.nanstd(cube2.vx_(), axis=0), np.nanstd(cube2.vy_(), axis=0),
                           N2d, N2d, #np.nanstd(cube2.errx_(), axis=0), np.nanstd(cube2.erry_(), axis=0),
                           date1=date2mjd(start), date2=date2mjd(end), source=source, sensor='', err2D=err2D, verbose=verbose)
    return c_m

def reset_duplicates(cube, verbose):
    ''' Manage the dates replication: fill with median value one layer and fill with nan the others  '''
    dates = np.array(cube.date1_()) + np.array(cube.offset_()) / 2

    if verbose: print('     Manage the dates replication...')
    sort_d_ind = np.argsort(dates)
    dates_eq_ind = [dates[sort_d_ind[i]] == dates[sort_d_ind[i + 1]] for i in range(len(sort_d_ind) - 1)]
    dates_eq_ind.append(False)
    loc_list = []
    for i in range(1, len(dates_eq_ind)):
        if dates_eq_ind[i - 1] == True:
            loc_list.append(sort_d_ind[i - 1]) # collect the list of layers to be avaraged
            if dates_eq_ind[i] == False:
                loc_list.append(sort_d_ind[i])
                # save the medians into the last layer
                cube.data[sort_d_ind[i]].vx = np.nanmedian( cube.vx_()[loc_list].filled(np.nan), 0)
                cube.data[sort_d_ind[i]].vy = np.nanmedian( cube.vy_()[loc_list].filled(np.nan), 0)
                cube.data[sort_d_ind[i]].errx = np.nanmedian( cube.errx_()[loc_list], 0)
                cube.data[sort_d_ind[i]].erry = np.nanmedian(cube.erry_()[loc_list], 0)
                # reset to -nan- all other layers
                for k in loc_list[:-1]:
                    cube.data[k].vx = np.ma.array(np.full([cube.ny, cube.nx], np.nan))
                    cube.data[k].vy = np.ma.array(np.full([cube.ny, cube.nx], np.nan))
                    cube.data[k].errx = np.full(cube.data[sort_d_ind[i]].errx.shape, np.nan)
                    cube.data[k].erry = np.full(cube.data[sort_d_ind[i]].erry.shape, np.nan)
                loc_list = []
    if verbose: print('     ... done')

    return cube

# ====== = ====== OTHERS ====== = ======

def compress_cubes(cube_filename,backup=False):
    c=cube_class()
    c.load(cube_filename)
    c.describe()
    if backup:
        if os.path.exists(cube_filename+'.BAK'):
            os.system('cp '+cube_filename+' '+cube_filename+'.BAK')

    c.write(replace=True,compress=True)
    c=None
