import numpy as np
import netCDF4
from osgeo import osr
import os
import time
from datetime import datetime, timedelta
from fparam import geo_param # /ST_RELEASE/COMMON/PYTHON/fparam.py
from mjd2date import mjd2date, date2mjd # /ST_RELEASE/UTILITIES/PYTHON/mjd2date.py
from progress_bar import progress_bar # # /ST_RELEASE/UTILITIES/PYTHON/progress_bar.py
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message='Warning: converting a masked element to nan.')

class cube_class:

    def __init__(self, nx=250, ny=250, proj='PS'):  

        # !!! CUBE ELEMENTS DIMENSIONS:
        # -Z-Y-X-

        self.filedir = ''
        self.filename = ''
        self.x = np.zeros(nx, dtype=np.float32)
        self.y = np.zeros(ny, dtype=np.float32)
        self.z = 0
        self.i = 0
        self.j = 0
        self.nx = nx
        self.ny = ny
        self.nz = 0
        self.shape = (self.nz, self.ny , self.nx)
        self.proj4 = ''
        self.srs = osr.SpatialReference()
        self.data = []
        self.version='v17jan2022'

        if proj == 'PS' or proj == 'PS_NORTH':
            self.geo = geo_param()
            self.geo.central_meridian = -45.0
            self.geo.secant_lat = 70.0
            self.srs.ImportFromEPSG(3413)
            self.proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
            self.mapping = {'ellipsoid': 'WGS84', 'false_easting': 0.0, 'false_northing': 0.0, \
                            'grid_mapping_name': 'polar_stereographic', 'latitude_of_projection_origin': 90.0, \
                            'standard_parallel': 70.0, 'straight_vertical_longitude_from_pole': -45.0}

        if proj == 'PS_SOUTH':
            self.geo = geo_param()
            self.geo.central_meridian = 0.0
            self.geo.secant_lat = -71.0
            self.srs.ImportFromEPSG(3031)
            self.proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
            self.mapping = {'ellipsoid': 'WGS84', 'false_easting': 0.0, 'false_northing': 0.0, \
                    'grid_mapping_name': 'polar_stereographic', 'latitude_of_projection_origin': -90.0, \
                    'standard_parallel': -71.0, 'straight_vertical_longitude_from_pole': 0.0}

        if proj == 'UTM':
            self.geo = geo_param(utm=True)
            self.mapping = { \
                'grid_mapping_name': 'transverse_mercator', 'spatial_ref': ''}
        self.geo.npix=int(self.geo.npix)
        self.geo.nrec = int(self.geo.nrec)
    
# ====== = ====== FILE INPUT/OUTPUT BASICS  ====== = ======

    def load(self, filepath, concat=False, verbose=False):  # {{{
        '''
        :param filepath: string. ! if relative path is provided (e.g. only a name), cube.filedir will be NOT filled and cube.save will crash
        :param concat: append the layers from the new cube to the current cube. The same space position is required.
        :param pick_month: keep the layers where center of date1-date2 period is in the given interval (including edges), [int,int]
                [6,8] - from June to August, [11,2] - from November to February
        :param sensor: replace in ALL layers -sensor- with given value
        :param verbose: print the intermediate into

        example of use:
        cube.load('path/path/name1.nc')
        cube.load('path/path/name2.nc', concat=True)
        '''
        if verbose:
            print(filepath)

        if not os.path.exists(filepath):
            print('ERROR: does not exists..')
            return False

        nc = netCDF4.Dataset(filepath)
        if not concat:

            if nc.Projection == 'UTM':
                self.__init__(proj='UTM')
            else:
                self.__init__()

            
            if os.path.dirname(filepath) == '':
                self.filedir = './'
            else:
                self.filedir = os.path.dirname(filepath) + '/'
            
            self.filename = os.path.basename(filepath)
            self.proj4 = nc.proj4
            self.x = nc.variables['x'][:]
            self.y = nc.variables['y'][:]

            try:
                nc.dimensions['z'] # does not generate error if the -z- dimension exist
                self.z = [z for z in range(nc.nz)]
            except:
                self.z = 'It looks like ELMER cube. -z- dimention does not exist.'
            
            self.nx = len(self.x)
            self.ny = len(self.y)

            try: # cubes v2.0
                self.i = nc.i
                self.j = nc.j
            except: # cubes v1.0
                try:
                    self.i = np.int32(self.filename.split('.nc')[0].split('_')[1].split('x')[1])
                    self.j = np.int32(self.filename.split('.nc')[0].split('_')[2].split('y')[1])
                except:
                    print('     >>> WARNING: -i-j- are not detected in the file metadata nor in the file name. Some operations will be unavailable. <<<')
                    self.i=-9999
                    self.j=-9999

            if 'ellipsoid' in nc.variables['mapping'].ncattrs():
                self.mapping['ellipsoid'] = nc.variables['mapping'].ellipsoid
                
            if nc.Projection == 'PS' or nc.Projection == 'PS_NORTH' or nc.Projection == 'PS_SOUTH':
                self.mapping['false_easting'] = nc.variables['mapping'].false_easting
                self.mapping['false_northing'] = nc.variables['mapping'].false_northing
                self.mapping['grid_mapping_name'] = nc.variables['mapping'].grid_mapping_name
                self.mapping['latitude_of_projection_origin'] = nc.variables['mapping'].latitude_of_projection_origin
                self.mapping['standard_parallel'] = nc.variables['mapping'].standard_parallel
                self.mapping['straight_vertical_longitude_from_pole'] = nc.variables[
                    'mapping'].straight_vertical_longitude_from_pole
                
            if nc.Projection == 'UTM':
                self.geo = geo_param(utm=True)
                print(self.geo.npix)
                tmp = nc.proj4.split('+')
                EPSG='326'
                for t in tmp: #parse proj4 string, ex: "+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs"
                    if 'zone' in t:
                        zone=t.split('=')[1]
                    if 'south' in t:
                        EPSG='327'
                EPSG=EPSG+zone
                if verbose: print(EPSG)

                or_srs = osr.SpatialReference()
                or_srs.ImportFromEPSG(int(EPSG))

                self.geo.projection_zone = or_srs.GetUTMZone()
                
                self.geo.false_easting = np.float64(or_srs.GetProjParm(osr.SRS_PP_FALSE_EASTING))
                self.geo.false_northing = np.float64(or_srs.GetProjParm(osr.SRS_PP_FALSE_NORTHING))
                self.geo.projection_k0 = np.float64(or_srs.GetProjParm(osr.SRS_PP_SCALE_FACTOR))
                self.geo.center_longitude = np.float64(or_srs.GetProjParm(osr.SRS_PP_CENTRAL_MERIDIAN))

                if self.geo.center_longitude > 180:
                    self.geo.center_longitude = self.geo.center_longitude - 360.

                self.geo.center_latitude = np.float64(or_srs.GetProjParm(osr.SRS_PP_LATITUDE_OF_CENTER))

                self.srs = or_srs
                if verbose: print('self.geo.projection_zone:', self.geo.projection_zone)

                if 'grid_mapping_name' in nc.variables['mapping'].ncattrs():
                    self.mapping['grid_mapping_name'] = nc.variables['mapping'].grid_mapping_name
                else:
                    self.mapping['grid_mapping_name'] = 'universal_transverse_mercator'

                if 'spatial_ref' in nc.variables['mapping'].ncattrs():
                    self.mapping['spatial_ref'] = nc.variables['mapping'].spatial_ref
                else:
                    self.mapping['spatial_ref'] = or_srs.ExportToWkt()
            
            self.geo.xmin = min(self.x)
            self.geo.ymax = max(self.y)
            self.geo.npix = len(self.x)
            self.geo.nrec = len(self.y)
            self.geo.posting = self.x[1] - self.x[0]
            self.geo.xposting = self.x[1] - self.x[0]
            self.geo.yposting = self.y[1] - self.y[0]
            self.geo.projection = nc.Projection

        if nc.variables['error_vx'][0].shape==(self.geo.nrec,self.geo.npix):
            print('2D error map detected !')

        
        vx = nc.variables['vx'][:,:,:]
        vy = nc.variables['vy'][:,:,:]
        
        date1 = list(nc.variables['date1'][:])
        date2 = list(nc.variables['date2'][:])
        
        for z in range(nc.nz):
            if nc.variables['error_vx'][0].shape==(): # error as unique value
                data = cube_element(nx=self.geo.npix, ny=self.geo.nrec, err2D=False)
            elif nc.variables['error_vx'][0].shape==(self.geo.nrec,self.geo.npix): # error as map
                data = cube_element(nx=self.geo.npix, ny=self.geo.nrec, err2D=True)
            else: raise ValueError ("     >>> WARNING: -Error- variable dimensions not recognized! Can't load file. <<< ")

            data.vx = vx[z,:,:]  # masked array type by default, NetCDF -missing_value- attribute is masked.
            data.vy = vy[z,:,:] # /!!!\ To get properly data them-self, do > b=a.filled(0). Or > b=a.filled(np.nan) to have NaN-s
            data.vx.mask = data.vx.mask | (data.vx == 0) | np.isnan(data.vx) # add the compatibility with source-cubes compiled from tifs
            data.vy.mask = data.vy.mask | (data.vy == 0) | np.isnan(data.vy)

            data.errx = nc.variables['error_vx'][z]
            data.erry = nc.variables['error_vy'][z]
            
            data.date1 = date1[z]
            data.date2 = date2[z]
            
            data.offset = data.date2 - data.date1
            
            data.source = str(nc.variables['source'][z].data, 'utf-8')
            
            if 'sensor' in nc.variables:
                data.sensor = str(nc.variables['sensor'][z].data, 'utf-8')
            
            data.source = data.source.strip('\x00')
            
            self.data.append(data)
            self.nz = self.nz + 1

        del vx, vy

        self.shape = (self.nz, self.ny, self.nx)

        nc.close() # }}}

    def write(self, replace=False, compress=True, verbose=False):  # {{{

        '''
        Save the cube variable to the NetCDF file using the -cube.filedir- and -cube.filename-
        If the filename.nc exist and -replace- is False, filename_1.nc is created.
        :param replace: replace if the file with given name exists
        :param compress: compress when saving. reduce the final size a lot, increase the cube.load time a little bit

        example of use:
        cube.write()
        '''

        nz = self.nz
        ny = self.ny
        nx = self.nx

        if verbose: print('    Save cube ', nz, ny, nx)
        
        filepath = os.path.join(self.filedir, self.filename)
        if verbose: print(filepath)
        
        if replace: 
            if os.path.exists(filepath):
                if verbose: print(filepath+'is present and is being updated.')
                os.system('rm ' + filepath) # force re-writing
        else:
            if os.path.exists(filepath): # save filename_1.nc
                print(filepath+' is present. To overwrite, use write(replace=True)')
                return

        
        dataset = netCDF4.Dataset(os.path.join(self.filedir, self.filename), 'w', format='NETCDF4')
        # create internal structure
        dataset.Conventions = 'CF-1.6'
        dataset.Title = 'Cube of Ice Velocity'
        dataset.Version = '2.0 (Dec2019)' # 0 is stored in NetCDF as -missing_value- and -fill_value- by default
        dataset.Author = 'J. Mouginot, R.Millan, A.Derkacheva'
        dataset.history = 'Created ' + time.ctime(time.time())
        dataset.Notes = 'Data were processed at the Department of Earth System Science, University of California, Irvine, USA and at the Institut des Geosciences de Environnement, Universite Grenoble Alpes, France '
        dataset.nx = self.nx
        dataset.ny = self.ny
        dataset.nz = self.nz
        dataset.i = self.i
        dataset.j = self.j

        dataset.createDimension('x', nx)
        dataset.createDimension('y', ny)
        dataset.createDimension('z', nz)
        dataset.createDimension('nletter', 300)
        dataset.createDimension('nletter2', 30)

        mapping = dataset.createVariable('mapping', 'S1')
        x = dataset.createVariable('x', 'f4', 'x')
        y = dataset.createVariable('y', 'f4', 'y')
        vx = dataset.createVariable('vx', 'f4', ('z', 'y', 'x'), fill_value=0, zlib=compress, least_significant_digit=2,complevel=3) #vx = dataset.createVariable('vx', 'f4', ('z', 'y', 'x'))
        vy = dataset.createVariable('vy', 'f4', ('z', 'y', 'x'), fill_value=0, zlib=compress, least_significant_digit=2,complevel=3)

        date1 = dataset.createVariable('date1', 'i4', 'z')
        date2 = dataset.createVariable('date2', 'i4', 'z')

        source = dataset.createVariable('source', 'S1', ('z', 'nletter'), zlib=compress,complevel=3)
        sensor = dataset.createVariable('sensor', 'S1', ('z', 'nletter2'), zlib=compress,complevel=3)
        creation_date = dataset.createVariable('creation_date', 'S1', ('z', 'nletter'), zlib=compress)
        if self.data[0].errx.shape==(): # error as unique value
            error_vx = dataset.createVariable('error_vx', 'f4', 'z', fill_value=0)
            error_vy = dataset.createVariable('error_vy', 'f4', 'z', fill_value=0)
        else: # error as map
            error_vx = dataset.createVariable('error_vx', 'f4', ('z', 'y', 'x'), fill_value=0,zlib=compress,complevel=3)
            error_vy = dataset.createVariable('error_vy', 'f4', ('z', 'y', 'x'), fill_value=0,zlib=compress,complevel=3)

        dataset.Projection = self.geo.projection
        dataset.proj4 = self.proj4

        if self.mapping['grid_mapping_name'] == 'polar_stereographic':
            mapping.grid_mapping_name = 'polar_stereographic'
            mapping.setncattr_string('ellipsoid', self.mapping['ellipsoid'])
            mapping.false_easting = self.mapping['false_easting']
            mapping.false_northing = self.mapping['false_northing']
            mapping.setncattr_string('grid_mapping_name', self.mapping['grid_mapping_name'])
            mapping.latitude_of_projection_origin = self.mapping['latitude_of_projection_origin']
            mapping.standard_parallel = self.mapping['standard_parallel']
            mapping.straight_vertical_longitude_from_pole = self.mapping['straight_vertical_longitude_from_pole']

        if self.mapping['grid_mapping_name'] == 'transverse_mercator' or self.mapping['grid_mapping_name'] == 'universal_transverse_mercator':
            mapping.grid_mapping_name = 'universal_transverse_mercator'
            mapping.utm_zone_number = self.geo.projection_zone
            mapping._CoordinateTransformType = "Projection"
            mapping._CoordinateAxisTypes = "GeoX GeoY"
            mapping.spatial_ref = self.srs.ExportToWkt()

        vx.setncattr_string('long_name', 'Ice velocity in x direction')
        vx.setncattr_string('standard_name', 'land_ice_x_velocity')
        vx.units = 'meter/year'
        vx.setncattr_string('grid_mapping', 'mapping')
        vx.missing_value = 0

        vy.setncattr_string('long_name', 'Ice velocity in y direction')
        vy.setncattr_string('standard_name', 'land_ice_y_velocity')
        vy.units = 'meter/year'
        vy.setncattr_string('grid_mapping', 'mapping')
        vy.missing_value = 0

        error_vx.setncattr_string('long_name', 'Estimated error for ice velocity in x direction')
        error_vx.setncattr_string('standard_name', 'vx_velocity_error')
        error_vx.missing_value = 0
        if self.data[0].errx.shape != ():
            error_vx.setncattr_string('grid_mapping', 'mapping')

        error_vy.setncattr_string('long_name', 'Estimated error for ice velocity in y direction')
        error_vy.setncattr_string('standard_name', 'vy_velocity_error')
        error_vy.missing_value = 0
        if self.data[0].erry.shape != ():
            error_vy.setncattr_string('grid_mapping', 'mapping')

        x.setncattr_string('long_name', 'Cartesian x-coordinate')
        x.setncattr_string('standard_name', 'projection_x_coordinate')
        x.units = 'm'
        x.axis = "X"

        y.setncattr_string('long_name', 'Cartesian y-coordinate')
        y.setncattr_string('standard_name', 'projection_y_coordinate')
        y.units = 'm'
        y.axis = "Y"

        date1.content = 'Date of the first acquisition'
        date1.units = 'day (Modified Julian Date)'
        date1.description = 'conversion to datetime in python : t = jdcal.jd2gcal(2400000.5, date); return datetime.date(t[0], t[1], t[2])'
        date2.content = 'Date of the second acquisition'
        date2.units = 'day (Modified Julian Date)'
        date1.description = 'conversion to datetime in python : t = jdcal.jd2gcal(2400000.5, date); return datetime.date(t[0], t[1], t[2])'

        creation_date.content = 'Date of last file modif. (update purposes)'

        # add data into the NetCDF
        x[:] = np.float32(self.x)
        y[:] = np.float32(self.y)
        vx[:,:,:] = self.vx_()
        vy[:,:,:] = self.vy_()
           
        for z in range(self.nz):
            
            if verbose: 
                if z % 10 == 0:progress_bar( z/self.nz*100. , msg = str(int(z/self.nz*100.))+'% layers saved in nc')
            
            date1[z] = min([self.data[z].date1, self.data[z].date2])
            date2[z] = max([self.data[z].date1, self.data[z].date2])
            
            error_vx[z] = self.data[z].errx
            error_vy[z] = self.data[z].erry
            
            tmp = netCDF4.stringtoarr(self.data[z].source, len(self.data[z].source))
            source[z, 0:len(tmp)] = netCDF4.stringtoarr(self.data[z].source, len(self.data[z].source))
            
            tmp = netCDF4.stringtoarr(self.data[z].sensor, len(self.data[z].sensor))
            sensor[z, 0:30] = netCDF4.stringtoarr(self.data[z].sensor[0:30], 30)
            # tmp = netCDF4.stringtoarr(self.data[i].creation_date, len(self.data[i].creation_date))
            # creation_date[i,0:len(tmp)] = netCDF4.stringtoarr(self.data[i].creation_date, len(self.data[i].creation_date))

        dataset.close()
        if verbose: print('     {0} data layers successfully saved in {1}/{2}'.format(len(self.data), self.filedir, self.filename))

    def describe(self):  # {{{
        ''' Print information about the cubes '''
        print('\n\n Cube description \n')
        print('Filename :', self.filename)
        print('Directory :',self.filedir)
        print('i, j (pixel position in full grid) :', self.i, self.j)
        print('nx, ny : ', self.nx, self.ny)
        print('Number of layers in data : ', self.nz)
        
        print('\n - Geographic description - ')
        self.geo.describe()
        print('Projection (proj4) : ', self.proj4)
        
        print('Overlap between adjacent cube is always 5 pixels.') # }}}

    def merge(self, cube2):  # {{{
        '''
        Merge spatially two cubes that have adjacent space position.
        :param cube2: path to the NetCDF file (str) or cube_class variable

        example of use:
        cube = cube.merge('path/path/filename2.nc')
        cube = cube.merge(cube2)
        '''

        if type(cube2)==str:
            cube2_path = cube2
            cube2 = cube_class()
            cube2.load(cube2_path)

        err_trygger=False
        if (self.data[0].errx.shape == () and cube2.data[0].errx.shape != ()) :
            print ("    WARNING! You try to merge the cubes with different errors' dimensions (unique value vs 2D map).")
            print ("    The both will be forced to the 2D map")
            err_trygger = 1
        elif (self.data[0].errx.shape != () and cube2.data[0].errx.shape == ()):
            print ("    WARNING! You try to merge the cubes with different errors' dimensions (unique value vs 2D map).")
            print ("    The both will be forced to the 2D map")
            err_trygger = 2
        elif (self.data[0].errx.shape != () and cube2.data[0].errx.shape != ()) :
            err_trygger = 3

        nx = self.nx
        ny = self.ny
        nx0 = self.nx
        ny0 = self.ny
        x = self.x
        y = self.y
        z = self.z
        # check if adjacent cube:
        if cube2.y[0] >= self.y[ny - 5] and cube2.y[cube2.ny - 5] <= self.y[0]:
            if cube2.x[0] <= self.x[nx - 5] and cube2.x[cube2.nx - 5] >= self.x[0]:

                # should I increase cube size ?
                if cube2.y[0] == self.y[ny - 5]:
                    y = np.ma.append(self.y[0:ny - 5], cube2.y[:])
                    ny = ny + cube2.ny - 5
                if cube2.y[cube2.ny - 5] == self.y[0]:
                    y = np.ma.append(cube2.y[0:cube2.ny - 5], self.y[:]) #y = np.ma.append(cube2.y[0:ny - 5], self.y[:])
                    ny = ny + cube2.ny - 5

                if cube2.x[0] == self.x[nx - 5]:
                    x = np.ma.append(self.x[0:nx - 5], cube2.x[:])
                    nx = nx + cube2.nx - 5
                if cube2.x[cube2.nx - 5] == self.x[0]:
                    x = np.ma.append(cube2.x[0:nx - 5], self.x[:])
                    nx = nx + cube2.nx - 5

                big = cube_class(nx=nx, ny=ny)
                big.i = cube2.i if cube2.x[0] < self.x[0] else self.i
                big.j = cube2.j if cube2.y[0] > self.y[0] else self.j  ## y pixel are inverse !
                big.filedir = self.filedir
                big.filename = 'bigc_x{:05}_y{:05}_'.format(big.i, big.j) + '_'.join(self.filename.split('_')[3:])
                big.proj4 = self.proj4
                big.srs = self.srs
                big.geo = self.geo
                big.geo.xmin = cube2.geo.xmin if cube2.geo.xmin < self.geo.xmin else self.geo.xmin
                big.geo.ymax = cube2.geo.ymax if cube2.geo.ymax > self.geo.ymax else self.geo.ymax
                big.geo.npix = nx
                big.geo.nrec = ny
                big.nx = nx
                big.ny = ny
                big.x = x
                big.y = y
                big.z = z
                if cube2.mapping['grid_mapping_name']=='transverse_mercator':
                    big.mapping['grid_mapping_name'] = 'transverse_mercator'
                    big.mapping['spatial_ref'] = cube2.srs.ExportToWkt()

                for cube_el in self.data:
                    if cube_el.errx.shape == () and not err_trygger :  # error as unique value
                        big_cube_element = cube_element(nx=nx, ny=ny, err2D=False)
                    else:  # error as map
                        big_cube_element = cube_element(nx=nx, ny=ny, err2D=True)
                    big_cube_element.vx[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.vx
                    big_cube_element.vy[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.vy
                    big_cube_element.source = cube_el.source
                    big_cube_element.sensor = cube_el.sensor
                    if err_trygger==1: # force unique value to the 2D map
                        big_cube_element.errx = np.full([ny,nx],cube_el.errx,dtype=np.float32)
                        big_cube_element.erry = np.full([ny,nx],cube_el.erry,dtype=np.float32)
                    elif err_trygger==2 or err_trygger==3:
                        big_cube_element.errx[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.errx
                        big_cube_element.erry[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.erry
                    else: # err_trygger=False, both cubes have 1D-errors
                        big_cube_element.errx = cube_el.errx
                        big_cube_element.erry = cube_el.erry
                    big_cube_element.date1 = cube_el.date1
                    big_cube_element.date2 = cube_el.date2
                    big_cube_element.offset = cube_el.date2 - cube_el.date1
                    big.data.append(big_cube_element)

                for cube_el in cube2.data:
                    # the same source = the pieces of the same image
                    if (cube_el.source in big.source_()) and (cube_el.source.count('\x00') != 300):
                        z = big.source_().index(cube_el.source)
                        big.data[z].vx[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.vx
                        big.data[z].vy[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.vy
                        if cube_el.errx.shape!=() : # source error as 2D
                            big.data[z].errx[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.errx
                            big.data[z].erry[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.erry
                        elif err_trygger==2: # source error as 1D, output error as 2D
                            big_cube_element.errx[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = np.full([ny, nx], cube_el.errx, dtype=np.float32)
                            big_cube_element.erry[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = np.full([ny, nx], cube_el.erry, dtype=np.float32)
                        # else: 1D error kept from the self-cube

                    # image which not yet in cube
                    else:
                        if cube_el.errx.shape==() and not err_trygger:  # source and output errors as unique value
                            big_cube_element = cube_element(nx=nx, ny=ny, err2D=False)
                            big_cube_element.errx = cube_el.errx
                            big_cube_element.erry = cube_el.erry
                        elif err_trygger==2: # source error as 1D, output error as 2D
                            big_cube_element = cube_element(nx=nx, ny=ny, err2D=True)
                            big_cube_element.errx[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = np.full([ny, nx], cube_el.errx, dtype=np.float32)
                            big_cube_element.erry[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = np.full([ny, nx], cube_el.erry, dtype=np.float32)
                        else:  # source error as 2D map
                            big_cube_element = cube_element(nx=nx, ny=ny, err2D=True)
                            big_cube_element.errx[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = cube_el.errx
                            big_cube_element.erry[cube2.j - big.j:cube2.j - big.j + 250, cube2.i - big.i:cube2.i - big.i + 250] = cube_el.erry
                        big_cube_element.vx[cube2.j - big.j:cube2.j - big.j + 250,
                        cube2.i - big.i:cube2.i - big.i + 250] = cube_el.vx
                        big_cube_element.vy[cube2.j - big.j:cube2.j - big.j + 250,
                        cube2.i - big.i:cube2.i - big.i + 250] = cube_el.vy
                        big_cube_element.source = cube_el.source
                        big_cube_element.sensor = cube_el.sensor
                        big_cube_element.date1 = cube_el.date1
                        big_cube_element.date2 = cube_el.date2
                        big_cube_element.offset = cube_el.date2 - cube_el.date1
                        big.data.append(big_cube_element)

                big.nz = len(big.data)
                print('Merged successfully. Current cube.filename ', big.filename)
                print('Current big_cube size nz-ny-nx: ', big.nz, big.ny, big.nx)
                self = big
                del big

            else: print('   Cubes are not adjacent. No successful merge could be done.')
        else: print('   Cubes are not adjacent. No successful merge could be done.')

        return self

    def insert(self, vx, vy, errx, erry, date1, date2, source, sensor, err2D=False, replace=False, verbose=True):
        ''' Insert non-cubic data into the cube.
            Data could be as one layer pack (z=1) or multy-layers lists/arrays
         :param vx, vy: velocity with variable dimensions -z-y-x-,  numpy array
         :param errx, erry: error as one value per layer -z- or 2D map -z-x-y- (set err2D=True), numpy array
         :param date1, date2: dates in MJD, integer or -z- array
         :param source, sensor: text descriptions, str or -z- list
         :param err2D : check True in the error representation is 2D map
         :param replace : check it to rewrite cube from zero, otherwise the data will be appended

         example of use: insert 3 layers filled with 0
         cube.insert( vx = np.full((3, cube2.ny, cube2.nx), 0), vy = np.full((3, cube2.ny, cube2.nx), 0),
                      errx = [0, 0, 0], erry = [0, 0, 0], date1=[110, 220, 330], date2=[111, 221, 331],
                      source = ['text1', 'text2', 'text3'], sensor=['', '', ''], err2D=False, verbose=True)
         '''

        if replace: # delete all source data
            while len(self.data) > 0:
                self.data.pop(0)
            self.nz = 0

        try: # arrays
            zz=len(date1)
            for z in range(len(date1)):
                # insert data
                data = cube_element(nx=self.geo.npix, ny=self.geo.nrec, err2D=err2D)
                try:
                    data.vx = np.ma.masked_invalid(vx[z,:,:])  # data are stored as numpy masked array
                    data.vy = np.ma.masked_invalid(vy[z,:,:])
                    if err2D:
                        data.errx = np.ma.masked_invalid(errx[z,:,:])
                        data.erry = np.ma.masked_invalid(erry[z,:,:])
                    else:
                        data.errx = np.ma.masked_invalid(errx[z])
                        data.erry = np.ma.masked_invalid(erry[z])
                except: raise ValueError(' ---> Velocity or errors are in unacceptable format (array case).')
                try:
                    data.date1 = date1[z]
                    data.date2 = date2[z]
                    data.offset = date2[z] - date1[z]
                except: raise ValueError(' ---> Dates or offset are in unacceptable format (array case).')
                try:
                    data.source = source[z]
                    data.source = data.source.strip('\x00')
                    data.sensor = sensor[z]
                    data.sensor = data.sensor.strip('\x00')
                except: raise ValueError(' ---> Source or sensor are in unacceptable format (array case).')
                # save to element to the cube
                self.data.append(data)
                self.nz = self.nz + 1
        except: # one layer
            data = cube_element(nx=self.geo.npix, ny=self.geo.nrec, err2D=err2D)
            try:
                data.vx = np.ma.masked_invalid(vx[:, :])  # data are stored as numpy masked array
                data.vy = np.ma.masked_invalid(vy[:, :])
                if err2D:
                    data.errx = np.ma.masked_invalid(errx[:, :])
                    data.erry = np.ma.masked_invalid(erry[:, :])
                else:
                    data.errx = np.ma.masked_invalid(errx)
                    data.erry = np.ma.masked_invalid(erry)
            except:raise ValueError(' ---> Velocity or errors are in unacceptable format (array case).')
            try:
                data.date1 = date1
                data.date2 = date2
                data.offset = date2-date1
            except:raise ValueError(' ---> Dates or offset are in unacceptable format (array case).')
            try:
                data.source = source
                data.source = data.source.strip('\x00')
                data.sensor = sensor
                data.sensor = data.sensor.strip('\x00')
            except:raise ValueError(' ---> Source or sensor are in unacceptable format (array case).')
            # save to element to the cube
            self.data.append(data)
            self.nz = self.nz + 1
        if verbose: print('     cube.nz = ', self.nz)

    def update(self, c2filepath, update=True, verbose=True): # {{{
        ''' Update cube1's layers with layers from -c2filepath- :
        if update=True:
            if the same source detected, values in the cube1 are replaced with cube2 values;
            if the cube2's layer not yet in the cube1, append it
        if update=False:
            only append the new layers, don't rewrite the existing data

        example of use:
        cube.update('/dir/dir/name.nc', False, False)   '''

        c2=cube_class()
        c2.load(c2filepath)

        # list for source matching: Track is in source
        # list1 = ['/'.join([s for s in self.source_()[z].split('/') if 'track' in s.lower()]) for z in range(self.nz)]
        # list2 = ['/'.join([s for s in c2.source_()[z].split('/') if 'track' in s.lower()]) for z in range(c2.nz)]
        list1 = ['/'.join(s.split('/')[-4:]) for s in self.source_()]
        list2 = ['/'.join(s.split('/')[-4:]) for s in c2.source_()]
        update_list = set(list1) & set(list2)
        # list for source matching: post-processed cubes case
        if len(list1)==0: update_list = set(self.source_()) & set(c2.source_())

        if len(update_list)==0: self.load(c2filepath, concat=True) # no layers to update
        else:
            for z2, source in enumerate(list2):
                #if '/'.join([s for s in source.split('/') if 'track' in s.lower()]) in update_list:
                if '/'.join(source.split('/')[-4:]) in update_list and update: # update an existing layer
                    z1=list1.index(source)
                    if verbose: print(' update: z=', z1, self.data[z1].source)
                    self.data[z1].vx = c2.data[z2].vx
                    self.data[z1].vy = c2.data[z2].vy
                    self.data[z1].errx = c2.data[z2].errx
                    self.data[z1].erry = c2.data[z2].erry
                    self.data[z1].source = c2.data[z2].source
                else: # append new layer
                    if verbose: print(' append: ', c2.data[z2].source)
                    self.insert(c2.data[z2].vx, c2.data[z2].vy, c2.data[z2].errx, c2.data[z2].erry,
                                c2.data[z2].date1, c2.data[z2].date2, c2.data[z2].source, c2.data[z2].sensor,
                                err2D=c2.data[z2].erry.shape != ()) # }}}

    def copy_structure(self,meta_only=False): # {{{
        '''
        Get a new cube from source cube with only the internal structure without data (such as nx-ny, x-y, geo, proj, etc)
        -z- dimension of the new cube will be 0.

        example of use:
        cube2 = cube.copy_structure()
        '''

        newcube = cube_class()
        newcube.filedir = self.filedir
        newcube.filename = self.filename
        newcube.x = self.x
        newcube.y = self.y
        newcube.z = self.z
        newcube.i = self.i
        newcube.j = self.j
        newcube.nx = self.nx
        newcube.ny = self.ny
        newcube.nz = 0
        newcube.proj4 = self.proj4
        newcube.srs = self.srs
        newcube.geo = self.geo
        return newcube

    def source_(self):
        return [d.source for d in self.data]
    def sensor_(self):
        return [d.sensor for d in self.data]
    def color_per_sensor_(self):
        sensors = self.sensor_()
        sources = self.source_()
        colors = []
        labels = []
        for i, sensor in enumerate(sensors):
            if sensor == 'S1A IW IW1 HH' or sensor == 'S1B IW IW1 HH' or sensor == 'S1BIWIW1HH' or sensor == 'S1AIWIW1HH' or sensor == 'S1' or sensor == 'S1A' or sensor == 'S1B':
                colors.append('green')
                labels.append('S1')
            elif sensor == 'L8' or sensor == 'landsat-8':
                colors.append('cornflowerblue')
                labels.append('L8')
            elif sensor == 'Sentinel-2' or sensor == 'SENTINEL-2' or sensor == 'S2':
                colors.append('purple')
                labels.append('S2')
            else:
                if 'LANDSAT' in sources[i]:
                    colors.append('cornflowerblue')
                    labels.append('L8')
                elif 'SENTINEL1' in sources[i]:
                    colors.append('green')
                    labels.append('S1')
                elif 'ERS' in  sources[i]:
                    colors.append('red')
                    labels.append('ERS')
                elif 'RADARSAT' in sources[i]:
                    colors.append('yellow')
                    labels.append('RADARSAT')
                elif 'alos' in sources[i] or 'ALOS' in sources[i]:
                    colors.append('orange')
                    labels.append('ALOS')
                elif 'envisat' in sources[i] or 'ENVISAT' in sources[i]:
                    colors.append('pink')
                    labels.append('ENVISAT')
                elif 'SENTINEL-2' in sources[i] or 'SENTINEL2' in sources[i]:
                    colors.append('purple')
                    labels.append('S2')
                else:
                    colors.append('black')
                    labels.append('OTHERS')
                    print(sources[i])

        return colors, labels

    def vx_(self, i=None, j=None):
        ''' Get -vx- data for the entire cube (i=None, j=None), or for a z-slice (given i or j), or for a pixel (given i and j)'''
        # Dimensions are z, y, x
        if j is not None:
            if j > self.ny: j = None
        if i is not None:
            if i > self.nx: i = None

        if i is None and j is None:
            vx = np.ma.zeros((self.nz, self.ny, self.nx))
            for z, d in enumerate(self.data): vx[z, :, :] = d.vx
        elif i is None and j is not None:
            vx = np.ma.zeros((self.nz, 1, self.nx))
            for z, d in enumerate(self.data): vx[z, 0, :] = d.vx[j, :]
        elif i is not None and j is None:
            vx = np.ma.zeros((self.nz, self.ny, 1))
            for z, d in enumerate(self.data): vx[z, :, 0] = d.vx[:, i]
        elif i is not None and j is not None:
            vx = np.ma.zeros((self.nz))
            for z, d in enumerate(self.data): vx[z] = d.vx[j, i]

        vx.set_fill_value(0)
        return vx
    def vy_(self, i=None, j=None):
        ''' Get -vy- data for the entire cube (i=None, j=None), or for a z-slice (given i or j), or for a pixel (given i and j)'''
        # Dimensions are z, y, x

        if j is not None:
            if j > self.ny: j = None
        if i is not None:
            if i > self.nx: i = None

        if i is None and j is None:
            vy = np.ma.zeros((self.nz, self.ny, self.nx))
            for z, d in enumerate(self.data): vy[z, :, :] = d.vy
        elif i is None and j is not None:
            vy = np.ma.zeros((self.nz, 1, self.nx))
            for z, d in enumerate(self.data): vy[z, 0, :] = d.vy[j, :]
        elif i is not None and j is None:
            vy = np.ma.zeros((self.nz, self.ny, 1))
            for z, d in enumerate(self.data): vy[z, :, 0] = d.vy[:, i]
        elif i is not None and j is not None:
            vy = np.ma.zeros((self.nz))
            for z, d in enumerate(self.data): vy[z] = d.vy[j, i]

        vy.set_fill_value(0)
        return vy
    def date1_(self,datetime=False):
        if datetime:
            from mjd2date import mjd2date
            return [mjd2date(d.date1) for d in self.data]
        else:
            return [d.date1 for d in self.data]
    def date2_(self,datetime=False):
        if datetime:
            from mjd2date import mjd2date
            return [mjd2date(d.date2) for d in self.data]
        else:
            return [d.date2 for d in self.data]
    def date_(self,datetime=False):
        '''
            return central date as modified julian date or datetime
        '''
        if not datetime:
            date1 = np.array( self.date1_())
            date2 = np.array( self.date2_())
        else:
            date1 = np.array( self.date1_(datetime=True))
            date2 = np.array( self.date2_(datetime=True))

        offset = (date2 - date1) / 2
        return date2 - offset
    def offset_(self):
        return [d.date2 - d.date1 for d in self.data]
    def errx_(self, i=None, j=None): #{{{
        ''' Error as unique-value-per-layer case: get -errx- data for all layer, i-j will be ignored.
        2D-map error case: Get -errx- data for the entire cube (i=None, j=None), or for a z-slice (given i or j), or for a pixel (given i and j)'''
        # Dimensions are z, y, x
        if self.data[0].errx.shape == ():
            return np.array([d.errx for d in self.data])
        else:
            if j is not None:
                if j > self.ny: j = None
            if i is not None:
                if i > self.nx: i = None

            if i is None and j is None:
                errx = np.ma.zeros((self.nz, self.ny, self.nx))
                for z, d in enumerate(self.data): errx[z, :, :] = d.errx
            elif i is None and j is not None:
                errx = np.ma.zeros((self.nz, 1, self.nx))
                for z, d in enumerate(self.data): errx[z, 0, :] = d.errx[j, :]
            elif i is not None and j is None:
                errx = np.ma.zeros((self.nz, self.ny, 1))
                for z, d in enumerate(self.data): errx[z, :, 0] = d.errx[:, i]
            elif i is not None and j is not None:
                errx = np.ma.zeros((self.nz))
                for z, d in enumerate(self.data): errx[z] = d.errx[j, i]
            return errx #}}}
    def erry_(self, i=None, j=None): #{{{
        ''' Error as unique-value-per-layer case: get -erry- data for all layer, i-j will be ignored.
        2D-map error case: Get -erry- data for the entire cube (i=None, j=None), or for a z-slice (given i or j), or for a pixel (given i and j)'''
        # Dimensions are z, y, x
        if self.data[0].erry.shape == ():
            return np.array([d.erry for d in self.data])
        else:
            if j is not None:
                if j > self.ny: j = None
            if i is not None:
                if i > self.nx: i = None

            if i is None and j is None:
                erry = np.ma.zeros((self.nz, self.ny, self.nx))
                for z, d in enumerate(self.data): erry[z, :, :] = d.erry
            elif i is None and j is not None:
                erry = np.ma.zeros((self.nz, 1, self.nx))
                for z, d in enumerate(self.data): erry[z, 0, :] = d.erry[j, :]
            elif i is not None and j is None:
                erry = np.ma.zeros((self.nz, self.ny, 1))
                for z, d in enumerate(self.data): erry[z, :, 0] = d.erry[:, i]
            elif i is not None and j is not None:
                erry = np.ma.zeros((self.nz))
                for z, d in enumerate(self.data): erry[z] = d.erry[j, i]
            return erry  # }}}

# ====== = ====== REDUCE DATA ====== = ======
    # If there are two methods -cube.meth- and -cube.meth2-,
    # the 1st-one will change the source -cube- and the 2nd-one will return a new variable -cube2-

    def remove_doublons(self):

        all_vx_sum=[]
        source=[]
        d2=[]
        for d in self.data:
            sumx = np.sum(d.vx)
            if (sumx in all_vx_sum) or (d.source in source):
                print('Found doublons :', d.source)
            else:
                d2.append(d)
                
            source.append(d.source)
            all_vx_sum.append(sumx)
        
        self.data = d2
        self.nz = len(d2)

    def pick_time3m(self, startday, endday, timebuffer, verbose=False):  # {{{
        '''Clip the z dimension according to the given dates (edges are included, -central_date- is checked per layer),
        with buffer dates (date1 > startday-buffertime and date2 < endday+buffertime)
        pick_time3m creates a new vube
        :param startyear: clip from date_central='dd-mm-yyyy' or datetime.date variable
        :param endyear: clip to date_central='dd-mm-yyyy' or datetime.date variable

        example of use:
        cube.pick_time3m('01-01-2015', '31-12-2019')
        '''
        startlimit = mjd2date(date2mjd(startday) - timebuffer)
        endlimit = mjd2date(date2mjd(endday) + timebuffer)


                #if not (startday <= mjd2date(self.data[i].date1 + self.data[i].offset / 2) <= endday) and (startlimit <= mjd2date(self.data[i].date1)) and (mjd2date(self.data[i].date2) <= endlimit):
                    # self.data.pop(i)
                    
        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('Selecting data from {0} to {1} ..'.format(startday, endday))

        if type(startday)==str:
            try:
                startday = datetime.strptime(startday, '%d-%m-%Y').date()
                endday = datetime.strptime(endday, '%d-%m-%Y').date()
            except:
                startday = datetime.strptime(startday, '%Y-%m-%d').date()
                endday = datetime.strptime(endday, '%Y-%m-%d').date()

        newcube = self.copy_structure()
        newcube.filedir = self.filedir
        newcube.filename = self.filename.replace('.nc', '_' + startday.strftime('%d%m%Y') + 'to' + endday.strftime('%d%m%Y') + '.nc')

        for d in self.data:
            if (startday <= mjd2date(d.date1 + d.offset / 2) <= endday) and (startlimit <= mjd2date(d.date1)) and (mjd2date(d.date2) <= endlimit):
                newcube.data.append(d)
        
        newcube.nz = len(newcube.data)

        print('Final cube z-size: ',len(newcube.data))
        return newcube #}}}
              
    def pick_time_(self, startday, endday, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given dates (edges are included, -central_date- is checked per layer).
        pick_time_ return a NEW cube variable, pick_time modifies the original one
        :param startyear: clip from date_central = 'dd-mm-yyyy' or datetime.date variable
        :param endyear: clip to date_central = 'dd-mm-yyyy' or datetime.date variable

        example of use:
        cube2 = cube1.pick_time2('01-01-2015', '31-12-2019')
        '''

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('Selecting data from {0} to {1} ..'.format(startday, endday))

        if type(startday)==str:
            try:
                startday = datetime.strptime(startday, '%d-%m-%Y').date()
                endday = datetime.strptime(endday, '%d-%m-%Y').date()
            except:
                startday = datetime.strptime(startday, '%Y-%m-%d').date()
                endday = datetime.strptime(endday, '%Y-%m-%d').date()

        newcube = self.copy_structure()
        newcube.filedir = self.filedir
        newcube.filename = self.filename.replace('.nc', '_' + startday.strftime('%d%m%Y') + 'to' + endday.strftime('%d%m%Y') + '.nc')

        for d in self.data:
            if (startday <= mjd2date(d.date1 + d.offset / 2) <= endday):
                newcube.data.append(d)
        
        newcube.nz = len(newcube.data)

        if verbose: print('Final cube z-size: ',len(newcube.data))
        return newcube
    def pick_time(self, startday, endday, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given dates (edges are included, -central_date- is checked per layer).
        pick_time modifies the original one, pick_time2 return a NEW cube variable
        :param startyear: clip from date_central='dd-mm-yyyy' or datetime.date variable
        :param endyear: clip to date_central='dd-mm-yyyy' or datetime.date variable

        example of use:
        cube.pick_time('01-01-2015', '31-12-2019')
        '''

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('{0} : {1} data layers picking...'.format(startday, endday))

        if type(startday)==str:
            try:
                startday = datetime.strptime(startday, '%d-%m-%Y').date()
                endday = datetime.strptime(endday, '%d-%m-%Y').date()
            except:
                startday = datetime.strptime(startday, '%Y-%m-%d').date()
                endday = datetime.strptime(endday, '%Y-%m-%d').date()

        i = 0
        while i < len(self.data):
            try:
                if not (startday <= mjd2date(self.data[i].date1 + self.data[i].offset / 2) <= endday):
                    self.data.pop(i)
                else:
                    i = i + 1
            except:
                print('No layers in this cube inside of time period {0} : {1}'.format(startday, endday))
                break
        print(len(self.data), ' layers were picked.')

        self.nz = len(self.data)
        self.filename = self.filename.replace('.nc', '_' + startday.strftime('%d%m%Y') + 'to' + endday.strftime('%d%m%Y') + '.nc') # }}}

    def pick_offset_(self, offset_MIN, offset_MAX, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given ofsets (edges are included).
        :param offset_MIN: integer
        :param offset_MAX: integer

        example of use:
        cube.pick_offset(5, 32)
        '''

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('{0}-{1} days offset layers picking...'.format(offset_MIN, offset_MAX))
        
        newcube = self.copy_structure()
        
        #data = np.copy(self.data)
        temporal_baseline = self.offset_()
        data = [d for i, d in self.data if (offset_MIN <= temporal_baseline[i] <= offset_MAX)]

        if verbose:
            if len(data) == 0:
                print('No layers in this cube with time offset {0}-{1} days'.format(offset_MIN, offset_MAX))
            else:
                print(len(data), ' layers were picked.')

        newcube.nz = len(data)
        newcube.data = data
        data=None

        return newcube
    def pick_offset(self, offset_MIN, offset_MAX, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given ofsets (edges are included).
        pick_time modifies the original one, pick_time2 return a NEW cube variable
        :param offset_MIN: integer
        :param offset_MAX: integer

        example of use:
        cube.pick_offset(5, 32)
        '''

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('{0}-{1} days offset layers picking...'.format(offset_MIN, offset_MAX))

        i = 0
        while i < len(self.data):
            if not (offset_MIN <= self.offset_()[i] <= offset_MAX):
                self.data.pop(i)
            else:
                i = i + 1

        if verbose:
            if i==0:
                print('No layers in this cube with time offset {0}-{1} days'.format(offset_MIN, offset_MAX))
            else:
                print(len(self.data), ' layers were picked.')

        self.nz = len(self.data)
            
    def pick_sensor(self, sensor, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given ofsets (edges are included).
        pick_sensor modifies the original one, pick_sensor2 return a NEW cube variable
        :param sensor: name of sensor, str. It should match any of spellings in -sensors_colormap- dictionary

        example of use:
        cube.pick_sensor('SENTINEL-1')
        '''

        from velocity_colormap import sensors_colormap # manage the different spellings, e.g. "S1" and "Sentinel-1"

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('Sensor=={0} layers picking...'.format(sensor))

        SENSORS = sensors_colormap()
        for sl in SENSORS.items():
            if sensor in sl[1]:
                key = sl[0]
                break
        i = 0
        while i < len(self.data):
            try:
                if not (self.data[i].sensor in SENSORS[key]):
                    self.data.pop(i)
                else:
                    i = i + 1
            except:
                print('No layers in this cube with time sensor=={0}'.format(sensor))
                break
        print(len(self.data), ' layers were picked.')

        self.nz = len(self.data)
        self.filename = '_'.join(self.filename.split('_')[0:3]) + '_{0}.nc'.format(sensor)  # }}}

    def pick_sensor2(self, sensor, verbose=False):  # {{{
        '''
        Clip the z dimension according to the given ofsets (edges are included).
        pick_sensor2 return a NEW cube variable, pick_sensor modifies the original one.
        :param sensor: name of sensor, str. It should match any of spellings in -sensors_colormap- dictionary

        example of use:
        cube2 = cube.pick_sensor2('SENTINEL-1')
        '''

        from velocity_colormap import sensors_colormap  # manage the different spellings, e.g. "S1" and "Sentinel-1"

        cube2=self.copy_structure()

        if verbose:
            print('Original cube z-size: ', len(self.data))
            print('Sensor={0} layers picking...'.format(sensor))

        SENSORS = sensors_colormap()
        for sl in SENSORS.items():
            if sensor in sl[1]:
                key = sl[0]
                break
        cube2.nz = 0
        for d in self.data:
                if d.sensor in SENSORS[key]:
                    cube2.insert(d.vx, d.vy, d.errx, d.erry,
                                 d.date1, d.date2, d.source,d.sensor,
                                 err2D=not(d.errx.shape==()), verbose=False)
        cube2.nz = len(cube2.data)
        print(cube2.nz, ' layers were picked.')

        cube2.filename = '_'.join(cube2.filename.split('_')[0:3]) + '_{0}.nc'.format(sensor)
        return cube2# }}}

    def sub_cube(self,ii=None, jj=None):
        ''' Extract from the cube a zone [i1:i2, j1:j2, z] INCLUDING the -ii-,-jj- edges in NORMAL enumerate way.
        1st index = 1 (not standard Python array mode).
        [6,10] -> from 6th element to 10th including, 5 items in total

        example of use:
        subc=c.sub_cube([6,10],[25,200]) '''

        # adjust the edges to Python enumeration mode
        ii[0] = ii[0] - 1
        jj[0] = jj[0] - 1
        # check source proj
        if '+proj=stere +lat_0=90 +lat_ts=70' in self.proj4: proj='PS'
        elif '+proj=stere +lat_0=-90 +lat_ts=-71' in self.proj4: proj='PS_SOUTH'
        else: proj='UTM'
        # cube to be filled
        sub = cube_class(nx=ii[1]-ii[0], ny=jj[1]-jj[0], proj=proj)

        # data slices
        vx = self.vx_()[:, jj[0]:jj[1], ii[0]:ii[1]]
        vy = self.vy_()[:, jj[0]:jj[1], ii[0]:ii[1]]
        if self.data[0].errx.shape==(): # a unique-per-layer error
            errx = self.errx_()
            erry = self.errx_()
        else : # 2D error
            errx = self.errx_()[:, jj[0]:jj[1], ii[0]:ii[1]]
            erry = self.erry_()[:, jj[0]:jj[1], ii[0]:ii[1]]

        sub.insert(vx=vx, vy=vy, errx=errx, erry=erry,
                   date1=self.date1_(), date2=self.date2_(), source=self.source_(), sensor=self.sensor_(),
                   err2D=(not self.data[0].errx.shape==()))

        sub.x = self.x[ii[0]:ii[1]]
        sub.y = self.y[jj[0]:jj[1]]
        sub.z = self.z
        sub.i = self.i + ii[0]
        sub.j = self.j + jj[0]
        sub.nz = len(sub.data)
        sub.filedir = self.filedir
        sub.filename = 'subc' + '_x{:05}_y{:05}_'.format(sub.i, sub.j) + '{0}by{1}pix.nc'.format(sub.nx, sub.ny)
        sub.proj4 = self.proj4
        sub.srs = self.srs
        if self.mapping['grid_mapping_name'] == 'transverse_mercator':
            sub.mapping['grid_mapping_name'] = 'transverse_mercator'
            sub.mapping['spatial_ref'] = sub.srs.ExportToWkt()
        sub.geo.xmin = min (sub.x)
        sub.geo.ymax = max (sub.y)
        sub.geo.nrec = sub.ny
        sub.geo.npix = sub.nx
        sub.geo.xposting = self.geo.xposting
        sub.geo.yposting = self.geo.yposting
        sub.geo.posting = self.geo.posting

        # print ('New cube shape nz-ny-nx : ', sub.shape())
        return sub

    def sieve_empty_layers(self, sieve_n=500,verbose=False):
        ''' 
        Delete layers where number of valid pixels is less than -sieve_n- value.
        
        example of use:
        cube.sieve_empty_layers(50)
        
        '''
        if verbose: print('Original cube.nz ', self.nz)

        for i, d in enumerate(self.data):
            if d.vx.count() < sieve_n:
                self.data.pop(i)

        self.nz = len(self.data)
        
        if verbose:
            if self.nz == 0:
                print (' /!\ Empty cube is returned')
            else: 
                print(len(self.data), ' layers were picked.')

# ====== = ====== UTILITIES ====== = ======
    def coord2pix(self, x, y):  # {{{
        ''' Return cube pixel position from cube projected coordinates.
        :param x, y: coordinate in the cube's projection
        example of use:
        ii, jj = cube.coord2pix(-242353, -23145)
        '''

        # manage memory rounding effects (JER: is it useful and true ??)
        x = x + 0.005
        y = y + 0.005

        # convert
        i = int(np.round( (x - self.geo.xmin) / self.geo.xposting))
        j = int(np.round( (y - self.geo.ymax) / self.geo.yposting))
        
        if  i < 0 or i >= self.nx or j < 0 or j >= self.ny : 
            return ValueError, ValueError
        else:
            return i, j
        
    def lonlat2pix(self, lon, lat,verbose=False):  # {{{
        ''' Return the local cube's pixels ID in python index mode (start with 0).
        :param lon, lat: coordinates in the cube's projection in longitude latitude
        example of use:
        ii, jj = cube.coord2pix(-26.5, 78.0)
        '''
        from ReprojectCoords import ReprojectCoords

        srs_ll = osr.SpatialReference()
        srs_ll.ImportFromEPSG(4326)

        if verbose:
            print('cube proj:',self.srs)
            print('latlon proj:',srs_ll)

        coordsx, coordsy = ReprojectCoords( [[lat,lon]] , srs_ll, self.srs)
        if verbose: 
            print(coordsx, coordsy)
        return self.coord2pix(coordsx[0], coordsy[0]) # }}}

    def pix2lonlat(self,i,j,center=False,verbose=False):
        '''
            Convert pixel position in cube into longitude and latitude coordinates

            :param i: pixel position along x-axis 
            :param j: pixel position along y-axis ( cube are (y,x) )
            :param center: True or False return the pixel's center coordinate (True) or left top corner (False)
            :param verbose: True or False

            :return longitude, latitude 
            
        '''
        from ReprojectCoords import ReprojectCoords
        xps, yps = self.pix2coord(i,j,center=center)
        
        srs_ll = osr.SpatialReference()
        srs_ll.ImportFromEPSG(4326)

        if verbose:
            print('cube proj:',self.srs)
            print('latlon proj:',srs_ll)

        lat, lon = ReprojectCoords( [[float(xps),float(yps)]], self.srs, srs_ll)

        return lon, lat

    def pix2coord(self, i, j, center=False):
        ''' Returm the geocoordinates for given cube's pixels ID (in python index mode - start with 0).
        :param i, j: pixel position
        :param center: return the pixel's center coordinate (True) or left top corner (False)

        example of use:
        xx, yy = cube.pix2coord(24,53, True)
        '''

        if not center:
            return self.x[i], self.y[j]
        else:
            return self.x[i] + round(self.geo.xposting / 2, 0), self.y[j] + round(self.geo.yposting / 2, 0)

    def search_cube(self, x, y, greenland=False):
        ''' Search the cube name that contains the given -x-y- point.
            Use only for non-merged spacely cubes! (=> -nx-ny- of a cube == -nx-ny of the grid)
           :param x, y: point geo coordinates, int or float
           :param greenland: the default geo_param is loaded, so the method could be used without any pre-loaded cube (empty -self-)
           :return name_pattern: "x00000_y00000", string

           example of use:
           name = cube.search_cube(-254466,-433245) '''

        if self.filename.split('_')[0]!='c':
            print(' /!\ Use the -search_cube- only for non-merged spacely cubes! (=> -nx-ny- of a cube == -nx-ny- of the grid)')

        if greenland:
            geo = geo_param()
            geo.load_greenland()
            xposting = 150
            yposting = -150
            xmin = geo.xmin
            ymax = geo.ymax
        else:
            xposting = self.geo.xposting
            yposting = self.geo.yposting
            # geo coordinates of global grig corner [0,0]
            xmin = self.geo.xmin - self.geo.xposting * self.i
            ymax = self.geo.ymax - self.geo.yposting * self.j

        # global pixel index
        i = ((x - xmin) / xposting) + 1  # python index mode correction
        j = ((y - ymax) / yposting) + 1
        # cube number
        n = np.ceil((i + (i / self.nx) * 5) / self.nx)
        m = np.ceil((j + (j / self.ny) * 5) / self.ny)
        # global index of cube corner
        gi = (n - 1) * (self.nx - 5)
        gj = (m - 1) * (self.ny - 5)
        name_pattern = 'x' + str(int(gi)).zfill(5) + '_y' + str(int(gj)).zfill(5)
        return name_pattern

    def deepcopy(self): #{{{
        ''' Use it to prevent the ordinar numpy copy behavoir: when "cube2 = cube1" used, a simple link is created,
                                                                so any modification of cube2 affects cube1.

            example of use:
            cube2 = cube1.deepcopy() '''
        import copy

        cube2 = cube_class()
        cube2.filedir = self.filedir
        cube2.filename = self.filename
        cube2.x = self.x
        cube2.y = self.y
        cube2.z = self.z
        cube2.i = self.i
        cube2.j = self.j
        cube2.nx = self.nx
        cube2.ny = self.ny
        cube2.nz = self.nz
        cube2.geo = self.geo
        cube2.srs = self.srs
        cube2.proj4 = self.proj4
        cube2.mapping = self.mapping

        # TODO switch to simpler
        # cube2.data = copy.deepcopy(self.data)
            
        for n in range(self.nz):
            d = copy.deepcopy(self.data[n])
            cube2.data.append(d)

        return cube2

    def toElmer(self, date=None, pick_z=None, saveto=''): #{{{
        ''' Save the given layer as ELMER-Ice compatible NetCDF (without z dimention and with flipped -y- coordinate)
        :param date: the sting type date as '01-01-2020', central layer date is checked (date2-date1/2)
                    If the flex=True, +-1 day will be searched to fix dates calculation variability
        :param z : the layer -z- could be given instead the date
        :param saveto: where to save the new NetCDF, '/dir/dir/name.nc' or '/deir/dir/'. Keep it empty to save in the same folder as the source

        example of use:
        cube.toElmer(date='03-05-2020')
        cube.toElmer(pick_z = 34)
        '''

        if date == None:
            date = mjd2date((self.data[pick_z].date1 + self.data[pick_z].date2) / 2)
        if saveto=='':
            saveto = self.filedir + self.filename.replace('.nc', '_date{0}.nc'.format(date)).replace('c_', 'elmer_')
        elif '.nc' not in saveto: # only folder is given
            saveto = saveto + '/' + self.filename.replace('.nc', '_date{0}.nc'.format(date)).replace('c_', 'elmer_')
        elif '/' not in saveto: # only core name is given
            saveto = self.filedire + saveto

        if pick_z==None: # date is given
            cube2 = self.pick_time2(mjd2date(date2mjd(date) - 1), mjd2date(date2mjd(date) + 1))
            if cube2.nz > 1: cube2.pick_time(date2mjd(date), date2mjd(date))
            cube2.data[0].vx.mask = np.logical_or(cube2.data[0].vx.mask, cube2.data[0].vy.mask)
            cube2.data[0].vy.mask = np.logical_or(cube2.data[0].vx.mask, cube2.data[0].vy.mask)
            cube2.write()
            com = 'ncwa -d z,0 -a z {0} {1} -O'.format(cube2.filedir + cube2.filename, cube2.filedir + 'tech.nc') # delete -z- dimension
            os.system(com)
            com = 'ncks -x -v source,sensor,creation_date {0} {0} -O'.format(cube2.filedir + 'tech.nc') # delete useless vars and dimentions "letter"
            os.system(com)
            os.remove(cube2.filedir + cube2.filename)  # delete tech cube
            if cube2.y[0]>cube2.y[1]:
                com = 'ncpdq -a -y {0} {1} -O'.format(cube2.filedir + 'tech.nc', saveto) # flip -y- coordinate
                os.system(com)
                os.remove(cube2.filedir + 'tech.nc')  # delete tech cube
            else:
                os.replace(cube2.filedir + 'tech.nc' , saveto)
            print('    Convert cube layer to ELMER NetCDF, save to ', saveto)

        else:
            self.data[pick_z].vx.mask = np.logical_or(self.data[pick_z].vx.mask, self.data[pick_z].vy.mask)
            self.data[pick_z].vy.mask = np.logical_or(self.data[pick_z].vx.mask, self.data[pick_z].vy.mask)
            com = 'ncwa -d z,{0} -a z {1} {2} -O'.format(pick_z, self.filedir + self.filename, self.filedir + 'tech.nc') # delete -z- dimension
            os.system(com)
            if self.y[0]>self.y[1]:
                com = 'ncpdq -a -y {0} {1} -O'.format(self.filedir + 'tech.nc', saveto) # flip -y- coordinate
                os.system(com)
                os.remove(self.filedir + 'tech.nc')  # delete tech cube
            else:
                os.replace(self.filedir + 'tech.nc' , saveto)
            print('    Convert cube layer to ELMER NetCDF, save to ', saveto)

    def toElmer_2(self, cube_std, date=None, pick_z=None, saveto=''):
        ''' Save the given layer as ELMER-Ice compatible NetCDF (without z dimention and with flipped -y- coordinate)
            vx + vy + stdx + stdy + num_valid_pix
        :param self: cube with speed, str with path
        :param sube_std: cube with STD in -vx-vy and num_valid_pix in -errx-erry-
        :param date: the sting type date as '01-01-2020', central layer date is checked (date2-date1/2)
                    If the flex=True, +-1 day will be searched to fix dates calculation variability
        :param z : the layer -z- could be given instead the date
        :param saveto: where to save the new NetCDF, '/dir/dir/name.nc' or '/deir/dir/'. Keep it empty to save in the same folder as the source

        example of use:
        cube.toElmer_2('./AvY_SM_2013to2019_jer_RUSpart_LwsSTD.nc', date='03-05-2020')
        cube.toElmer_2('./AvY_SM_2013to2019_jer_RUSpart_LwsSTD.nc',pick_z = 34)
        '''

        def do_it(pick_z, c_v, c_std, saveto):
            tech = os.path.dirname(c_std) + '/tech_std.nc'

            # get only pick_z layer and delete -z- dimension : vv_cube
            com = 'ncwa -d z,{0} -a z {1} {2} -O'.format(pick_z, c_v, os.path.dirname(c_v) + '/tech_vv.nc')
            os.system(com)
            # get only pick_z layer and delete -z- dimension : std_cube
            com = 'ncwa -d z,{0} -a z {1} {2} -O'.format(pick_z, c_std, tech)
            os.system(com)
            # delet erry in std cube_cube
            com = 'ncks -O -x -v error_vy,source,sensor,creation_date {0} {0}'.format(tech)
            os.system(com)
            # rename Var and Attr in std cube_cube
            com = 'ncrename -h -O -v vx,stdx -v vy,stdy -v error_vx,pixn {0}'.format(tech)
            os.system(com)
            # rewrite std attributes
            com = 'ncatted -h -a long_name,stdx,m,string,"STD on vx when 2015-2019 speed averaged" -a standard_name,stdx,m,string,"STD_vx" -a units,stdx,m,string," " {0}'.format(
                tech)
            os.system(com)
            com = 'ncatted -h -a long_name,stdy,m,string,"STD on vy when 2015-2019 speed averaged" -a standard_name,stdy,m,string,"STD_vy" -a units,stdy,m,string," " {0}'.format(
                tech)
            os.system(com)
            com = 'ncatted -h -a long_name,pixn,m,string,"Number of valid pixels when 2015-2019 speed averaged" -a standard_name,pixn,m,string,"valid_pix_N" -a units,pixn,m,string," " {0}'.format(
                tech)
            os.system(com)
            # append vx-vy to the std
            com = 'ncks -v vx,vy -A {0} {1}'.format(os.path.dirname(c_v) + '/tech_vv.nc', tech)
            os.system(com)
            # flip Y
            com = 'ncpdq -a -y {0} {1} -O'.format(tech, saveto)
            os.system(com)

            # delete tech cubes
            os.remove(tech)
            os.remove(os.path.dirname(c_v) + '/tech_vv.nc')




        if date == None:
            date = mjd2date((self.data[pick_z].date1 + self.data[pick_z].date2) / 2)
        if saveto=='':
            saveto = self.filedir + self.filename.replace('.nc', '_date{0}.nc'.format(date)).replace('c_', 'elmer_')
        elif '.nc' not in saveto: # only folder is given
            saveto = saveto + '/' + self.filename.replace('.nc', '_date{0}.nc'.format(date)).replace('c_', 'elmer_')
        elif '/' not in saveto: # only core name is given
            saveto = self.filedire + saveto

        if pick_z==None: # date is given
            cube2 = self.pick_time2(mjd2date(date2mjd(date) - 1), mjd2date(date2mjd(date) + 1))
            if cube2.nz > 1: cube2.pick_time(date2mjd(date), date2mjd(date))
            cube2.data[0].vx.mask = np.logical_or(cube2.data[0].vx.mask, cube2.data[0].vy.mask)
            cube2.data[0].vy.mask = np.logical_or(cube2.data[0].vx.mask, cube2.data[0].vy.mask)
            cube2.write(replace=True)

            do_it(pick_z=0, c_v=cube2.filedir+cube2.filename, c_std=cube_std, saveto=saveto)

            print('    Convert cube layer to ELMER STD_EXTENDED NetCDF, save to ', saveto)

        else:
            self.data[pick_z].vx.mask = np.logical_or(self.data[pick_z].vx.mask, self.data[pick_z].vy.mask)
            self.data[pick_z].vy.mask = np.logical_or(self.data[pick_z].vx.mask, self.data[pick_z].vy.mask)
            self.write(replace=True)

            do_it(pick_z=pick_z, c_v=self.filedir + self.filename, c_std=cube_std, saveto=saveto)


            print('    Convert cube layer to ELMER STD_EXTENDED NetCDF, save to ', saveto)


# ====== = ====== DATA VISUALISATION ====== = ======
    def plot(self,i,j,dmin=datetime(2013,1,1),dmax=datetime(2021,1,1),vmin=0,vmax=5000):

        '''
        i : position x in pixel
        j : position y in pixel

        param dmin : date min value
        param dmax : date max value

        param vmin : v min value
        param vmax : v max value
        '''
        import matplotlib.pyplot as plt

        # get time data
        date1 = self.date1_(datetime=True)
        date2 = self.date2_(datetime=True)
        date= np.asarray(date1) + (np.asarray(date2)-np.asarray(date1))/2.
        # 
        vx = np.nanmedian(self.vx_()[:,j-2:j+2,i-2:i+2],axis=(1,2))
        vy = np.nanmedian(self.vy_()[:,j-2:j+2,i-2:i+2],axis=(1,2))

        v = np.sqrt(vx**2+vy**2)
        errv = np.sqrt(self.errx_()**2+self.erry_()**2)
        colors, labels = self.color_per_sensor_()
            
        plt.scatter(date,v,c=colors,marker='.')

        plt.xlim(dmin,dmax)
        plt.ylim(vmin,vmax)
        plt.ylabel('Ice speed (m/yr)')
        lon,lat = self.pix2lonlat(i,j)
        plt.title('{lon:5.2f}$^\circ$W - {lat:5.2f}$^\circ$N'.format(lon=lon[0],lat=lat[0]))

        #plt.title('{lon:.2f} - {lat:.2f}'.format(lon=lon,lat=lat))

        legend_elements = [ plt.scatter([0], [0], marker='.', color='green', label='S1'), 
                        plt.scatter([0], [0], marker='.', color='purple', label='S2'), 
                        plt.scatter([0], [0], marker='.', color='cornflowerblue', label='L8'),
                        plt.scatter([0], [0], marker='.', color='red', label='ERS'),
                        plt.scatter([0], [0], marker='.', color='yellow', label='RADARSAT'),
                        plt.scatter([0], [0], marker='.', color='orange', label='ALOS'),
                        plt.scatter([0], [0], marker='.', color='pink', label='ENVISAT')]

        plt.legend(handles=legend_elements)

    def plot_vv_space(self, z, show=True, save=False):  # {{{
        import matplotlib.pyplot as plt
        from velocity_colormap import velocity_colormap
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning, message='Warning: converting a masked element to nan.')

        ''' Plot the simple preview for velocity magnitude at given -z- layer as 2D map (with standard log colormap 1-3000).
        :param z: the z index of the layer to plot (count starts from 0), integer
        :param show: open the image window (don't block the code but each plot window should be closed manually at the end)
        :param save: save as png in the source cube folder '''

        # import colormap
        cmap = velocity_colormap()
        # data should be flipped by -y- axe to be north-oriented
        vvmap = np.sqrt(self.data[z].vx ** 2 + self.data[z].vy ** 2)
        if not type(self.z)==str: # not ELMER cube with flipped -y- coordinate
            vvmap = vvmap[::-1]
        else:
            print ('    >>> WARNING: ELMER cube detected. Coordinate -y- is correctly managed at the image but is fliped in the file. <<<')
        # force the 1-to-3000m/y data range for unified color visualisation
        vvmap = np.clip(vvmap, 1, 3000)
        vvmap[0, 0] = 1
        vvmap[0, 1] = 3000
        # Scale the data by log and by 0-to-255 for visualisation with defined log colormap
        vvmap = np.log(vvmap)
        vvmap = (255 * (vvmap - np.nanmin(vvmap)) / (np.nanmax(vvmap) - np.nanmin(vvmap)))
        # PLOT DATA : solution with a workaround
        plt.figure(figsize=(7, 7))
        fig = plt.imshow(np.full_like(vvmap.data, np.nan), cmap=cmap, origin='lower') # make the empty dataframe with square pixels space
        ax = plt.gca()
        ax.pcolormesh(vvmap, cmap=cmap) # apply the plot function with correct interpolation management
        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)

        plt.title('{file},layer {i}, {date1} : {date2}'.format(file=self.filename, i=z,
                                                               date1=mjd2date(self.data[z].date1),
                                                               date2=mjd2date(self.data[z].date2)), fontsize=10)
        plt.tight_layout()

        if show:
            plt.show(block=False)

        if save:
            savepath = os.path.join(self.filedir, self.filename.split('.nc')[0]) + '_layer' + str(z) + '.png'
            if os.path.exists(savepath): savepath.replace('.png', '_1.png')
            print('Saving figure to {0}...'.format(savepath))
            plt.savefig(savepath, dpi=150, bbox_inches='tight', pad_inches=0.1)
            plt.close()  # }}}

    def plot_vv_time(self, x, y, geo_or_ind='geo', show=True, save=False,
                     subplt=False, show_offset=False, filter_holes=False, show_error=False,
                     show_line=True, color=None, show_label=True,
                     window=5):         #  {{{
        import matplotlib.pyplot as plt
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning, message='Warning: converting a masked element to nan.')

        '''
        Plot the velocity magnitude at the given point as time profile.
        :param geo_or_ind: geo coordinates -x-y- or internal cube pixel indexes -i-j- are given, string 'geo' / 'ind'
        :param show: open the image window (don't block the code but each plot window should be closed manually)
        :param save: save as png in the source cube folder
        :param subplt: add data to the existing current plot
        :param show_offset: plot the ticks of date1-date2 interval (not nice to use with the long x-ax time range, ~ more than 1 year)
        :param filter_holes: filter the dates without velocity data to have a continuous plotline
        :param show_error: plot velocity error bar
        :param show_line: drow the line connected velocity points (according to the filer_holes impact)
        :param show_label: put (True) or hide(False) the item in the legend, or write some specific text (str)
        :param window: the averaging (be median) window, 5*5 pix by default
        '''

        if geo_or_ind == 'ind':
            i = x
            j = y #j = self.ny - y - 1
        else:
            # transform geo coordinates to the pixel indexes
            i = self.coord2pix(x, y)[0]
            j = self.coord2pix(x, y)[1]
            if i==ValueError:
                raise IndexError ("Coordinates out of the cube extent. Use 'ind' to use pixel-based values")

        w = round(window / 2)
        # get time data
        date1 = np.array([mjd2date(d) for d in self.date1_()])
        date2 = np.array([mjd2date(d) for d in self.date2_()])
        offset_bar = (date2 - date1) / 2
        date_c = date2 - offset_bar
        # get velocity data
        vx = np.nanmedian(np.nanmedian(self.vx_()[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
        vy = np.nanmedian(np.nanmedian(self.vy_()[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
        vv = np.round(np.sqrt(vx ** 2 + vy ** 2), 2)
        vv = np.ma.masked_invalid(vv)
        if np.isnan(np.nanmin(vv.data)) or (not np.any(vv.data)):  # edge indexes, force window=1 mode
            vx = self.vx_()[:, j, i]
            vy = self.vy_()[:, j, i]
            vv = np.round(np.sqrt(vx ** 2 + vy ** 2), 2)
        if self.data[0].errx.shape != ():  # 2D error
            errx = np.nanmedian(np.nanmedian(self.errx_()[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
            erry = np.nanmedian(np.nanmedian(self.errx_()[:, j - w: j + w + 1, i - w: i + w + 1], 1), 1)
            error = np.round((errx * abs(vx) + erry * abs(vy)) / vv, 2)
                    ### old mode: np.round(np.sqrt(errx ** 2 + erry ** 2) / 2, 2)
        else:  # 1D error
            error = np.round((self.errx_() * abs(vx) + self.erry_() * abs(vy)) / vv, 2)

        if not vv.any():
            print('NO VALID VELOCITY DATA at x:y {0}:{1} // cube {2}'.format(x, y, self.filename))
        else:
            # combine data to sort them by date
            data = np.array([[date_c[n], offset_bar[n], vv[n], error[n]] for n in range(self.nz)])
            data = np.array(sorted(data, key=lambda date: date[0]))

            if filter_holes:
                data_v = []
                for d in data:
                    if d.all(): data_v.append(d)
                data = np.array(data_v)

            if show_label == False:
                label = None
            elif type(show_label) == str:
                label = show_label
            else:
                label = '[{0}, {1}]'.format(x, y)

            if show_line == False or show_line == None:
                linestyle = ''
            else:
                linestyle = '-'

            # create figure or get current ax
            if subplt:
                ax = plt.gca()
                ymin = vv.min() - error[np.where(vv == np.min(vv))][0] - 20
                ymax = vv.max() + error[np.where(vv == np.max(vv))][0] + 20
                ax.set_ylim(min(max(ymin, 0),ax.get_ylim()[0]), ymax)
            else:
                fig, ax = plt.subplots(1, 1, figsize=(12, 4))
                ymin = vv.min() - error[np.where(vv == np.min(vv))][0] - 20
                ymax = vv.max() + error[np.where(vv == np.max(vv))][0] + 20
                ax.set_ylim(max(ymin, 0), ymax)
                ax.set_ylabel('Velocity magnitude, m/y')
                plt.grid(True, 'major', ls='--', lw=0.5, c='k', alpha=0.3)
            # plot lines, points, offset, error bars
            p = ax.plot(data[:, 0], data[:, 2], linestyle=linestyle, label=label,
                        zorder=1, marker='o', lw=0.7, markersize=2, color=color)
            if show_offset:
                color = p[0].get_color()
                ax.errorbar(data[:, 0], data[:, 2], xerr=data[:, 1], alpha=0.3, fmt=',', c=color, zorder=1)
            if show_error:
                color = p[0].get_color()
                ax.errorbar(data[:, 0], data[:, 2], yerr=data[:, 3], alpha=0.3, fmt=',', c=color, zorder=1)

            plt.legend(loc='best')

            if show:
                plt.show(block=False)

            if save:
                saveto = os.path.join(self.filedir, self.filename.split('.nc')[0]) + '_x{0}y{1}.png'.format(str(x), str(y))
                if os.path.exists(saveto): saveto = saveto.replace('.png', '_1.png')
                print('Saving figure to {0}...'.format(saveto))
                plt.savefig(saveto, dpi=150, bbox_inches='tight', pad_inches=0.1)
                plt.close()

# ====== = ====== DATA PROCESSING ====== = ======
    # functions' cores are replaced to /ST_RELEASE/MOSAIC/PYTHON/cube_processing.py

    def filter_space(self, space_window=3):
        '''Median space filter with moving window of given size. Velocity magnitude is filtered and vx-vy are restored.
            example of use:
            cube.filter_space() '''
        from src.processing.archive.cube_processing import filter_space
        filter_space(self, space_window)

    def filter_time_stat(self, vvar, threshold=1, replace=True, extended_out=False): #{{{
        '''
        Simple statistical filter: throw out the vx/vy there abs(time_median - value) > max(pixel_std, stdmax)
        :param vvar: filtering criteria, fraction of the pixel's mean velocity that is acceptable
                    as +- variation from the v[z] value if it is > than STD*threshold (e.g. 1STD, 3STD).
                       Not more than 0.5 could be recomended.
        :param threshold: STD*-threshold- (e.g. 1STD, 3STD)
        :param replace: True = replace the value in the source cube, False = np.array outputs
        :param extended_out: True = include statistics matrix to outputs
        :return:
                overwrited cube with filtered velocity ( -> replace = True)
                vxf, vyf ( -> replace = False, -> extended_out= False)
                median and STD maps: medvx, medvy, stdvx, stdvy ( -> extended_out = True)

        example of use:
        cube.filter_time_stat(0.4)
        vxf, vyf, medvx, medvy, stdvx, stdvy =  cube.filter_time_stat(0.4, replace = False, extended_out = True)
        '''

        from src.processing.archive.cube_processing import filter_time_stat
        if replace:
            filter_time_stat(self, vvar, threshold, replace, extended_out)
        else:
            if extended_out:
                vxf, vyf, medvx, medvy, stdvx, stdvy = filter_time_stat(self, vvar, threshold, replace, extended_out)
                return vxf, vyf, medvx, medvy, stdvx, stdvy
            else:
                vxf, vyf = filter_time_stat(self, vvar, threshold, replace, extended_out)
                return vxf, vyf  # }}}

    def filter_time_RollM(self, RollWindow=15, RollWindowA=60, RollMinN=3, window=1, verbose=False):
        ''' Apply Rolling MEDIAN filter to velocity magnitude & direction through the pixels (time-dimension), restore vx-vy.
        Output error is local rolling STD of the magnitude.
        :param RollWindowV: the rolling window size in days for the velocity magnitude, int
        :param RollWindowA: the rolling window size in days for the direction (30 days min recommended), int
        :param RollMinN: the minimum quantity of valid values in a window to do calculations, int
        :param window: a space window to do calculation (1 = one pixel, 3 = mean in a 3*3 space box, etc), int

        example of use:
        cube.filter_time_RollM()
        '''
        from src.processing.archive.cube_processing import filter_time_RollM
        filter_time_RollM(self, RollWindow, RollWindowA, RollMinN, window, verbose)

    def filter_time_lowess(self, lowess_window=20, save2cube='', drop_duplicates=True, mode='replace_orig', verbose=False): #{{{
        '''Apply LOWESS filter to Vx-Vyю
            Does NOT re-estimate errors, keeps the source-data values
        :param lowess_window: proxi of number of surrounding points that are taken into accounting for the regression line calculation
        :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
        :param drop_duplicates: before do LOWESS, keep only 1 value per day if many layers are presented (median is taken)
        :param mode: only one mode is implemented currently
                     'replace_orig' -> replace pixels value with lowess estimation, only for non-masked source pixels
        example of use:
        cube_lws = cube.filter_time_lowess(save2cube='lowess with 20pnt')
        vx_f, vy_f  =  cube.filter_time_lowess()
        '''

        from src.processing.archive.cube_processing import filter_time_lowess
        if save2cube != '':
            c_lws = filter_time_lowess(self, lowess_window, save2cube, drop_duplicates, mode, verbose)
            return c_lws
        else:
            print(' RETURN: 2 NUMPY ARRAYS')
            vx_f, vy_f = filter_time_lowess(self, lowess_window, save2cube, drop_duplicates, mode, verbose)
            return vx_f, vy_f
            #}}}

    def filter_time_spline(self, dates_out, smooth=0.05, save2cube='', verbose=False, i=None, j=None):
        ''' Apply the CUBIC spline fitting
            Does NOT re-estimate errors, fill it with 0
            cube_out.date1 = cube_out.date2 = dates_out
            :param dates_out : date_out as Modif Julian date (np.array) or pandas interval text ('D', 'W', 'SM', etc)
            :param smooth: between 0 and 1 (values between 0.001 and 0.0001 seem to work when weghts are not set),
                            try 0.05 as reference
            :param save2cube: if you want a cube variable as the output, fill it with the str description for the -source- metadata
            :param i, j: if you want apply the fitting only at one-pixel, xy-yf np.array as output

            need : pip install csaps

            example of use:
            c_spl = c.filter_time_spline('W', 0.05, save2cube='spline 0.05', verbose=True)
            c_spl = c.filter_time_spline(dates_out=[111, 222, 333], 0.05, save2cube='spline 0.05', verbose=True)
            vx_f, vy_f = c.filter_time_spline('W', 0.05, save2cube='', verbose=True)
        '''

        from src.processing.archive.cube_processing import filter_time_spline
        if save2cube != '':
            c_spl = filter_time_spline(self, dates_out, smooth, save2cube, verbose, i, j)
            return c_spl
        else:
            vx_f, vy_f = filter_time_spline(self, dates_out, smooth, save2cube, verbose, i, j)
            return vx_f, vy_f

    def cube2vel(self, plot_preview=False): #{{{
        ''':return:
               status
               filtered weighted mean maps: vxo, vyo
               filtered weighted standard deviation maps: stdxo stdyo
               filtered count map: nno
               median maps: medxo medyo
               standard deviation maps: stdxo0 stdyo0
               err maps : errxo erryo
               count map : nno0 '''
        from src.processing.archive.cube_processing import cube2vel
        status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = \
            cube2vel(self, plot_preview)
        return status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 #}}}

    def cube2vel_py(self, time_filter=False, plot_preview=False, save2cube='', verbose=False): #{{{
        '''
        Calculate weighted mean velocity and some statistics.
        :param space_filter: True -> apply space median filter with -space_window- mooving window
        :param time_filter: apply time filter filter_time_stat with -stat_stdmax- threshold (similar to Fortran)
        :param plot_preview: plot the mean velocity, weighted std and error, layers N, bool
        :param save2cube: place a description for cube -source- attribute to save the result in cube
        :return:
               status
               filtered weighted mean maps: vxo, vyo
               filtered weighted standard deviation maps: stdxo stdyo
               filtered count map: nno
               median maps: medxo medyo
               standard deviation maps: stdxo0 stdyo0
               err maps : errxo erryo
               count map : nno0
        '''
        from src.processing.archive.cube_processing import cube2vel_py
        if save2cube != '':
            c2 = cube2vel_py(self, time_filter, plot_preview, save2cube,  verbose)
            return c2
        else:
            status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 = \
                cube2vel_py(self, time_filter, plot_preview, save2cube,  verbose)
            return status, vxo, vyo, stdxo, stdyo, nno, medxo, medyo, stdxo0, stdyo0, errxo, erryo, nno0 #}}}

    def time_regularisation(self, time_interval, time_step, func='mean', weighted = False, verbose=False):
        ''' Regularize z time step (date1-date2). When doing, each velocity layer is attributed to the central date of the date1-date2 interval.
        Does NOT fill data gaps. Generates the empty layers for time intervals without data.
        :param cube: cube_class variable
        :param time_interval: start-end dates, ['01-01-2015', '31-12-2019']
        :param time_step: 'D': day, 'W': week, 'SM': half month, 'M': month, 'Y': year
                        use the number+letter to derive specific interval, such as '6M'=half of year
        :param func: function of data aggregation, 'mean' or 'median' or 'std'
        :param weighted: False = simple numpy calculations, True = weighted calculations with cube2vel
        :return: cube_class variables

        example of use:
        c_week = c.time_regularisation(['01-01-2015', '31-12-2019'], 'SM', func='mean')
        '''

        from src.processing.archive.cube_processing import time_regularisation
        c_m = time_regularisation(self, time_interval, time_step, func, weighted, verbose) #}}}
        return c_m

    def average_year(self, time_step, time_interval=-1, func='median', weighted=False, verbose=False):
        ''' Calculate the typical average velocity by given -time_steps- of an entire year
            using the source layers from the given -time_interval-.
            For "mean" and "median" modes -error- = STD/sqrt(number_of_valid_pixels)
            For "std" mode -error- = number_of_valid_pixels
            Does NOT fill data gaps or empty time intervals.

            (*) use it with "std" mode to get the second cube required in cube.toElmer_2

        :param cube: cube_class variable
        :param time_step: 'D': day, 'W': week, 'SM': 2 weeks, 'M': month, 'Y': year
                          use the number+letter to derive your interval, such as '6M'=half of year
        :param time_interval: start-end dates in dd-mm format or "-1" for the entire year, ['21-04', '31-12']
        :param func: function of data aggregation, 'mean' or 'median' or 'std'
        :param weighted: False = simple numpy calculations, True = weighted calculations with cube2vel
        :return: cube_class variables

        example of use:
        c_AvY = c.average_year('SM', time_interval = ['01-01', '31-12'])
        c_AvY = c.average_year('SM', func='std', verbose=True)
        '''

        from src.processing.archive.cube_processing import average_year
        c_a = average_year(self, time_step, time_interval, func, weighted, verbose) #}}}
        return c_a


# ========= = ========== = ========== = ======== = ========== = ==========
class cube_element:  # class cube_element {{{

    def __init__(self, nx=250, ny=250, err2D=False):
        self.vx = np.ma.array(np.full([ny, nx], np.nan, dtype=np.float32),fill_value = 0)
        self.vy = np.ma.array(np.full([ny, nx], np.nan, dtype=np.float32), fill_value = 0)
        if not err2D:
            self.errx = np.float32(0.)
            self.erry = np.float32(0.)
        else:
            self.errx = np.ma.array(np.full([ny, nx], 0, dtype=np.float32),fill_value = 0)
            self.erry = np.ma.array(np.full([ny, nx], 0, dtype=np.float32),fill_value = 0)
        self.date1 = 0
        self.date2 = 0
        self.offset = 0
        self.source = ''  # UNIQ IDENTIFIER
        self.sensor = ''
        # self.creation_date = 0.   # }}}

