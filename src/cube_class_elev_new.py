import numpy as np
from netCDF4 import Dataset
import netCDF4
from osgeo import osr
import os
import time
from datetime import datetime, timedelta
from fparam import geo_param

from mjd2date import mjd2date, date2mjd # /ST_RELEASE/UTILITIES/PYTHON/mjd2date.py
from progress_bar import progress_bar # # /ST_RELEASE/UTILITIES/PYTHON/progress_bar.py
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message='Warning: converting a masked element to nan.')

class cube_class:

    def __init__(self, nx=250, ny=250, proj='PS'):  # {{{
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
        self.proj4 = ''
        self.srs = osr.SpatialReference()
        self.data = []
        if proj == 'PS' or proj == 'PS_NORTH':
            self.geo = geo_param()
            self.srs.ImportFromEPSG(3413)
            self.proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
            self.mapping = {'ellipsoid': 'WGS84', 'false_easting': 0.0, 'false_northing': 0.0, \
                            'grid_mapping_name': 'polar_stereographic', 'latitude_of_projection_origin': 90.0, \
                            'standard_parallel': 70.0, 'straight_vertical_longitude_from_pole': -45.0}
        if proj == 'PS_SOUTH':
            self.geo = geo_param()
            self.srs.ImportFromEPSG(3031)
            self.proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
            self.mapping = {'ellipsoid': 'WGS84', 'false_easting': 0.0, 'false_northing': 0.0, \
                    'grid_mapping_name': 'polar_stereographic', 'latitude_of_projection_origin': -90.0, \
                    'standard_parallel': -71.0, 'straight_vertical_longitude_from_pole': 0.0}
        if proj == 'UTM':
            print('Cube class UTM')
            self.geo = geo_param(utm=True)
            self.mapping = { \
                'grid_mapping_name': 'transverse_mercator', 'spatial_ref': ''}
        self.geo.npix=int(self.geo.npix)
        self.geo.nrec = int(self.geo.nrec)

# ====== = ====== FILE INPUT/OUTPUT BASICS  ====== = ======
    # }}}


    def write(self, replace=False,compress=False):  # {{{
        ''' Save the cube variable to the NetCDF file. If the filename.nc exist and -replace- is False, filename_1.nc is created.
        example of use:
        cube.write()
        '''

        if type(self.z) == str: # ELMER cube, -z- dimension does not exist
            raise KeyError("    >>>>>> ELMER cube detected (without -z- dimension). You can't use cube.write() <<<<<<< ")

        nz = self.nz
        ny = self.ny
        nx = self.nx

        print('    Save cube ', nz, ny, nx)
        filepath = os.path.join(self.filedir, self.filename)
        print(filepath)
        #if replace: os.system('rm ' + filepath) # force re-writing

        #else:
        #   if os.path.exists(filepath): # save filename_1.nc
        #       self.filename = self.filename.replace('.nc', '_1.nc')
        #       filepath = os.path.join(self.filedir, self.filename)
        #       if os.path.exists(filepath): # remove old *_1.nc
        #3           os.system('rm ' + filepath)

        #if not os.path.exists(filepath):
        dataset = Dataset(os.path.join(self.filedir, self.filename), 'w', format='NETCDF4')

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
        h = dataset.createVariable('h', 'f4', ('z', 'y', 'x'), fill_value=0, zlib=compress, least_significant_digit=2,complevel=1)
        errh = dataset.createVariable('errh', 'f4', ('z', 'y', 'x'), fill_value=0, zlib=compress, least_significant_digit=2,complevel=1)
        m = dataset.createVariable('m', 'f4', ('y', 'x'))

        date1 = dataset.createVariable('date1', 'i4', 'z')
        date2 = dataset.createVariable('date2', 'i4', 'z')
        source = dataset.createVariable('source', 'S1', ('z', 'nletter'), zlib=compress,complevel=1)
        sensor = dataset.createVariable('sensor', 'S1', ('z', 'nletter2'), zlib=compress,complevel=1)
        creation_date = dataset.createVariable('creation_date', 'S1', ('z', 'nletter'), zlib=compress)

        dataset.Projection = self.geo.projection
        dataset.proj4 = self.proj4
        print(dataset)

        if self.mapping['grid_mapping_name'] == 'polar_stereographic':
            mapping.grid_mapping_name = 'polar_stereographic'
            mapping.setncattr_string('ellipsoid', self.mapping['ellipsoid'])
            mapping.false_easting = self.mapping['false_easting']
            mapping.false_northing = self.mapping['false_northing']
            mapping.setncattr_string('grid_mapping_name', self.mapping['grid_mapping_name'])
            mapping.latitude_of_projection_origin = self.mapping['latitude_of_projection_origin']
            mapping.standard_parallel = self.mapping['standard_parallel']
            mapping.straight_vertical_longitude_from_pole = self.mapping['straight_vertical_longitude_from_pole']

        if self.mapping['grid_mapping_name'] == 'transverse_mercator':
            mapping.grid_mapping_name = 'transverse_mercator'
            mapping.spatial_ref = self.srs.ExportToWkt()

        h.setncattr_string('long_name', 'Surface elevation')
        h.setncattr_string('standard_name', 'elevation')
        h.units = 'meter'
        h.setncattr_string('grid_mapping', 'mapping')
        h.missing_value = 0

        errh.setncattr_string('long_name', 'Error Melt rate')
        errh.setncattr_string('standard_name', 'elevation')
        errh.units = 'meter'
        errh.setncattr_string('grid_mapping', 'mapping')
        errh.missing_value = 0
        
        m.setncattr_string('long_name', 'Ice mask')
        m.setncattr_string('standard_name', 'mask')
        m.units = 'n/a'
        m.setncattr_string('grid_mapping', 'mapping')
        m.missing_value = 0
        
        x.setncattr_string('long_name', 'Cartesian x-coordinate')
        x.setncattr_string('standard_name', 'projection_x_coordinate')
        x.units = 'meter'
        x.axis = "X"

        y.setncattr_string('long_name', 'Cartesian y-coordinate')
        y.setncattr_string('standard_name', 'projection_y_coordinate')
        y.units = 'meter'
        y.axis = "Y"

        date1.setncattr_string('long_name', 'Time of acquisition')
        date1.setncattr_string('standard_name', 'time_acquisition') 
        date1.content = 'Date of acquisition'
        date1.units = 'day (Modified Julian Date)'
        date1.axis = "Z"

        date2.setncattr_string('long_name', 'Time of acquisition')
        date2.setncattr_string('standard_name', 'time_acquisition') 
        date2.content = 'Date of acquisition'
        date2.units = 'day (Modified Julian Date)'
        date2.axis = "Z"

        creation_date.content = 'Date of last file modif. (update purposes)'

        # add data into the NetCDF
        x[:] = np.float32(self.x)
        y[:] = np.float32(self.y)
        h[:,:,:] = self.h_()
        errh[:,:,:] = self.errh_()
        m[:,:] = self.m

        #try:
        for z in range(self.nz):
            progress_bar( z/self.nz*100. , msg = str(z/self.nz*100.)+'% layers saved in nc')
            #vx[z, :, :] = self.data[z].vx # fill_value in masked_array = 0, in nc_variable = 0
            #vy[z, :, :] = self.data[z].vy
            #h[z, :, :] = self.data[z].h

            date1[z] = self.data[z].date1
            date2[z] = self.data[z].date2

            tmp = netCDF4.stringtoarr(self.data[z].source, len(self.data[z].source))
            source[z, 0:len(tmp)] = netCDF4.stringtoarr(self.data[z].source, len(self.data[z].source))
            tmp = netCDF4.stringtoarr(self.data[z].sensor, len(self.data[z].sensor))
            sensor[z, 0:30] = netCDF4.stringtoarr(self.data[z].sensor[0:30], 30)
            # tmp = netCDF4.stringtoarr(self.data[i].creation_date, len(self.data[i].creation_date))
            # creation_date[i,0:len(tmp)] = netCDF4.stringtoarr(self.data[i].creation_date, len(self.data[i].creation_date))

            #dataset.close()
            print('     {0} data layers successfully saved in {1}/{2}'.format(len(self.data), self.filedir, self.filename))
            
        #except:
        #    print(dataset)
        #    dataset.close()
        #    os.remove(filepath)
        #    print ('    >>>>>>>> WARNING! ', filepath)
        #    print ('    >>>>>>>> WARNING! Saving was crashed, some data does not fit NetCDF structure or datatypes. ')
    # }}}
    def describe(self):  # {{{

        print('Filename :', self.filename)
        print('i, j (pixel position in full grid) :', self.i, self.j)
        print('Projection (proj4) : ', self.proj4)
        print('geo.xmin :', self.geo.xmin)
        print('geo.ymax :', self.geo.ymax)
        print('geo.posting :', self.geo.posting)
        print('geo.xposting :', self.geo.xposting)
        print('geo.yposting :', self.geo.yposting)
        print('geo.npix (x) :', self.geo.npix)
        print('geo.nrec (y) :', self.geo.nrec)
        print('Number of layers in data : ', self.nz)
        print('Overlap between adjacent cube is always 5 pixels.')
    # }}}
    def shape(self): # {{{
        return np.array([self.nz, self.ny, self.nx])
    # }}}
    def merge(self, cube2):  # {{{
        '''
        Merge spatially two cubes with adjasent space position.
        :param cube2: path to the NetCDF file (str) or cube_class variable

        example of use:
        cube = cube.merge('path/path/filename2.nc')
        cube = cube.merge(cube2)
        '''

        if type(cube2)==str:
            cube2_path = cube2
            cube2 = cube_class()
            cube2.load(cube2_path)

        nx = self.nx
        ny = self.ny
        nx0 = self.nx
        ny0 = self.ny
        x = self.x
        y = self.y
        z = self.z
        # check if adjacent cube:
        if cube2.y[0] >= self.y[ny - 5] or cube2.y[cube2.ny - 5] <= self.y[0]:
            if cube2.x[0] <= self.x[nx - 5] or cube2.x[cube2.nx - 5] >= self.x[0]:
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
                #big.i = cube2.i if cube2.x[0] < self.x[0] else self.i
                if cube2.x[0] < self.x[0]:
                    big.i = cube2.i
                else:
                    big.i= self.i

                if cube2.y[0]> self.y[0] or cube2.y[0] < self.y[0]:
                    big.j = cube2.j
                else:
                    print('HERE')
                    big.j = self.j
                big.j=self.j
                #big.j = self.j
                #big.j = cube2.j if cube2.y[0] > self.y[0] else self.j  ## y pixel are inverse !
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
                print(nx,ny)
                for cube_el in self.data:
                    big_cube_element = cube_element(nx=nx, ny=ny)
                    #big_cube_element.h[cube2.j - big.j:cube2.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.h
                    big_cube_element.h[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.h
                    big_cube_element.errh[self.j - big.j:self.j - big.j + ny0,self.i - big.i:self.i - big.i + nx0] = cube_el.errh
                    big_cube_element.source = cube_el.source
                    big_cube_element.sensor = cube_el.sensor
                    big_cube_element.date1 = cube_el.date1
                    big.data.append(big_cube_element)

                for ind in range(0,len(cube2.data)): # in cube2.data:
                    # the same source = the pieces of the same image
                    cube_el=cube2.data[ind]

                    if (cube_el.date1 in big.date1_()):
                        z = big.date1_().index(cube_el.date1)
                        #big.data[z].h[self.j - big.j : self.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.h
                        big.data[z].h[cube2.j - big.j : cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.h
                        big.data[z].errh[cube2.j - big.j : cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.errh
                    
                    else: #1D error kept from the self-cube
                        #big_cube_element.h[self.j - big.j:self.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.h
                        big_cube_element.h[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.h
                        big_cube_element.errh[cube2.j - big.j:cube2.j - big.j + 250,cube2.i - big.i:cube2.i - big.i + 250] = cube_el.errh
                        big_cube_element.source = cube_el.source
                        big_cube_element.sensor = cube_el.sensor
                        big_cube_element.date1 = cube_el.date1
                        big_cube_element.offset = cube_el.date1
                        big.data.append(big_cube_element)

                big.nz = len(big.data)
                print('Merged successfully. Current cube.filename ', big.filename)
                print('Current big_cube size nz-ny-nx: ', big.nz, big.ny, big.nx)
                self = big
                del big

            else: print('   Cubes are not adjacent. No successful merge could be done.')
        else: print('   Cubes are not adjacent. No successful merge could be done.')

        return self

    def insert(self, h, date1, source, sensor, err2D=False, replace=False, verbose=True):
        '''Add external python-like data to the cube. Data as one layer pack (z=1) or multy-layers list/array
         :param vx, vy, h: velocity -z-y-x-,  numpy array
         :param errx, erry: error as one value per layer -z- or 2D map -z-x-y- (set err2D=True), numpy array
         :param date1, date2: dates in MJD, integer or -z- array
         :param source, sensor: text descriptions, str or -z- list
         :param err2D : check True in the error representation is 2D map
         :param replace : check it to rewrite cube from zero, otherwise the data will be appended
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

                    data.h = np.ma.masked_invalid(h[z,:,:])
                    data.errh = np.ma.masked_invalid(errh[z,:,:])

                except: raise ValueError(' ---> Velocity or errors are in unacceptable format (array case).')
                try:
                    data.date1 = date1[z]
                    data.offset = date1[z]
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
                data.errh = np.ma.masked_invalid(h[:, :])
                data.errh = np.ma.masked_invalid(errh[:, :])

            except:raise ValueError(' ---> Velocity or errors are in unacceptable format (array case).')
            try:
                data.date1 = date1
                data.offset = date1
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
           # }}}

    def load(self, filepath, concat=False, pick_month=False, sensor=None, verbose=False):  # {{{
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
            self.__init__()

            self.filedir = os.path.dirname(filepath) + '/'
            self.filename = os.path.basename(filepath)
            self.proj4 = nc.proj4
            self.x = nc.variables['x'][:]
            self.y = nc.variables['y'][:]
            self.m= nc.variables['m'][:,:]
            self.nx = len(self.x)
            self.ny = len(self.y)
            self.i = nc.i
            self.j = nc.j

            self.geo.xmin = min(self.x)
            self.geo.ymax = max(self.y)
            self.geo.npix = len(self.x)
            self.geo.nrec = len(self.y)
            self.geo.posting = self.x[1] - self.x[0]
            self.geo.xposting = self.x[1] - self.x[0]
            self.geo.yposting = self.y[1] - self.y[0]
            self.geo.projection = nc.Projection

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

        try:
            # nc.dimensions['z'] # if "normal" cube with -z- dimension, does not generate an error. load it.
            h = nc.variables['h'][:,:,:]
            errh = nc.variables['errh'][:,:,:]
            date1 = list(nc.variables['date1'][:])
            date2 = list(nc.variables['date2'][:])

            for z in range(nc.nz):
                data = cube_element(nx=self.geo.npix, ny=self.geo.nrec)
                data.h = h[z,:,:]  # masked array type by default, NetCDF -missing_value- attribute is masked.
                data.errh = errh[z,:,:]  # masked array type by default, NetCDF -missing_value- attribute is masked.
                #data.h = data.h.mask | (data.h == 0) # add the compatibility with source-cubes compiled from tifs
                data.date1 = date1[z]
                data.date2 = date2[z]
                data.offset = data.date1
                data.source = str(nc.variables['source'][z].data, 'utf-8')
                if 'sensor' in nc.variables:
                    data.sensor = str(nc.variables['sensor'][z].data, 'utf-8')
                # data.creation_date = nc.variables['creation_date'][z]
                data.source = data.source.strip('\x00')
                if sensor is not None:
                    data.sensor = sensor
                else:
                    data.sensor = data.sensor.strip('\x00')
                self.data.append(data)
                self.nz = self.nz + 1

        except: # ELMER compatible cube (without z dimension, 1 layer)
            print('     >>> WARNING: ELMER cube detected (without -z- dimensions) <<< ')
            print('     >>> WARNING: coordinate -y- is fliped here! <<< ')

            data.h = nc.variables['h'][:]  # masked array type by default, NetCDF -missing_value- attribute is masked.
            #data.h.mask = data.h.mask | (data.h == 0)  # add the compatibility with source-cubes compiled from tifs
            data.date1 = nc.variables['date1'][:]
            data.offset = data.date1
            data.source = str(nc.variables['source'][:].data, 'utf-8')
            if 'sensor' in nc.variables:
                data.sensor = str(nc.variables['sensor'][:].data, 'utf-8')
            data.source = data.source.strip('\x00')
            data.sensor = data.sensor.strip('\x00')
            self.data.append(data)
            self.nz = 1

        self.shape = (self.nz, self.ny, self.nx)
        nc.close() # }}}

    def update(self, c2filepath, update=True, verbose=True):
        ''' Update cube1's layers with layers from -c2filepath- :
        if update=True:
            if the same source detected, values in the cube1 are replaced with cube2 values;
            if the cube2's layer not yet in the cube1, append it
        if update=False:
            only append the new layers, don't rewrite the existing data

        syntax example:
        cube.update('/dir/dir/name.nc', False, False)   '''
        # {{{
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
                    self.data[z1].h = c2.data[z2].h
                    self.data[z1].errh = c2.data[z2].errh
                    self.data[z1].source = c2.data[z2].source
                else: # append new layer
                    if verbose: print(' append: ', c2.data[z2].source)
                    self.insert(c2.data[z2].h, c2.data[z2].errh,
                                c2.data[z2].date1,c2.data[z2].date2,  c2.data[z2].source, c2.data[z2].sensor) # }}}

    def copy_structure(self):

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


# ====== = ====== CONVERT CUBES DATA TO LIST OR ARRAY ====== = ======
    def source_(self):
        return [d.source for d in self.data]

    def sensor_(self):
        return [d.sensor for d in self.data]

    def h_(self, i=None, j=None):
        ''' Get -h- data for the whole cube (i=None, j=None), for a z-slice (for given i or j), or for a pixel -i-j- '''
        # Dimensions are z, y, x

        if j is not None:
            if j > self.ny: j = None
        if i is not None:
            if i > self.nx: i = None

        if i is None and j is None:
            h = np.ma.zeros((self.nz, self.ny, self.nx))
            for z, d in enumerate(self.data): h[z, :, :] = d.h
        elif i is None and j is not None:
            h = np.ma.zeros((self.nz, 1, self.nx))
            for z, d in enumerate(self.data): h[z, 0, :] = d.h[j, :]
        elif i is not None and j is None:
            h = np.ma.zeros((self.nz, self.ny, 1))
            for z, d in enumerate(self.data): h[z, :, 0] = d.h[:, i]
        elif i is not None and j is not None:
            h = np.ma.zeros((self.nz))
            for z, d in enumerate(self.data): h[z] = d.h[j, i]

        return h

    def errh_(self, i=None, j=None):
        ''' Get -h- data for the whole cube (i=None, j=None), for a z-slice (for given i or j), or for a pixel -i-j- '''
        # Dimensions are z, y, x

        if j is not None:
            if j > self.ny: j = None
        if i is not None:
            if i > self.nx: i = None

        if i is None and j is None:
            errh = np.ma.zeros((self.nz, self.ny, self.nx))
            for z, d in enumerate(self.data): errh[z, :, :] = d.errh
        elif i is None and j is not None:
            errh = np.ma.zeros((self.nz, 1, self.nx))
            for z, d in enumerate(self.data): errh[z, 0, :] = d.errh[j, :]
        elif i is not None and j is None:
            errh = np.ma.zeros((self.nz, self.ny, 1))
            for z, d in enumerate(self.data): errh[z, :, 0] = d.errh[:, i]
        elif i is not None and j is not None:
            errh = np.ma.zeros((self.nz))
            for z, d in enumerate(self.data): errh[z] = d.errh[j, i]

        return errh
    def date1_(self):
        return [d.date1 for d in self.data]

    def date2_(self):
        return [d.date2 for d in self.data]

    def offset_(self):
        return [d.date1 for d in self.data]


# ====== = ====== REDUCE DATA ====== = ======
    def pick_time2(self, startday, endday):  # {{{
        '''
        :param startyear: clip from date_central = 'dd-mm-yyyy' or datetime.date variable
        :param endyear: clip to date_central = 'dd-mm-yyyy' or datetime.date variable
        '''
        print('Original cube z-size: ', len(self.data))
        print('Selecting data from {0} to {1} ..'.format(startday, endday))

        if type(startday)==str:
            try:
                startday = datetime.strptime(startday, '%d-%m-%Y').date()
                endday = datetime.strptime(endday, '%d-%m-%Y').date()
            except:
                startday = datetime.strptime(startday, '%Y-%m-%d').date()
                endday = datetime.strptime(endday, '%Y-%m-%d').date()
        
        #newcube = cube_class()
        newcube = self.copy_structure()
        newcube.filedir = self.filedir
        newcube.filename = self.filename.replace('.nc', '_' + startday.strftime('%d%m%Y') + 'to' + endday.strftime('%d%m%Y') + '.nc')
        # newcube.x = self.x
        # newcube.y = self.y
        # newcube.i = self.i
        # newcube.j = self.j
        # newcube.nx = self.nx
        # newcube.ny = self.ny
        # newcube.nz = self.nz
        # newcube.proj4 = self.proj4
        # newcube.srs = self.srs

        for d in self.data:
            if (startday <= mjd2date(d.date1 + d.offset / 2) <= endday):
                newcube.data.append(d)
        
        newcube.nz = len(newcube.data)

        print('Final cube z-size: ',len(newcube.data))
        return newcube

    def pick_time(self, startday, endday):  # {{{
        '''
        :param startyear: clip from date_central='dd-mm-yyyy' or datetime.date variable
        :param endyear: clip to date_central='dd-mm-yyyy' or datetime.date variable
        :param savedir: if you want to change the source data directory (the stored file will have a different name than
                        original so no overwriting anyway
        '''

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

    def pick_offset(self, offset_MIN, offset_MAX):  # {{{
        print('Original cube z-size: ', len(self.data))
        print('{0}-{1} days offset layers picking...'.format(offset_MIN, offset_MAX))

        i = 0
        while i < len(self.data):
            if not (offset_MIN <= self.offset_()[i] <= offset_MAX):
                self.data.pop(i)
            else:
                i = i + 1
        if i==0:
            print('No layers in this cube with time offset {0}-{1} days'.format(offset_MIN, offset_MAX))
        else:
            print(len(self.data), ' layers were picked.')

        self.nz = len(self.data)
        self.filename = '_'.join(self.filename.split('_')[0:3]) + '_{0}-{1}d.nc'.format(offset_MIN, offset_MAX)  # }}}

    def pick_sensor(self, sensor):  # {{{
        from velocity_colormap import sensors_colormap

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


    def sub_cube(self,ii=None, jj=None):
        ''' Extract from the cube a zone [i1:i2, j1:j2, z] INCLUDING THE -ii-,-jj- EDGES in NORMAL enumerate way.
        1st index = 1 (not standard Python array mode).
        [6,10] -> from 6th element to 10th including, 5 items in total
        Syntax : subc=c.sub_cube([6,10],[25,200])'''

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
        h = self.h_()[:, jj[0]:jj[1], ii[0]:ii[1]]
        errh = self.errh_()[:, jj[0]:jj[1], ii[0]:ii[1]]

        sub.insert(h=h,errh=errh,
                   date1=self.date1_(), source=self.source_(), sensor=self.sensor_())

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

        print ('New cube shape nz-ny-nx : ', sub.shape())
        return sub

    def sieve_empty_layers(self, sieve_n=500):
        ''' Delete the layers where number of valid pixels is less than sieve_n. '''
        nz_orig = self.nz
        print('Original cube.nz ', nz_orig)
        tot = self.nx * self.ny

        i = 0
        while i < len(self.data):
            if tot - np.count_nonzero(self.data[i].h.mask) < sieve_n:
                self.data.pop(i)
            else:
                i = i + 1

        self.nz = len(self.data)
        self.filename = self.filename.replace('.nc','_sieve{0}.nc'.format(sieve_n))
        if self.nz == 0:
            print (' /!\ Empty cube is returned')
        else: print(len(self.data), ' layers were picked.')

# ====== = ====== UTILITIES ====== = ======
    def coord2pix(self, x, y):  # {{{
        ''' Returm local cube pixels in python index mode [0:249]'''
        # manage memory rounding effects
        x = x + 0.005
        y = y + 0.005
        # convert
        i = int(np.ceil((x - self.geo.xmin) / self.geo.xposting)) - 1  # [1,250] -> python index mode [0:249]
        j = int(np.ceil((y - self.geo.ymax) / self.geo.yposting))
        # check
        m = 1
        if not 0 <= i <= self.nx-1: print('Given -x- is out of the range of the cube'); m = -1
        if not 0 <= j <= self.ny-1: print('Given -y- is out of the range of the cube'); m = -1
        if m == 1:
            return i, j
        else:
            return ValueError, ValueError # }}}

    def pix2coord(self, i, j, center=False):
        if not center:  # get the pixel left top corner
            return self.x[i], self.y[j]
        else:  # get the pixel center
            return self.x[i] + round(self.geo.xposting / 2, 0), self.y[j] + round(self.geo.yposting / 2, 0)

    def search_cube(self, x, y, greenland=False):
        ''' Search the cube that contains the given -x-y- point. Use only for non-merged spacely cubes! (=> -nx-ny- of a cube == -nx-ny of the greed)
           :param x, y: point geo coordinates, numbers
           :param greenland: the default geo_param is loaded, so the method could be used without any pre-loaded cube (empty -self-)
           :return name_pattern: "x00000_y00000", string '''
        if self.filename.split('_')[0]!='c':
            print(' /!\ Use the -search_cube- only for non-merged spacely cubes! (=> -nx-ny- of a cube == -nx-ny of the greed)')

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

    def deepcopy(self):
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

        for n in range(self.nz):
            d = copy.deepcopy(self.data[n])
            cube2.data.append(d)

        return cube2



# ========= = ========== = ========== = ======== = ========== = ==========
class cube_element:  # class cube_element {{{

    def __init__(self, nx=250, ny=250):
        self.h = np.ma.array(np.full([ny, nx], np.nan, dtype=np.float32))
        self.errh = np.ma.array(np.full([ny, nx], np.nan, dtype=np.float32))
        self.date1 = 0
        self.date2 = 0
        self.offset = 0
        self.source = ''  # UNIQ IDENTIFIER
        self.sensor = ''
        self.creation_date = 0.   # }}}

