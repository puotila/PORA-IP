#!/usr/bin/env python
"""
Plot Polar ORA-IP annual mean profiles (T and S),
and Hiroshi Sumata's and WOA13 1995-2004 profiles.
"""

import os
import sys
import re
import copy
import glob
import string
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime
from netcdftime import utime
from seawater import dens0

class ProfVar(object):
    """ Temperature (T), salinity (S) or density (R) vertical
        profile variable.
    """
    def __init__(self,name,level_bounds):
        self.name = name
        self.level_bounds = level_bounds
        nz = level_bounds.shape[0]
        self._fillValue = 99
        self.data = np.ma.masked_equal(self._fillValue*np.ones((nz,)),self._fillValue)
        # upper level depth
        self.uz = level_bounds[:,0]
        # lower level depth
        self.lz = level_bounds[:,1]
        # middle level depth
        self.mz = np.mean(level_bounds,axis=1)

class Product(object):
    """ ORAIP annual means of T and S for the basin-average profile.
    """
    def __init__(self,basin='Antarctic',\
                 syr=1993,eyr=2010,\
                 path='/home/uotilap/tiede/ORA-IP/annual_mean/',\
                 level_bounds=None):
        self.path = path
        self.basin = basin
        self.syr, self.eyr = syr, eyr
        self.LevelBounds = {}
        if level_bounds is None:
            self.LevelBounds['T'] = np.array([[   0, 100],[ 100,  300],\
                                              [ 300, 700],[ 700, 1500],\
                                              [1500,3000]])
            self.LevelBounds['S'] = np.array([[   0, 100],[ 100,  300],\
                                              [ 300, 700],[ 700, 1500],\
                                              [1500,3000]])
        else:
            self.LevelBounds = level_bounds
        self.bathymetry = self.readWOA13Bathymetry()
        for vname in ['T','S']:
            setattr(self,vname,ProfVar(vname,self.LevelBounds[vname]))
        self.nclatname, self.nclonname    = 'lat', 'lon'
        self.nctimename, self.ncdepthname = 'time', 'depth'
        self.linestyle = '-'
        self.lettercolor = 'white'
        self.edgecolor = 'black'

    def getNetCDFfilename(self,varname,ulb,llb):
        return self.fpat % (varname,self.dsyr,self.deyr,ulb,llb)

    def getNetCDFfilepointer(self,fn,path=None):
        if path is None:
            path = self.path
        if not os.path.exists(path):
            # we are on voima
            path = '/lustre/tmp/uotilap/ORA-IP/annual_mean'
        if not os.path.exists(path):
            # last resort is cwd
            path = './'
        try:
           fp = nc.Dataset(os.path.join(path,fn))
           print "Reading %s" % os.path.join(path,fn)
        except:
           print "Cant read %s!" % os.path.join(path,fn)
           sys.exit(0)
        return fp

    def getFillValue(self,ncvar):
        if hasattr(ncvar,'_FillValue'):
            FillValue = ncvar._FillValue
        elif hasattr(ncvar,'missing_value'):
            FillValue = ncvar.missing_value
        else:
            FillValue = None
        return FillValue

    def readVarProfile(self,varname):
        """ varname is either T or S
        """
        data = []
        for li, lb in enumerate(self.LevelBounds[varname]):
            fn    = self.getNetCDFfilename(varname,0,lb[1])
            ldata = self.readOneFile(fn,self.ncvarname[varname],lb[1])
            if lb[0]==0.:
                udata = 0.0*ldata
            else:
                fn    = self.getNetCDFfilename(varname,0,lb[0])
                udata = self.readOneFile(fn,self.ncvarname[varname],lb[0])
            data.append((ldata - udata)/(lb[1] - lb[0])) # [t,y,x] variable values from level averages
        return np.ma.array(data) # [z,t,y,x]

    def maskBadSalinity(self,data):
        # get rid of bad mdata values
        # basically if S in the lower layer is smaller than in the upper one,
        # mask it.
        for z in range(1,data['S'].shape[0]):
            zmask = np.ma.make_mask(data['S'][z]<data['S'][z-1])
            for varname in ['S','T']:
                zdata = data[varname][z]
                dmask = np.ma.mask_or(zdata.mask,zmask)
                data[varname][z] = np.ma.array(zdata,mask=dmask)
        return data

    def averageBasinAndTime(self,data):
        """ data needs to be [z,t,y,x]
        """
        for varname in ['S','T']:
            pdata = getattr(self,varname)
            pdata.data = np.ma.mean(data[varname],\
                         axis=tuple(range(1,data[varname].ndim))) # time and basin average

    def readProfile(self):
        """ varname is either T or S
        """
        data = {'T':[],'S':[]} # [z,t,y,x]
        for varname in ['S','T']:
            data[varname] = self.readVarProfile(varname)
        data = self.maskBadSalinity(data)
        self.averageBasinAndTime(data)

    def readLatLon(self,fp):
        lat = np.array(fp.variables[self.nclatname][:])
        lon = np.array(fp.variables[self.nclonname][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        return lon, lat

    def getDates(self,fp):
        time = fp.variables[self.nctimename]
        if hasattr(time,'calendar'):
            cdftime = utime(time.units,calendar=time.calendar.lower())
        else:
            cdftime = utime(time.units)
        return [cdftime.num2date(t) for t in time[:]]

    def findClosestLocation(self,lon,lat):
        ix = np.where(np.abs(lon-self.plon)==np.min(np.abs(lon-self.plon)))[0][0]
        iy = np.where(np.abs(lat-self.plat)==np.min(np.abs(lat-self.plat)))[0][0]
        return ix, iy

    def findBasinIndex(self, lon, lat):
        lon, lat  = np.meshgrid(lon, lat)
        # mask which is one in the basin and masked elsewhere
        # so a field multiplied by the mask retains its values
        # in the basin and gets masked outside
        basinmask = np.ma.masked_all(lat.shape)
        if self.basin=='Antarctic':
            # Antarctic shelf/deep regions
            iy, ix = np.where(((lon>330) & (lon<=360) & (lat<=-60)) | ((lon>0) & (lon<=35) & (lat<=-60))    |
                              ((lon>35) & (lon<=68) & (lat<=-61)) | ((lon>68) & (lon<=95) & (lat<=-60))     |
                              ((lon>95) & (lon<=110) & (lat<=-62)) | ((lon>110) & (lon<=160) & (lat<=-64))  |
                              ((lon>160) & (lon<=235) & (lat<=-66)) | ((lon>235) & (lon<=280) & (lat<=-68)) |
                              ((lon>280) & (lon<=300) & (lat<=-66)) | ((lon>300) & (lon<=315) & (lat<=-64)) |
                              ((lon>315) & (lon<=330) & (lat<=-62)))
        elif self.basin=='Arctic':
            iy, ix = np.where(lat>80)
        elif self.basin=='Nansen':
            iy, ix = np.where(((lon<135) & (lat>80))| ((lon>315) & (lat>80)))
        elif self.basin=='Canadian':
            iy, ix = np.where((lon>=135) & (lon<=315) & (lat>80))
        else:
            print "%s basin has not been defined!" % self.basin
            sys.exit(0)
        basinmask[iy,ix]=1
        return basinmask

    def readWOA13Bathymetry(self,bfile='landsea_01.msk'):
        """ Can be downloaded from
            https://www.nodc.noaa.gov/OC5/woa13/masks13.html
        """
        depth = range(0,105,5)+range(125,525,25)+range(550,2050,50)+range(2100,9200,100)
        dat = np.loadtxt(bfile,skiprows=2,delimiter=',')
        b2d = np.reshape([depth[int(i)-1] for i in dat[:,2]],(180,360))
        return np.ma.masked_values(b2d,0)

    def readOneFile(self,fn,ncvarname,maxdpth=0.):
        """
        Read data from a netCDF file and return its temporal mean
        within given year range [syr, eyr] for the basin average.
        """
        fp = self.getNetCDFfilepointer(fn)
        lon, lat = self.readLatLon(fp)
        basinmask = self.findBasinIndex(lon,lat)
        # mask too shallow regions from depth integrals as their
        # values are too small
        bathymask = np.ma.make_mask(self.bathymetry<maxdpth)
        # combine basin and bathymasks
        fldmask = np.ma.mask_or(basinmask.mask,bathymask)
        dates    = self.getDates(fp)
        ldata = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                ncvar = fp.variables[ncvarname]
                FillValue = self.getFillValue(ncvar)
                data = np.ma.squeeze(np.ma.array(ncvar[i]))
                if FillValue is not None:
                    data  = np.ma.masked_values(data,FillValue)
                data = np.ma.array(data,mask=fldmask)
                ldata.append(data) # do not average across the basin
        fp.close()
        return np.ma.array(ldata)

    def getLayeredDepthProfile(self,varname,depth,data):
        """
        Average 3D hires profile data according to level_bounds
        """
        if data.ndim>1:
            print "Warning: basin average has %d dimensions!" % data.ndim
        ldata    = []
        for li, lb in enumerate(self.LevelBounds[varname]):
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            data1 = data[iz]
            data2 = data1[~np.isnan(data1)]
            ldata.append(np.ma.mean(data2))
        return np.ma.array(ldata,mask=np.isnan(ldata))

class CGLORS(Product):
    def __init__(self,basin,syr,eyr):
        super( CGLORS, self).__init__(basin,syr,eyr)
        self.dset = 'CGLORS'
        self.dsyr, self.deyr = 1989, 2014
        self.fpat = 'CGLORS025v5/%s_ORCA025-%s_%d.nc'
        self.ncvarname = {'T':'votemper',\
                          'S':'vosaline'}
        self.linecolor = self.scattercolor = 'lightgreen'
        self.legend = 'C-GLORS025v5'

    def getDates(self,fp,year,months=range(1,13)):
        return [datetime(year,month,15) for month in months]

    def getNetCDFfilename(self,varname,year):
        if year in range(1989,1993):
            grid, ncdepthname = 'WOA', 'deptht'
        else:
            grid, ncdepthname = '1x1', 'dep'
        if varname=='T' and year==2010:
            grid, ncdepthname = 'WOA', 'deptht'
        fn = self.fpat % (self.ncvarname[varname],grid,year)
        return fn, ncdepthname

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin averaged values.
        """
        for varname in ['S','T']:
            pdata = getattr(self,varname)
            ncvarname = self.ncvarname[varname]
            tdata = []
            for year in range(self.syr,self.eyr+1):
                fn, ncdepthname = self.getNetCDFfilename(varname,year)
                fp = self.getNetCDFfilepointer(fn)
                lon, lat = self.readLatLon(fp)
                dates    = self.getDates(fp,year)
                depth    = np.array(fp.variables[ncdepthname])
                basinmask = self.findBasinIndex(lon,lat)
                ncvar    = fp.variables[ncvarname]
                FillValue = self.getFillValue(ncvar)
                for i,date in enumerate(dates[:]):
                    if date.year in range(self.syr,self.eyr+1):
                       data = np.ma.array(fp.variables[ncvarname][i])
                       if FillValue is None:
                           data *= basinmask
                       else:
                           data  = np.ma.masked_values(data, FillValue)*basinmask
                       # basin average
                       data_ba = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                       tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
                fp.close()
            pdata.data = np.ma.mean(tdata,axis=0) # temporal average

class ECDA(Product):
    def __init__(self,basin,syr,eyr):
        super( ECDA, self).__init__(basin,syr,eyr)
        self.dset = 'ECDA'
        self.dsyr, self.deyr = 1993, 2011
        self.fpat = "ECDA_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'yellow'
        self.lettercolor = 'black'
        self.legend = 'ECDA3'

    def getDates(self,fp):
        if fp.variables.has_key(self.nctimename):
            time = fp.variables[self.nctimename]
        else:
            time = fp.variables[self.nctimename.upper()]
        if hasattr(time,'calendar'):
            cdftime = utime(time.units,calendar=time.calendar.lower())
        else:
            cdftime = utime(time.units)
        return [cdftime.num2date(t) for t in time[:]]

class GloSea5(Product):
    def __init__(self,basin,syr,eyr):
        super( GloSea5, self).__init__(basin,syr,eyr)
        self.dset = 'GloSea5'
        self.dsyr, self.deyr = 1993, 2014
        self.fpat = "GloSea5_GO5_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'blue'
        self.nctimename = 'time_counter'
        self.legend = 'GloSea5-GO5'

    def readLatLon(self,fp):
        lat = np.arange(-89.5,90.5,1.)
        lon = np.arange(0.5,360.5,1.)
        return lon, lat

class MOVEG2(Product):
    def __init__(self,basin,syr,eyr):
        super( MOVEG2, self).__init__(basin,syr,eyr)
        self.dset = self.legend = 'MOVEG2'
        self.dsyr, self.deyr = 1993, 2012
        self.fpat = "MOVEG2_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'cyan'
        self.lettercolor = 'black'

    def getDates(self,fp):
        time = fp.variables[self.nctimename]
        if hasattr(time,'calendar'): #T
            cdftime = utime(time.units,calendar=time.calendar)
            dates = [cdftime.num2date(t) for t in time[:]]
        else: #S
            m = re.search('months since\s+(\d+)-(\d+)-(\d+)',time.units)
            year0, month0, day0 = [int(s) for s in m.groups()]
            dates = [datetime(year0+int(t/12),month0+int(t%12),day0) for t in time[:]]
        return dates

class UoR(Product):
    def __init__(self,basin,syr,eyr):
        super( UoR, self).__init__(basin,syr,eyr)
        self.dset = 'UoR'
        self.dsyr, self.deyr = 1989, 2010
        self.fpat = "UoR_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'lightblue'
        self.lettercolor = 'black'
        self.legend = 'UR025.4'

class EN4(Product):
    def __init__(self,basin,syr,eyr):
        super( EN4, self).__init__(basin,syr,eyr)
        self.dset  = 'EN4'
        self.legend = 'EN4.2.0.g10'
        self.dsyr, self.deyr = 1950, 2015
        self.fpat = "EN4.2.0.g10_int%s_annmean_%dto%d_%d-%dm.nc"
        self.ncvarname = {'T':'t_int_',\
                          'S':'s_int_'}
        self.linecolor = self.scattercolor = 'pink'
        self.bathymetry = self.bathymetry[7:,:]

    def readProfile(self):
        """ varname is either T or S
        """
        data = {'T':[],'S':[]} # [z,t,x,y]
        for varname in ['S','T']:
            mdata = []
            for li, lb in enumerate(self.LevelBounds[varname]):
                ncvarname = "%s%d" % (self.ncvarname[varname],lb[1])
                fn    = self.getNetCDFfilename(varname,0,lb[1])
                ldata = self.readOneFile(fn,ncvarname,lb[1])
                if lb[0]==0.:
                    udata = 0.0*ldata
                else:
                    ncvarname = "%s%d" % (self.ncvarname[varname],lb[0])
                    fn    = self.getNetCDFfilename(varname,0,lb[0])
                    udata = self.readOneFile(fn,ncvarname,lb[0])
                mdata.append((ldata - udata)/(lb[1] - lb[0])) # [t,x,y] variable values from level averages
            data[varname] = np.ma.array(mdata) # [z,t,x,y]
        data = self.maskBadSalinity(data)
        self.averageBasinAndTime(data)

class GECCO2(Product):
    def __init__(self,basin,syr,eyr):
        super( GECCO2, self).__init__(basin,syr,eyr)
        self.dset  = 'GECCO2'
        self.dsyr, self.deyr = 1948, 2012
        self.fpat = "GECCO2_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'S_0_%d'}
        self.linecolor = self.scattercolor = 'darkred'
        self.legend = 'GECCO2'
        self.FillValue = -1e34

    def getGECCO2SalinityDates(self,fp):
        return [datetime(year,1,1) for year in range(self.dsyr,self.deyr+1)]

    def readOneSalinityField(self,fp,varname,maxdpth,basinmask,i):
        # mask too shallow regions from depth integrals as their
        # values are too small
        bathymask = np.ma.make_mask(self.bathymetry<maxdpth)
        # combine basin and bathymasks
        fldmask = np.ma.mask_or(basinmask.mask,bathymask)
        ncvarname = self.ncvarname[varname] % maxdpth
        ncvar = fp.variables[ncvarname]
        data = np.ma.squeeze(np.ma.array(ncvar[i]))
        ldata  = np.ma.array(data,mask=fldmask)
        return np.ma.masked_values(ldata,self.FillValue) # do not average across the basin

    def readGECCO2Profile(self):
        data = {'T':[],'S':[]} # [z,t,y,x]
        data['S'] = self.readGECCO2SalinityProfile()
        # read temperature profile
        data['T'] = self.readVarProfile('T')
        data = self.maskBadSalinity(data)
        self.averageBasinAndTime(data)

    def readGECCO2TemperatureProfile(self,varname='T'):
        tdata = [] # [z,t,y,x]
        for li, lb in enumerate(self.LevelBounds[varname]):
            ncvarname = self.ncvarname[varname]
            fn    = self.getNetCDFfilename(varname,0,lb[1])
            ldata = self.readOneFile(fn,ncvarname,lb[1])
            if lb[0]==0.:
                udata = 0.0*ldata
            else:
                ncvarname = self.ncvarname[varname]
                fn    = self.getNetCDFfilename(varname,0,lb[0])
                udata = self.readOneFile(fn,ncvarname,lb[0])
            tdata.append((ldata - udata)/(lb[1] - lb[0])) # [t,x,y] variable values from level averages
        return np.ma.array(tdata)

    def readGECCO2SalinityProfile(self,varname='S'):
        fn = 'GECCO2_intS_annmean_1948to2011_all_layers_r360x180.nc'
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        dates     = self.getGECCO2SalinityDates(fp)
        basinmask = self.findBasinIndex(lon,lat)
        pdata     = getattr(self,varname)
        tdata     = []
        for i,date in enumerate(dates):
            if date.year in range(self.syr,self.eyr+1):
                ldata = []
                for li, lb in enumerate(self.LevelBounds[varname]):
                    lldata = self.readOneSalinityField(fp,varname,lb[1],basinmask,i)
                    if lb[0]==0.:
                        ludata = 0.0*lldata
                    else:
                        ludata = self.readOneSalinityField(fp,varname,lb[1],basinmask,i)
                    ldata.append((lldata*lb[1] - ludata*lb[0])/(lb[1] - lb[0])) # [z,y,x]
                tdata.append(ldata) # [t,z,y,x]
        fp.close()
        tdata = np.ma.masked_values(tdata,self.FillValue)
        return np.ma.swapaxes(tdata,0,1) # do not temporal average, -> [z,t,y,x]

class GLORYS2V4(Product):
    def __init__(self,basin,syr,eyr):
        super( GLORYS2V4, self).__init__(basin,syr,eyr)
        self.dset  = 'GLORYS'
        self.fpat  = 'GSOP_GLORYS2V4_ORCA025_%s.nc'
        self.ncvarname = {'T':'z%dheatc',\
                          'S':'z%dsaltc'}
        self.ncdepthname = 'zdepth'
        self.linecolor = self.scattercolor = 'orange'
        self.legend = 'GLORYS2v4'

    def getNetCDFfilename(self,varname):
        if varname=='T':
            fn = self.fpat % ('HC')
        else:
            fn = self.fpat % ('SC')
        return fn

    def readOneField(self,fp,varname,maxdpth,basinmask,i):
        # mask too shallow regions from depth integrals as their
        # values are too small
        bathymask = np.ma.make_mask(self.bathymetry<maxdpth)
        # combine basin and bathymasks
        fldmask = np.ma.mask_or(basinmask.mask,bathymask)
        ncvarname = self.ncvarname[varname] % maxdpth
        ncvar = fp.variables[ncvarname]
        FillValue = self.getFillValue(ncvar)
        data = np.ma.squeeze(np.ma.array(ncvar[i]))
        if FillValue is not None:
           data  = np.ma.masked_values(data,FillValue)
        ldata  = np.ma.array(data,mask=fldmask)
        return np.ma.array(ldata) # do not average across the basin

    def readVarProfile(self,varname):
        """ varname is either T or S
        """
        fn = self.getNetCDFfilename(varname)
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        dates     = self.getDates(fp)
        depth     = np.array(fp.variables[self.ncdepthname])
        basinmask = self.findBasinIndex(lon,lat)
        data = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                mdata = []
                for li, lb in enumerate(self.LevelBounds[varname]):
                    lldata = self.readOneField(fp,varname,lb[1],basinmask,i)
                    if lb[0]==0.:
                        ludata = 0.0*lldata
                    else:
                        ludata = self.readOneField(fp,varname,lb[0],basinmask,i)
                    mdata.append((lldata - ludata)/(lb[1] - lb[0]))
                data.append(np.ma.array(mdata)) # [t,z,y,x]
        fp.close()
        return np.ma.swapaxes(data,0,1) # [t,z,y,x] -> [z,t,y,x]

class TOPAZ(Product):
    def __init__(self,basin,syr,eyr):
        super( TOPAZ, self).__init__(basin,syr,eyr)
        self.dset  = 'TOPAZ'
        self.dsyr, self.deyr = 1993, 2013
        self.fpat  = "TP4_r360x180_%s_%04d_%02d.nc"
        self.ncvarname = {'T':'temperature',\
                          'S':'salinity'}
        self.nclatname = 'latitude'
        self.nclonname = 'longitude'
        self.linecolor = self.scattercolor = 'green'
        self.legend    = 'TOPAZ4'

    def getNetCDFfilename(self,varname,year,month):
        if varname=='T':
            fn = self.fpat % ('temp',year,month)
        else:
            fn = self.fpat % ('salt',year,month)
        return fn

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin-averaged profile.
        """
        for varname in ['S','T']:
            pdata = getattr(self,varname)
            ncvarname = self.ncvarname[varname]
            tdata = []
            for year in range(self.syr,self.eyr+1):
                for month in range(1,13):
                    fn = self.getNetCDFfilename(varname,year,month)
                    fp = self.getNetCDFfilepointer(fn)
                    lon, lat = self.readLatLon(fp)
                    depth    = np.array(fp.variables[self.ncdepthname])
                    basinmask= self.findBasinIndex(lon,lat)
                    ncvar    = fp.variables[ncvarname]
                    FillValue = self.getFillValue(ncvar)
                    data     = np.ma.array(fp.variables[ncvarname])
                    data     = np.ma.masked_values(data, FillValue)*basinmask
                    data_ba  = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                    tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
            pdata.data = np.ma.mean(tdata,axis=0)

class MultiModelMean(Product):
    def __init__(self,basin):
        super( MultiModelMean, self).__init__(basin)
        self.dset = self.legend = 'MMM'
        self.linecolor = self.scattercolor = self.edgecolor = 'lightgrey'
        self.linestyle = ':'

    def calcMultiModelMean(self,products,vname):
        setattr(getattr(self,vname),'data',\
                np.ma.mean([getattr(getattr(p,vname),'data')\
                for p in products],axis=0))

class Sumata(Product):
    def __init__(self,basin,syr=1980,eyr=2015):
        super( Sumata, self).__init__(basin,syr,eyr)
        self.dset = self.legend = 'Sumata'
        self.dsyr, self.deyr = syr, eyr
        self.fpat = "ts-clim/hiroshis-clim/archive_v12_QC2_3_DPL_checked_2d_season_all-remapbil-oraip.nc"
        self.ncvarname = {'T':'temperature',\
                          'S':'salinity'}
        self.linecolor = self.scattercolor = 'black'

    def getNetCDFfilename(self):
        return self.fpat

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin-averaged profile.
        """
        fn = self.getNetCDFfilename()
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        depth     = np.array(fp.variables[self.ncdepthname])
        basinmask = self.findBasinIndex(lon,lat)
        for varname in ['S','T']:
            ncvarname = self.ncvarname[varname]
            ncvar     = fp.variables[ncvarname]
            FillValue = self.getFillValue(ncvar)
            pdata = getattr(self,varname)
            tdata = []
            # if seasonal average is not representative we may
            # need to plot seasons separately
            for i in range(4):
                data    = np.ma.masked_values(ncvar[i],FillValue)*basinmask
                data_ba = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
            pdata.data = np.ma.mean(tdata,axis=0)
        fp.close()

class WOA13(Sumata):
    def __init__(self,basin,syr=1995,eyr=2012):
        super( WOA13, self).__init__(basin,syr,eyr)
        self.dset = self.legend = 'WOA13'
        self.fpat = "ts-clim/woa13/woa13-clim-1995-2012-season.nc"
        self.scattercolor = 'white'
        self.linestyle = '--'
        self.lettercolor = 'black'

class ORAP5(Product):
    def __init__(self,basin,syr,eyr):
        super( ORAP5, self).__init__(basin,syr,eyr)
        self.dset  = 'ORAP5'
        self.dsyr, self.deyr = 1993, 2012
        self.fpat  = "%s3D_orap5_1m_%d-%d_r360x180.nc"
        self.ncvarname = {'T':'votemper',\
                          'S':'vosaline'}
        self.ncdepthname = 'deptht'
        self.nctimename = 'time_counter'
        self.linecolor = self.scattercolor = 'red'
        self.legend = 'ORAP5'

    def getDates(self,fp):
        time = fp.variables[self.nctimename]
        m = re.search('month.*since\s+(\d+)-(\d+)-(\d+)',time.units)
        year0, month0, day0 = [int(s) for s in m.groups()]
        return [datetime(year0+int(t/12),month0+int(t%12),day0) for t in time[:]]

    def getNetCDFfilename(self,varname):
        if varname=='T':
            fn = self.fpat % ('temperature',self.dsyr,self.deyr)
        else:
            fn = self.fpat % ('salinity',self.dsyr,self.deyr)
        return fn

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin-averaged profile.
        """
        for varname in ['S','T']:
            fn = self.getNetCDFfilename(varname)
            fp = self.getNetCDFfilepointer(fn)
            lon, lat  = self.readLatLon(fp)
            dates     = self.getDates(fp)
            depth     = np.array(fp.variables[self.ncdepthname])
            basinmask = self.findBasinIndex(lon,lat)
            ncvarname = self.ncvarname[varname]
            ncvar     = fp.variables[ncvarname]
            FillValue = self.getFillValue(ncvar)
            pdata     = getattr(self,varname)
            tdata     = []
            for i,date in enumerate(dates[:]):
                if date.year in range(self.syr,self.eyr+1):
                    data    = np.ma.masked_values(ncvar[i],FillValue)*basinmask
                    data_ba = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                    tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
            fp.close()
            pdata.data = np.ma.mean(tdata,axis=0)

class MOVEG2i(Product):
    def __init__(self,basin,syr,eyr):
        super( MOVEG2i, self).__init__(basin,syr,eyr)
        self.dset  = 'MOVEG2i'
        self.dsyr, self.deyr = 1980, 2012
        self.fpat  = "MOVEG2i_%s3d_%d-%d.nc"
        self.ncvarname = {'T':'temp',\
                          'S':'sal'}
        self.ncdepthname = 'level'
        self.nctimename = 'time'
        self.linecolor = self.scattercolor = 'cyan'
        self.legend = 'MOVE-G2i'

    def getDates(self,fp):
        time = fp.variables[self.nctimename]
        m = re.search('month.*since\s+(\d+)-(\d+)-(\d+)',time.units)
        year0, month0, day0 = [int(s) for s in m.groups()]
        return [datetime(year0+int(t/12),month0+int(t%12),day0) for t in time[:]]

    def getNetCDFfilename(self,varname):
        if varname=='T':
            fn = self.fpat % ('temp',self.dsyr,self.deyr)
        else:
            fn = self.fpat % ('sal',self.dsyr,self.deyr)
        return fn

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin-averaged profile.
        """
        for varname in ['S','T']:
            fn = self.getNetCDFfilename(varname)
            fp = self.getNetCDFfilepointer(fn)
            lon, lat  = self.readLatLon(fp)
            dates     = self.getDates(fp)
            depth     = np.array(fp.variables[self.ncdepthname])
            basinmask = self.findBasinIndex(lon,lat)
            ncvarname = self.ncvarname[varname]
            ncvar     = fp.variables[ncvarname]
            FillValue = self.getFillValue(ncvar)
            pdata     = getattr(self,varname)
            tdata     = []
            for i,date in enumerate(dates[:]):
                if date.year in range(self.syr,self.eyr+1):
                    data    = np.ma.masked_values(ncvar[i],FillValue)*basinmask
                    # Note that data is (nz,ny,nx) and basinmask (ny,nx)
                    # their multiplication data*basinmask returns (nz,ny,nz) where
                    # each data[nz] is multiplied by basinmask, clever eh?
                    data_ba = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                    tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
            fp.close()
            pdata.data = np.ma.mean(tdata,axis=0)

class SODA331(Product):
    def __init__(self,basin,syr,eyr):
        super( SODA331, self).__init__(basin,syr,eyr)
        self.dset = 'SODA3.3.1'
        self.dsyr, self.deyr = 1980, 2015
        self.fpat = '%s3D_SODA3.3.1/%s3D_SODA_3_3_1_%d.nc'
        self.ncvarname = {'T':'temp',\
                          'S':'salt'}
        self.nclatname = 'latitude'
        self.nclonname = 'longitude'
        self.ncdepthname = 'depth'
        self.linecolor = self.scattercolor = 'purple'
        self.legend = 'SODA3.3.1'

    def getDates(self,fp,year,months=range(1,13)):
        return [datetime(year,month,1) for month in months]

    def getNetCDFfilename(self,varname,year):
        if varname=='T':
            fn = self.fpat % ('temperature','temperature',year)
        else:
            fn = self.fpat % ('salinity','salinity',year)
        return fn

    def readProfile(self):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        for the basin-averaged profile.
        """
        for varname in ['S','T']:
            pdata = getattr(self,varname)
            ncvarname = self.ncvarname[varname]
            tdata = []
            for year in range(self.syr,self.eyr+1):
                fn = self.getNetCDFfilename(varname,year)
                fp = self.getNetCDFfilepointer(fn)
                lon, lat = self.readLatLon(fp)
                dates    = self.getDates(fp,year)
                depth    = np.array(fp.variables[self.ncdepthname])
                basinmask= self.findBasinIndex(lon,lat)
                ncvar    = fp.variables[ncvarname]
                FillValue = self.getFillValue(ncvar)
                for i,date in enumerate(dates[:]):
                    if date.year in range(self.syr,self.eyr+1):
                       data    = np.ma.masked_values(ncvar[i],FillValue)*basinmask
                       data_ba = np.ma.mean(data,axis=tuple(range(1, data.ndim)))
                       tdata.append(self.getLayeredDepthProfile(varname,depth,data_ba))
                fp.close()
            pdata.data = np.ma.mean(tdata,axis=0)

class Products(object):
    """ Container for ORA-IP products
    """
    def __init__(self,productobjs,basin='Antarctic',\
                 syr=1993,eyr=2010):
        self.products = []
        self.syr, self.eyr = syr, eyr
        self.title = self.basin = basin
        for pobj in productobjs:
            self.products.append(pobj(basin,syr,eyr))
        self.mmm = MultiModelMean(basin)
        if basin not in ['Antarctic']:
            self.sumata = Sumata(basin)
        self.woa13 = WOA13(basin)
        self.xlabel = {'T':"Temperature [$^\circ$C]",\
                       'S':"Salinity [psu]"}
        self.ylabel = "depth [m]"
        self.ymin   = 3100 # m
        # Products per 3 panels:
        self.ProductPanels = {'CGLORS':0,'GECCO2':0,'GLORYS':0,'GloSea5':0,\
                              'ORAP5':1,'SODA3.3.1':1,'TOPAZ':1,'UoR':1,\
                              'ECDA':2,'MOVEG2i':2,'EN4':2}
        self.pretitle = ['(a)','(b)','(c)']
        modstr = '_'.join([p.dset for p in self.products])
        self.fileout = "%s_%04d-%04d_%s" % \
                       (modstr,syr,eyr,basin)

    def readProfiles(self):
        """ Reads now both T and S profiles
        """
        for product in self.products:
            if product.dset=='GECCO2':
                product.readGECCO2Profile()
            else:
                product.readProfile()
        if self.basin not in ['Antarctic']:
            self.sumata.readProfile()
        self.woa13.readProfile()

    def getMultiModelMean(self,vname):
        """ EN4 is not a part of MMM!
        """
        products = [product for product in self.products if product.dset not in ['EN4']]
        self.mmm.calcMultiModelMean(products,vname)

    def getDataRange(self,vname):
        if self.basin in ['Antarctic']:
            prdlist = self.products+[self.woa13]
        else:
            prdlist = self.products+[self.sumata]+[self.woa13]
        data = np.ma.hstack([np.ma.array(getattr(getattr(p,vname),'data')) \
            for p in prdlist])
        return np.ma.min(data), np.ma.max(data)

    def calcDensityMap(self):
        # Calculate how many gridcells we need in the x and y dimensions
        tmin, tmax = self.getDataRange('T')
        smin, smax = self.getDataRange('S')
        xdim = int(round((smax-smin)/0.1+1,0))
        ydim = int(round((tmax-tmin)/0.1+1,0))
        # Create empty grid of zeros
        dens = np.zeros((ydim,xdim))
        # Create temp and salt vectors of appropiate dimensions
        ti = np.linspace(1,ydim-1,ydim)*0.1+tmin
        si = np.linspace(1,xdim-1,xdim)*0.1+smin
        # Loop to fill in grid with densities
        for j in range(0,int(ydim)):
            for i in range(0, int(xdim)):
                dens[j,i]=dens0(si[i],ti[j])
        # Substract 1000 to convert to sigma-t
        return si, ti, dens - 1000

    def plotOneProfile(self,product,vname,ax):
        y = getattr(getattr(product,vname),'lz')
        x = getattr(getattr(product,vname),'data')
        lne = ax.plot(np.ma.hstack((x[0],x)),\
                      np.hstack((0,y)),\
                      lw=2,linestyle=product.linestyle,\
                      drawstyle='steps-post',color=product.linecolor)[0]
        return lne

    def getRefProductList(self):
        # First Sumata and WOA13 climatologies, and MMM
        if self.basin in ['Antarctic']:
            prdlist = [self.woa13,self.mmm]
        else:
            prdlist = [self.sumata,self.woa13,self.mmm]
        return prdlist

    def plotDepthProfile(self,vname):
        xmin, xmax = self.getDataRange(vname)
        fig = plt.figure(figsize=(8*2,10))
        axs = [plt.axes([0.10, 0.1, .2, .8]),\
               plt.axes([0.40, 0.1, .2, .8]),\
               plt.axes([0.70, 0.1, .2, .8])]
        lnes = [[],[],[]]
        lgds = [[],[],[]]
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            prdlist  = self.getRefProductList()
            for product in prdlist:
                x = getattr(getattr(product,vname),'data')
                if x.all() is np.ma.masked:
                    """ all values are masked
                    """
                    continue
                lne.append(self.plotOneProfile(product,vname,ax))
                lgd.append(product.legend)
        # then individual models
        for product in self.products:
            x = getattr(getattr(product,vname),'data')
            if x.all() is np.ma.masked:
                """ all values are masked
                """
                continue
            panelno = self.ProductPanels[product.dset]
            ax, lne, lgd = axs[panelno],lnes[panelno],lgds[panelno]
            lne.append(self.plotOneProfile(product,vname,ax))
            lgd.append(product.legend)
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            ax.invert_yaxis()
            ax.set_yticks([0,100,300,700,1500,3000])
            ax.set_ylim(self.ymin,0)
            if vname=='T':
                ax.set_xlim(xmin-np.abs(0.1*xmin),xmax+np.abs(0.1*xmax))
            else:
                ax.set_xlim(xmin-np.abs(0.01*xmin),xmax+np.abs(0.01*xmax))
                #ax.set_xlim([32, 36])
            ax.set_ylabel(self.ylabel)
            ax.set_title("%s %s" % (self.pretitle[panelno],self.title))
            ax.set_xlabel(self.xlabel[vname])
            if vname=='T':
                leg = ax.legend(lne,tuple(lgd),ncol=1,bbox_to_anchor=(0.3, 0.25))
            else:
                leg = ax.legend(lne,tuple(lgd),ncol=1,bbox_to_anchor=(0.3, 0.25))
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(1.0)
            leg.get_frame().set_alpha(1.0)
        #plt.show()
        plt.savefig('./basin_avg/'+vname+'_'+self.fileout+'.pdf')

    def plotTSProfile(self):
        fig = plt.figure(figsize=(8*2,10))
        axs = [plt.axes([0.10, 0.1, .2, .8]),\
               plt.axes([0.40, 0.1, .2, .8]),\
               plt.axes([0.70, 0.1, .2, .8])]
        lnes = [[],[],[]]
        lgds = [[],[],[]]
        si,ti,dens = self.calcDensityMap()
        scatterlw, scattersize = 1, 200
        for panelno, ax in enumerate(axs):
            CS = ax.contour(si,ti,dens, linestyles='dashed', colors='k')
            ax.clabel(CS, fontsize=12, inline=1, fmt='%2.1f') # Label every second level
            lne, lgd = lnes[panelno],lgds[panelno]
            prdlist  = self.getRefProductList()
            for product in prdlist:
                x = getattr(getattr(product,'T'),'data')
                if x.all() is np.ma.masked:
                    """ all values are masked
                    """
                    continue
                y = getattr(getattr(product,'T'),'data')
                x = getattr(getattr(product,'S'),'data')
                lne.append(ax.scatter(x,y,lw=scatterlw,s=scattersize,\
                           color=product.scattercolor,edgecolor=product.edgecolor))
                for i in range(len(y)):
                    ax.annotate("%d" % (i+1), (x[i],y[i]),\
                                va='center',ha='center',\
                                color=product.lettercolor)
                lgd.append(product.legend)
        # then individual models
        for product in self.products:
            x = getattr(getattr(product,'T'),'data')
            if x.all() is np.ma.masked:
                """ all values are masked
                """
                continue
            panelno = self.ProductPanels[product.dset]
            ax, lne, lgd = axs[panelno],lnes[panelno],lgds[panelno]
            y = getattr(getattr(product,'T'),'data')
            x = getattr(getattr(product,'S'),'data')
            lne.append(ax.scatter(x,y,lw=scatterlw,s=scattersize,\
                       color=product.scattercolor))
            for i in range(len(y)):
                ax.annotate("%d" % (i+1), (x[i],y[i]),\
                            va='center',ha='center',\
                            color=product.lettercolor)
            lgd.append(product.legend)
            # mark depth levels with numbers
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            ax.set_ylim(np.min(ti)-0.2,np.max(ti)+0.2)
            ax.set_xlim(np.min(si)-0.2,np.max(si)+0.2)
            ax.set_ylabel(self.xlabel['T'])
            ax.set_title("%s %s" % (self.pretitle[panelno],self.title))
            ax.set_xlabel(self.xlabel['S'])
            ax.legend(lne,tuple(lgd),ncol=1,\
                      bbox_to_anchor=(1.0, 0.2))
        #plt.show()
        plt.savefig('./basin_avg/TS_'+self.fileout+'.pdf')

if __name__ == "__main__":
    #for basin in ['Antarctic','Arctic','Nanse','Canadian']:
    for basin in ['Canadian','Nansen']:
        if basin in ['Antarctic']:
            prset = Products([UoR,GloSea5,MOVEG2i,GECCO2,EN4,\
                              ECDA,ORAP5,GLORYS2V4,CGLORS,SODA331],basin)
        else:
            prset = Products([UoR,GloSea5,MOVEG2i,GECCO2,EN4,\
                              ECDA,ORAP5,GLORYS2V4,CGLORS,SODA331,TOPAZ],basin)
        #prset = Products([GECCO2],basin)
        prset.readProfiles()
        for vname in ['T','S']:
            prset.getMultiModelMean(vname)
            prset.plotDepthProfile(vname)
        prset.plotTSProfile()
    print "Finnished!"
