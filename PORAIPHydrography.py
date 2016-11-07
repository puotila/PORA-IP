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
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime
from netcdftime import utime
from seawater import dens0

class ProfVar(object):
    """ Temperature (T) or salinity (S) vertical
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
        self.mz = np.average(level_bounds,axis=1)

class Product(object):
    """ ORAIP annual means of T and S at (plon,plat).
    """
    def __init__(self,plon,plat,\
                 syr=1993,eyr=2010,\
                 path='/home/uotilap/tiede/ORA-IP/annual_mean/',\
                 level_bounds=None):
        self.path = path
        self.plon, self.plat = plon, plat
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

    def readProfile(self,varname):
        """ varname is either T or S
        """
        pdata = getattr(self,varname)
        for li, lb in enumerate(self.LevelBounds[varname]):
            fn    = self.getNetCDFfilename(varname,0,lb[1])
            ldata = self.readOneFile(fn,self.ncvarname[varname])
            if lb[0]==0.:
                udata = 0.0*ldata
            else:
                fn    = self.getNetCDFfilename(varname,0,lb[0])
                udata = self.readOneFile(fn,self.ncvarname[varname])
            pdata.data[li] = (ldata - udata)/(lb[1] - lb[0])

    def readLatLon(self,fp):
        lat = np.array(fp.variables[self.nclatname][:])
        lon = np.array(fp.variables[self.nclonname][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        return lon, lat

    def getDates(self,fp):
        time = fp.variables[self.nctimename]
        if hasattr(time,'calendar'):
            cdftime = utime(time.units,calendar=time.calendar)
        else:
            cdftime = utime(time.units)
        return [cdftime.num2date(t) for t in time[:]]

    def findClosestLocation(self,lon,lat):
        ix = np.where(np.abs(lon-self.plon)==np.min(np.abs(lon-self.plon)))[0][0]
        iy = np.where(np.abs(lat-self.plat)==np.min(np.abs(lat-self.plat)))[0][0]
        return ix, iy

    def readOneFile(self,fn,ncvarname):
        """
        Read data from a netCDF file and return its temporal mean
        within given year range [syr, eyr] at the closest location to
        (plon, plat).
        """
        fp = self.getNetCDFfilepointer(fn)
        lon, lat = self.readLatLon(fp)
        ix, iy   = self.findClosestLocation(lon,lat)
        dates    = self.getDates(fp)
        ldata = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                ncvar = fp.variables[ncvarname]
                if ncvar.ndim==4:
                    ldata.append(np.ma.squeeze(np.ma.array(ncvar[i,:,iy,ix])))
                else:
                    ldata.append(np.ma.array(ncvar[i,iy,ix]))
        fp.close()
        return np.ma.mean(np.ma.array(ldata))

    def getLayeredDepthProfile(self,varname,depth,data):
        """
        Average 3D hires profile data according to level_bounds
        """
        ldata    = []
        for li, lb in enumerate(self.LevelBounds[varname]):
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            ldata.append(np.ma.average(data[iz]))
        return np.ma.array(ldata,mask=np.isnan(ldata))

class CGLORS(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( CGLORS, self).__init__(plon,plat,syr,eyr)
        self.dset = 'CGLORS'
        self.dsyr, self.deyr = 1989, 2014
        self.fpat = 'CGLORS025v5/%s_ORCA025-%s_%d.nc'
        self.ncvarname = {'T':'votemper',\
                          'S':'vosaline'}
        self.linecolor = self.scattercolor = 'lightgreen'

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

    def readProfile(self,varname):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        at the closest location to (plon, plat).
        """
        pdata = getattr(self,varname)
        ncvarname = self.ncvarname[varname]
        tdata = []
        for year in range(self.syr,self.eyr+1):
            fn, ncdepthname = self.getNetCDFfilename(varname,year)
            fp = self.getNetCDFfilepointer(fn)
            lon, lat = self.readLatLon(fp)
            dates    = self.getDates(fp,year)
            depth    = np.array(fp.variables[ncdepthname])
            ix, iy   = self.findClosestLocation(lon,lat)
            ncvar    = fp.variables[ncvarname]
            for i,date in enumerate(dates[:]):
                if date.year in range(self.syr,self.eyr+1):
                    data = np.ma.array(fp.variables[ncvarname][i,:,iy,ix])
                    tdata.append(self.getLayeredDepthProfile(varname,depth,data))
            fp.close()
        pdata.data = np.ma.average(tdata,axis=0)

class ECDA(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( ECDA, self).__init__(plon,plat,syr,eyr)
        self.dset = 'ECDA'
        self.dsyr, self.deyr = 1993, 2011
        self.fpat = "ECDA_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'yellow'
        self.lettercolor = 'black'

    def getDates(self,fp):
        if fp.variables.has_key(self.nctimename):
            time = fp.variables[self.nctimename]
        else:
            time = fp.variables[self.nctimename.upper()]
        if hasattr(time,'calendar'):
            cdftime = utime(time.units,calendar=time.calendar)
        else:
            cdftime = utime(time.units)
        return [cdftime.num2date(t) for t in time[:]]

class GloSea5(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( GloSea5, self).__init__(plon,plat,syr,eyr)
        self.dset = 'GloSea5'
        self.dsyr, self.deyr = 1993, 2014
        self.fpat = "GloSea5_GO5_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'blue'
        self.nctimename = 'time_counter'

    def readLatLon(self,fp):
        lat = np.arange(-89.5,90.5,1.)
        lon = np.arange(0.5,360.5,1.)
        return lon, lat

class MOVEG2(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( MOVEG2, self).__init__(plon,plat,syr,eyr)
        self.dset = 'MOVEG2'
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
    def __init__(self,plon,plat,syr,eyr):
        super( UoR, self).__init__(plon,plat,syr,eyr)
        self.dset = 'UoR'
        self.dsyr, self.deyr = 1989, 2010
        self.fpat = "UoR_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'vertically_integrated_salinity'}
        self.linecolor = self.scattercolor = 'lightblue'
        self.lettercolor = 'black'

class EN4(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( EN4, self).__init__(plon,plat,syr,eyr)
        self.dset  = 'EN4'
        self.dsyr, self.deyr = 1950, 2015
        self.fpat = "EN4.2.0.g10_int%s_annmean_%dto%d_%d-%dm.nc"
        self.ncvarname = {'T':'t_int_',\
                          'S':'s_int_'}
        self.linecolor = self.scattercolor = 'pink'

    def readProfile(self,varname):
        """ varname is either T or S
        """
        pdata = getattr(self,varname)
        for li, lb in enumerate(self.LevelBounds[varname]):
            ncvarname = "%s%d" % (self.ncvarname[varname],lb[1])
            fn    = self.getNetCDFfilename(varname,0,lb[1])
            ldata = self.readOneFile(fn,ncvarname)
            if lb[0]==0.:
                udata = 0.0*ldata
            else:
                ncvarname = "%s%d" % (self.ncvarname[varname],lb[0])
                fn    = self.getNetCDFfilename(varname,0,lb[0])
                udata = self.readOneFile(fn,ncvarname)
            pdata.data[li] = (ldata - udata)/(lb[1] - lb[0])

class GECCO2(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( GECCO2, self).__init__(plon,plat,syr,eyr)
        self.dset  = 'GECCO2'
        self.dsyr, self.deyr = 1948, 2012
        self.fpat = "GECCO2_int%s_annmean_%dto%d_%d-%dm_r360x180.nc"
        self.ncvarname = {'T':'vertically_integrated_temperature',\
                          'S':'S_0_%d'}
        self.linecolor = self.scattercolor = 'darkred'

    def getGECCO2SalinityDates(self,fp):
        return [datetime(year,1,1) for year in range(self.dsyr,self.deyr+1)]

    def readGECCO2SalinityProfile(self,varname='S'):
        fn = 'GECCO2_intS_annmean_1948to2011_all_layers_r360x180.nc'
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        dates     = self.getGECCO2SalinityDates(fp)
        ix, iy    = self.findClosestLocation(lon,lat)
        pdata     = getattr(self,varname)
        tdata     = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                ldata = []
                for li, lb in enumerate(self.LevelBounds[varname]):
                    ncvarname = self.ncvarname[varname] % lb[1]
                    lldata = np.ma.array(fp.variables[ncvarname][i,iy,ix])
                    if lb[0]==0.:
                        ludata = 0.0*lldata
                    else:
                        ncvarname = self.ncvarname[varname] % lb[0]
                        ludata = np.ma.array(fp.variables[ncvarname][i,iy,ix])
                    ldata.append((lldata*lb[1] - ludata*lb[0])/(lb[1] - lb[0]))
                tdata.append(ldata)
        fp.close()
        pdata.data = np.ma.average(tdata,axis=0)

class GLORYS2V4(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( GLORYS2V4, self).__init__(plon,plat,syr,eyr)
        self.dset  = 'GLORYS'
        self.fpat  = 'GSOP_GLORYS2V4_ORCA025_%s.nc'
        self.ncvarname = {'T':'z%dheatc',\
                          'S':'z%dsaltc'}
        self.ncdepthname = 'zdepth'
        self.linecolor = self.scattercolor = 'orange'

    def getNetCDFfilename(self,varname):
        if varname=='T':
            fn = self.fpat % ('HC')
        else:
            fn = self.fpat % ('SC')
        return fn

    def readProfile(self,varname):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        at the closest location to (plon, plat).
        """
        fn = self.getNetCDFfilename(varname)
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        dates     = self.getDates(fp)
        depth     = np.array(fp.variables[self.ncdepthname])
        ix, iy    = self.findClosestLocation(lon,lat)
        pdata     = getattr(self,varname)
        tdata     = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                ldata = []
                for li, lb in enumerate(self.LevelBounds[varname]):
                    ncvarname = self.ncvarname[varname] % lb[1]
                    lldata = np.ma.array(fp.variables[ncvarname][i,iy,ix])
                    if lb[0]==0.:
                        ludata = 0.0*lldata
                    else:
                        ncvarname = self.ncvarname[varname] % lb[0]
                        ludata = np.ma.array(fp.variables[ncvarname][i,iy,ix])
                    ldata.append((lldata - ludata)/(lb[1] - lb[0]))
                tdata.append(ldata)
        fp.close()
        pdata.data = np.ma.average(tdata,axis=0)

class TOPAZ(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( TOPAZ, self).__init__(plon,plat,syr,eyr)
        self.dset  = 'TOPAZ'
        self.dsyr, self.deyr = 1993, 2013
        self.fpat  = "TP4_r360x180_%s_%04d_%02d.nc"
        self.ncvarname = {'T':'temperature',\
                          'S':'salinity'}
        self.nclatname = 'latitude'
        self.nclonname = 'longitude'
        self.linecolor = self.scattercolor = 'green'

    def getNetCDFfilename(self,varname,year,month):
        if varname=='T':
            fn = self.fpat % ('temp',year,month)
        else:
            fn = self.fpat % ('salt',year,month)
        return fn

    def readProfile(self,varname):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        at the closest location to (plon, plat).
        """
        pdata = getattr(self,varname)
        ncvarname = self.ncvarname[varname]
        tdata = []
        for year in range(self.syr,self.eyr+1):
            for month in range(1,13):
                fn = self.getNetCDFfilename(varname,year,month)
                fp = self.getNetCDFfilepointer(fn)
                lon, lat = self.readLatLon(fp)
                depth    = np.array(fp.variables[self.ncdepthname])
                ix, iy   = self.findClosestLocation(lon,lat)
                ncvar    = fp.variables[ncvarname]
                data     = np.ma.array(fp.variables[ncvarname][:,iy,ix])
                fp.close()
                tdata.append(self.getLayeredDepthProfile(varname,depth,data))
        pdata.data = np.ma.average(tdata,axis=0)

class MultiModelMean(Product):
    def __init__(self,plon,plat):
        super( MultiModelMean, self).__init__(plon,plat)
        self.dset = 'MMM'
        self.linecolor = self.scattercolor = self.edgecolor = 'lightgrey'
        self.linestyle = ':'

    def calcMultiModelMean(self,products,vname):
        setattr(getattr(self,vname),'data',\
                np.ma.average([getattr(getattr(p,vname),'data')\
                for p in products],axis=0))

class Sumata(Product):
    def __init__(self,plon,plat,syr=1980,eyr=2015):
        super( Sumata, self).__init__(plon,plat,syr,eyr)
        self.dset = 'Sumata'
        self.dsyr, self.deyr = syr, eyr
        self.fpat = "ts-clim/hiroshis-clim/archive_v12_QC2_3_DPL_checked_2d_season_all-remapbil-oraip.nc"
        self.ncvarname = {'T':'temperature',\
                          'S':'salinity'}
        self.linecolor = self.scattercolor = 'black'

    def getNetCDFfilename(self):
        return self.fpat

    def readProfile(self,varname):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        at the closest location to (plon, plat).
        """
        fn = self.getNetCDFfilename()
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        depth     = np.array(fp.variables[self.ncdepthname])
        ix, iy    = self.findClosestLocation(lon,lat)
        ncvarname = self.ncvarname[varname]
        ncvar     = fp.variables[ncvarname]
        # if seasonal average is not representative we may
        # need to plot seasons separately
        data     = np.ma.average(ncvar[:,:,iy,ix],axis=0)
        fp.close()
        pdata = getattr(self,varname)
        pdata.data = self.getLayeredDepthProfile(varname,depth,data)

class WOA13(Sumata):
    def __init__(self,plon,plat,syr=1995,eyr=2012):
        super( WOA13, self).__init__(plon,plat,syr,eyr)
        self.dset = 'WOA13'
        self.fpat = "ts-clim/woa13/woa13-clim-1995-2012-season.nc"
        self.scattercolor = 'white'
        self.linestyle = '--'
        self.lettercolor = 'black'

class ORAP5(Product):
    def __init__(self,plon,plat,syr,eyr):
        super( ORAP5, self).__init__(plon,plat,syr,eyr)
        self.dset  = 'ORAP5'
        self.dsyr, self.deyr = 1993, 2012
        self.fpat  = "%s3D_orap5_1m_%d-%d_r360x180.nc"
        self.ncvarname = {'T':'votemper',\
                          'S':'vosaline'}
        self.ncdepthname = 'deptht'
        self.nctimename = 'time_counter'
        self.linecolor = self.scattercolor = 'red'

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

    def readProfile(self,varname):
        """ varname is either T or S
        Read data from a netCDF file and return its temporal mean
        at the closest location to (plon, plat).
        """
        fn = self.getNetCDFfilename(varname)
        fp = self.getNetCDFfilepointer(fn)
        lon, lat  = self.readLatLon(fp)
        dates     = self.getDates(fp)
        depth     = np.array(fp.variables[self.ncdepthname])
        ix, iy    = self.findClosestLocation(lon,lat)
        ncvarname = self.ncvarname[varname]
        ncvar     = fp.variables[ncvarname]
        pdata     = getattr(self,varname)
        tdata     = []
        for i,date in enumerate(dates[:]):
            if date.year in range(self.syr,self.eyr+1):
                data = np.ma.array(fp.variables[ncvarname][i,:,iy,ix])
                tdata.append(self.getLayeredDepthProfile(varname,depth,data))
        fp.close()
        pdata.data = np.ma.average(tdata,axis=0)

class Products(object):
    """ Container for ORA-IP products
    """
    def __init__(self,productobjs,lon,lat,\
                 syr=1993,eyr=2010):
        self.products = []
        self.lon, self.lat = lon, lat
        self.syr, self.eyr = syr, eyr
        self.title = "%4.1f$^\circ$N, %5.1f$^\circ$E" % \
                     (lat,lon)
        for pobj in productobjs:
            self.products.append(pobj(lon,lat,syr,eyr))
        self.mmm = MultiModelMean(lon,lat)
        self.sumata = Sumata(lon,lat)
        self.woa13 = WOA13(lon,lat)
        self.xlabel = {'T':"Temperature [$^\circ$C]",\
                       'S':"Salinity [psu]"}
        self.ylabel = "depth [m]"
        self.ymin   = 3100 # m
        # Products per 3 panels:
        self.ProductPanels = {'CGLORS':0,'GECCO2':0,'GLORYS':0,'GloSea5':0,\
                              'ORAP5':1,'TOPAZ':1,'UoR':1,\
                              'ECDA':2,'MOVEG2':2,'EN4':2}
        self.pretitle = ['(a)','(b)','(c)']
        modstr = '_'.join([p.dset for p in self.products])
        self.fileout = "%s_%04d-%04d_%04.1fN_%05.1fE" % \
                       (modstr,syr,eyr,lat,lon)

    def readProfiles(self,varname):
        for product in self.products:
            if varname=='S' and product.dset=='GECCO2':
                product.readGECCO2SalinityProfile()
            else:
                product.readProfile(varname)
        self.sumata.readProfile(varname)
        self.woa13.readProfile(varname)

    def getMultiModelMean(self,vname):
        self.mmm.calcMultiModelMean(self.products,vname)

    def getDataRange(self,vname):
        data = np.ma.hstack([np.ma.array(getattr(getattr(p,vname),'data')) \
            for p in self.products+[self.sumata]+[self.woa13]])
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

    def plotDepthProfile(self,vname):
        xmin, xmax = self.getDataRange(vname)
        fig = plt.figure(figsize=(8*2,10))
        axs = [plt.axes([0.10, 0.1, .2, .8]),\
               plt.axes([0.40, 0.1, .2, .8]),\
               plt.axes([0.70, 0.1, .2, .8])]
        lnes = [[],[],[]]
        lgds = [[],[],[]]
        #if self.lat>60:
        #    omproducts = [self.sumata,self.woa13,self.mmm]
        #else:
        #    omproducts = [self.woa13,self.mmm]
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            # First Sumata and WOA13 climatologies, and MMM
            for product in [self.sumata,self.woa13,self.mmm]:
                x = getattr(getattr(product,vname),'data')
                if x.all() is np.ma.masked:
                    """ all values are masked
                    """
                    continue
                lne.append(self.plotOneProfile(product,vname,ax))
                lgd.append(product.dset)
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
            lgd.append(product.dset)
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            ax.invert_yaxis()
            ax.set_yticks([0,100,300,700,1500,3000])
            ax.set_ylim(self.ymin,0)
            if vname=='T':
                ax.set_xlim(xmin-np.abs(0.1*xmin),xmax+np.abs(0.1*xmax))
            else:
                ax.set_xlim(xmin-np.abs(0.01*xmin),xmax+np.abs(0.01*xmax))
            ax.set_ylabel(self.ylabel)
            ax.set_title("%s %s" % (self.pretitle[panelno],self.title))
            ax.set_xlabel(self.xlabel[vname])
            if vname=='T':
                leg = ax.legend(lne,tuple(lgd),ncol=1,bbox_to_anchor=(1.2, 0.5))
            else:
                leg = ax.legend(lne,tuple(lgd),ncol=1,bbox_to_anchor=(0.6, 0.5))
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(1.0)
            leg.get_frame().set_alpha(1.0)
        #plt.show()
        plt.savefig(vname+'_'+self.fileout+'.pdf')

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
            # First climatologies (Sumata and WOA13) then MMM
            for product in [self.sumata,self.woa13,self.mmm]:
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
                lgd.append(product.dset)
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
            lgd.append(product.dset)
            # mark depth levels with numbers
        for panelno, ax in enumerate(axs):
            lne, lgd = lnes[panelno],lgds[panelno]
            ax.set_ylim(np.min(ti)-0.2,np.max(ti)+0.2)
            ax.set_xlim(np.min(si)-0.2,np.max(si)+0.2)
            ax.set_ylabel(self.xlabel['T'])
            ax.set_title("%s %s" % (self.pretitle[panelno],self.title))
            ax.set_xlabel(self.xlabel['S'])
            ax.legend(lne,tuple(lgd),ncol=1,\
                      bbox_to_anchor=(0.5, 0.9))
        #plt.show()
        plt.savefig('TS_'+self.fileout+'.pdf')

if __name__ == "__main__":
    # Central Arctic close to the North Pole
    lon, lat = 10., 88.
    # Nansen Basin that Marika says is good
    lon, lat = 100., 84.
    # Amerasian Basin that Marika says is good
    lon, lat = 215., 80.
    #prset = Products([UoR],lon,lat)
    prset = Products([UoR,GloSea5,MOVEG2,GECCO2,EN4,\
                      ECDA,ORAP5,TOPAZ,GLORYS2V4,CGLORS],lon,lat)
    for vname in ['T','S']:
        prset.readProfiles(vname)
        prset.getMultiModelMean(vname)
        prset.plotDepthProfile(vname)
    prset.plotTSProfile()
    print "Finnished!"
