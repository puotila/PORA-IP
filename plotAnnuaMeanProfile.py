#!/usr/bin/env python
"""
Plot ORA-IP annual mean profile (T or S)
and WOA13 1995-2004 profile.
"""

import os
import sys
sys.path.insert(1,'/lustre/tmp/uotilap/marnelam/netCDF4-1.2.2-py2.7-linux-x86_64.egg')
#sys.path.insert(1,'/opt/Python/2.7/lib/python2.7/site-packages')
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

# Global variables
Alphabets = list(string.ascii_lowercase)
LevelBounds = {}
LevelBounds['T'] = np.array([[   0, 100],[ 100,  300],\
                             [ 300, 700],[ 700, 1500],\
                             [1500,3000],[3000, 6000]])
LevelBounds['S'] = np.array([[   0,  10],[  10,  100],\
                             [ 100, 300],[ 300,  700],\
                             [ 700,1500],[1500, 3000],\
                                         [3000, 4000]])
ModelLineColors = {"GECCO2":"darkred",\
                  # "GMAO":"blue",\
                  #"MOVECORE":"red",\
                  #"ORAS3":"red",\
                   "ORAP5":"red",\
                  #"CFSR":"green",\
                   "TP4":"green",\
                   "GLORYS2V1":"orange",\
                   "G2V3":"orange",\
                   "GSOP_GLORYS2V4":"orange",\
                   "MOVEG2":"cyan",\
                  #"PEODAS":"cyan",\
                   "CGLORS":"lightgreen",\
                   "GloSea":"blue",\
                   "GloSea5_GO5":"blue",\
                  #"K7ODA":"pink",\
                  #"ORAS4":"lightgreen",\
                   "SODA":"brown",\
                   "ECDA":"yellow",\
                  #"PECDAS":"yellow",\
                   "UoR":"lightblue",\
                   "UR025":"lightblue",\
                   "EN3":"pink",\
                   "EN3v2a":"pink",\
                  }
# Models without Arctic: "ECCOJPL","GODAS","MOVEC","K7"
# Models without sea ice: SODA
# Not a model: EN3

class WOA13profile(object):
    """ WOA13 decadal means in 1 deg grid.
    """
    def __init__(self,vname,plon,plat,dset='WOA13',\
                 months=range(1,13),syr=1995,eyr=2004,\
                 path='/lustre/tmp/uotilap/marnelam/WOA13/'):
        self.dset = dset
        self.vname = vname # T or S
        self.lon = plon
        self.lat = plat
        self.months = months
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        if vname=='T':
            ncname = 't_an'
        elif vname=='S':
            ncname = 's_an'
        else:
            print "%s has not been implemented! Exiting..." % vname
            sys.exit(1)
        data = []
        fn = os.path.join(path,"%s_merged_%s_%04d-%04d.nc" % \
                               (dset.lower(),vname.lower(),syr,eyr))
        fp = nc.Dataset(fn)
        depth = np.array(fp.variables['depth'][:])
        lat = np.array(fp.variables['lat'][:])
        iy = np.where(np.abs(lat-plat)==np.min(np.abs(lat-plat)))[0][0]
        lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        ix = np.where(np.abs(lon-plon)==np.min(np.abs(lon-plon)))[0][0]
        time = np.array(fp.variables['time'][:])+1
        for t in time:
            if t in months:
                data.append(fp.variables[ncname][t-1,:,iy,ix])
        fp.close()
        tavg_data = np.ma.average(data,axis=0)
        # vertical layer averaging
        self.data = []
        for lb in self.level_bounds:
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            self.data.append(tavg_data[iz].mean())
        self.depth = np.hstack((self.level_bounds[:,0],4000))

class ORAIPprofile(object):
    """ ORAIP annual means in 1 deg grid.
    """
    def __init__(self,vname,plon,plat,dset='oraip',\
                 syr=1993,eyr=2009,\
                 path='/lustre/tmp/uotilap/ORA-IP/annual_mean/'):
                 #path='/lustre/tmp/marnelam/reanalyses/heat_content/'):
        print "Reading %s" % dset
        self.dset = dset
        if dset=='EN3' and vname=='S':
            self.dset = dset = 'EN3v2a'
        self.vname = vname # T or S
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        self.depth = np.hstack((self.level_bounds[:,0],4000))
        if self.dset in ['GECCO2'] and vname=='S':
            #self.data = np.ma.masked_less_equal(self.readGECCO2salinity(),32)
            self.data = self.readGECCO2salinity()
            return
        if vname=='T':
            self.ncname = 'vertically_integrated_temperature'
        elif vname=='S':
            self.ncname = 'vertically_integrated_salinity'
        else:
            print "%s has not been implemented! Exiting..." % vname
            sys.exit(1)
        data = [[] for i in self.level_bounds[:,0]]
        for li, lb in enumerate(self.level_bounds):
            fns = glob.glob(os.path.join(path,"%s_int%s_annmean_????to????_%d-%dm_r360x180.nc" % \
                                              (dset,vname,lb[0],lb[1])))
            if len(fns) and os.path.exists(fns[0]):
                    ldata = self.readOneFile(fns[0],lb)
            else:
                # level is missing, need to calculate it from two other
                # levels. E.g. 100-300m is 0-300m minus 0-100m
                fn = glob.glob(os.path.join(path,"%s_int%s_annmean_????to????_%d-%dm_r360x180.nc" % \
                                                 (dset,vname,0,lb[0])))[0]
                ldatau = self.readOneFile(fn,[0,lb[0]])
                if dset in ['GloSea5_GO5'] and lb[1] in [6000]:
                    fn = glob.glob(os.path.join(path,"%s_int%s_annmean_????to????_%d-%s_r360x180.nc" % \
                                                     (dset,vname,0,'bottom')))[0]
                else:
                    fn = glob.glob(os.path.join(path,"%s_int%s_annmean_????to????_%d-%dm_r360x180.nc" % \
                                                     (dset,vname,0,lb[1])))[0]
                ldatal = self.readOneFile(fn,[0,lb[1]])
                ldata = ldatal - ldatau
            data[li].append(ldata/np.diff(self.level_bounds)[li])
        if vname=='T':
            self.data = np.ma.masked_equal(np.ma.squeeze(data),0)
        else: #S
            #self.data = np.ma.masked_less_equal(np.ma.squeeze(data),32)
            self.data = np.ma.squeeze(data)
        #self.data = np.ma.squeeze(data)

    def readGECCO2salinity(self,dsyr=1948,deyr=2011):
        """
        GECCO2_intS_annmean_1948to2011_all_layers_r360x180.nc
        """
        fn = "%s_int%s_annmean_%04dto%04d_all_layers_r360x180.nc" % \
               (self.dset,self.vname,dsyr,deyr)
        fp = nc.Dataset(os.path.join(self.path,fn))
        print "Reading %s." % fn
        lat = np.array(fp.variables['lat'][:])
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        years = np.arange(dsyr,deyr+1)
        time = fp.variables['time']
        # as 0-10m is missing, assume it is the same than 0-100m
        # values are level means NOT level integrals!
        # also assuming that S_K in the netcdf file is the mean salinity from
        # surface to bottom
        ncnames = ['S_0_100','S_0_100','S_0_300','S_0_700','S_0_1500','S_0_3000','S_0_4000']
        # level thicknesses
        zdpths = [100,100,300,700,1500,3000,4000]
        data, zdata = [], []
        for i,t in enumerate(time[:]):
            year = years[i]
            if year in range(self.syr,self.eyr+1):
                zdata.append(np.ma.array(fp.variables[ncnames[0]][i,iy,ix]))
        data.append(np.ma.mean(zdata)) # 0-100m
        data.append(np.ma.mean(zdata)) # 0-100m
        for li,ncname in enumerate(ncnames[2:]):
            zdata = []
            for i,t in enumerate(time[:]):
                year = years[i]
                if year in range(self.syr,self.eyr+1):
                    zdatal = np.ma.array(fp.variables[ncname][i,iy,ix])*zdpths[li+2]
                    zdatau = np.ma.array(fp.variables[ncnames[li+1]][i,iy,ix])*zdpths[li+1]
                    zdata.append((zdatal - zdatau)/(zdpths[li+2]-zdpths[li+1]))
            data.append(np.ma.mean(zdata))
        fp.close()
        return np.ma.array(data)

    def readOneFile(self,fn,lb):
        #fpat = ".+/%s_int%s_annmean_(\d+)to(\d+)_%d-%dm_r360x180.nc" % \
        #       (self.dset,self.vname,lb[0],lb[1])
        fpat = ".+/%s_int%s_annmean_(\d+)to(\d+)_%d-.+_r360x180.nc" % \
               (self.dset,self.vname,lb[0])
        m = re.match(fpat,fn)
        dsyr, deyr = [int(i) for i in m.groups()]
        fp = nc.Dataset(fn)
        print "Reading %s." % fn
        if fp.variables.has_key('latitude'):
            lat = np.array(fp.variables['latitude'][:])
        elif fp.variables.has_key('lat'):
            lat = np.array(fp.variables['lat'][:])
        else:
            lat = np.arange(-89.5,90.)
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        if fp.variables.has_key('longitude'):
            lon = np.array(fp.variables['longitude'][:])
        elif fp.variables.has_key('lon'):
            lon = np.array(fp.variables['lon'][:])
        else:
            lon = np.arange(.5,360.)
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        if fp.variables.has_key('TIME'):
            time = fp.variables['TIME']
        elif fp.variables.has_key('time_counter'):
            time = fp.variables['time_counter']
        else:
            time = fp.variables['time']
        m1 = re.search('months since\s+(\d+)-(\d+)-(\d+)',time.units)
        m2 = re.search('month since\s+(\d+)-(\d+)-(\d+)',time.units)
        if m1 or m2:
            if m1:
                m = m1
            else:
                m = m2
            year0, month0, day0 = [int(s) for s in m.groups()]
            dates = [datetime(year0+int(t/12),month0+int(t%12),day0) for t in time[:]]
        else:
            if hasattr(time,'calendar'):
                cdftime = utime(time.units,calendar=time.calendar)
            else:
                cdftime = utime(time.units)
            dates = [cdftime.num2date(t) for t in time[:]]
        ldata = []
        for i,t in enumerate(time[:]):
            date = dates[i]
            if date.year in range(self.syr,self.eyr+1):
                if self.dset in ['K7ODA'] or \
                  (self.dset in ['MOVEG2'] and self.vname=='S'):
                    ldata.append(np.ma.array(fp.variables[self.ncname][i,:,iy,ix]).squeeze())
                else:
                    ldata.append(np.ma.array(fp.variables[self.ncname][i,iy,ix]))
        fp.close()
        return np.ma.mean(np.ma.array(ldata))

class TOPAZprofile(object):
    """ TP4 annual means in 1 deg grid.
        Monthly files of form: TP4_r360x180_temp|salt_1999_12.nc
        temperature|salinity(depth, latitude, longitude)
    """
    def __init__(self,vname,plon,plat,dset='TP4',\
                 syr=1993,eyr=2009,\
                 path='/lustre/tmp/uotilap/ORA-IP/annual_mean/'):
        print "Reading %s" % dset
        self.dset = dset
        self.vname = vname # T or S
        if self.vname=='T':
            fvarstr = 'temp'
            self.ncname = 'temperature'
        else:
            fvarstr = 'salt'
            self.ncname = 'salinity'
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        fn = "%s_r360x180_%s_%04d_%02d.nc" % (self.dset,fvarstr,syr,1)
        fp = nc.Dataset(fn)
        depth = np.array(fp.variables['depth'][:])
        lat = np.array(fp.variables['latitude'][:])
        lon = np.array(fp.variables['longitude'][:])
        fp.close()
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        self.months = range(1,13)
        years = range(syr,eyr+1)
        ldata = []
        for y in years:
            for m in self.months:
                fn = "%s_r360x180_%s_%04d_%02d.nc" % (self.dset,fvarstr,y,m)
                fp = nc.Dataset(fn)
                ldata.append(np.ma.array(fp.variables[self.ncname][:,iy,ix])) # depth, lat, lon
                fp.close()
        tavg_data = np.ma.mean(ldata,axis=0)
        self.data = []
        for lb in self.level_bounds:
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            self.data.append(tavg_data[iz].mean())
        self.depth = np.hstack((self.level_bounds[:,0],4000))

class ORAP5profile(object):
    """ ORAP5 annual means in 1 deg grid.
        temperature|salinity3D_orap5_1m_1993-2012_r360x180.nc
        votemper|vosaline(time_counter, deptht, lat, lon)
    """
    def __init__(self,vname,plon,plat,dset='ORAP5',\
                 syr=1993,eyr=2009,\
                 path='/lustre/tmp/uotilap/ORA-IP/annual_mean/'):
        print "Reading %s" % dset
        self.dset = dset
        self.vname = vname # T or S
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        if self.vname=='T':
            fvarstr = 'temperature'
            self.ncname = 'votemper'
        else:
            fvarstr = 'salinity'
            self.ncname = 'vosaline'
        fn = "%s3D_%s_1m_1993-2012_r360x180.nc" % (fvarstr,self.dset.lower())
        fp = nc.Dataset(fn)
        depth = np.array(fp.variables['deptht'][:])
        lat = np.array(fp.variables['lat'][:])
        lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        time = fp.variables['time_counter']
        m = re.search('months since\s+(\d+)-(\d+)-(\d+)',time.units)
        if m:
            year0, month0, day0 = [int(s) for s in m.groups()]
            dates = [datetime(year0+int(t/12),month0+int(t%12),day0) for t in time[:]]
        else:
            print "Cannot convert time to dates!"; sys.exit(1)
        ldata = []
        for i,t in enumerate(time[:]):
            date = dates[i]
            if date.year in range(self.syr,self.eyr+1):
                ldata.append(np.ma.array(fp.variables[self.ncname][i,:,iy,ix]).squeeze())
        fp.close()
        tavg_data = np.ma.mean(ldata,axis=0)
        self.data = []
        for lb in self.level_bounds:
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            self.data.append(tavg_data[iz].mean())
        self.depth = np.hstack((self.level_bounds[:,0],4000))

class GSOP_GLORYS2V4profile(object):
    """ GSOP_GLORYS2V4 annual means in 1 deg grid.
        GSOP_GLORYS2V3_ORCA025_H|SC.nc
        votemper|vosaline(time, lat, lon)
    """
    def __init__(self,vname,plon,plat,dset='GSOP_GLORYS2V4',\
                 syr=1993,eyr=2009,\
                 path='/lustre/tmp/uotilap/ORA-IP/annual_mean/'):
        print "Reading %s" % dset
        self.dset = dset
        self.vname = vname # T or S
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        self.depth = np.hstack((self.level_bounds[:,0],4000))
        if self.vname=='T':
            fn = "%s_ORCA025_HC.nc" % (self.dset)
            varname = 'heatc'
        else:
            fn = "%s_ORCA025_SC.nc" % (self.dset)
            varname = 'saltc'
        fp = nc.Dataset(fn)
        lat = np.array(fp.variables['lat'][:])
        lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        time = fp.variables['time']
        cdftime = utime(time.units,calendar=time.calendar)
        dates = [cdftime.num2date(t) for t in time[:]]
        data = [[] for i in self.level_bounds[:,0]]
        for li, lb in enumerate(self.level_bounds):
            ncnameu = "z%d%s" % (lb[0],varname)
            if lb[1]==6000:
                ncnamel = "zbot%s" % (varname)
            else:
                ncnamel = "z%d%s" % (lb[1],varname)
            ldata  = []
            for i,t in enumerate(time[:]):
                date = dates[i]
                if date.year in range(self.syr,self.eyr+1):
                    ldatal = np.ma.array(fp.variables[ncnamel][i,iy,ix])
                    if ncnameu in ['z0heatc','z0saltc']:
                        ldata = np.ma.mean(ldatal)
                    else:
                        ldatau = np.ma.array(fp.variables[ncnameu][i,iy,ix])
                        ldata = np.ma.mean(ldatal) - np.ma.mean(ldatau)
            data[li].append(ldata/np.diff(self.level_bounds)[li])
        fp.close()
        if vname=='T':
            self.data = np.ma.masked_equal(np.ma.squeeze(data),0)
        else: #S
            self.data = np.ma.masked_less_equal(np.ma.squeeze(data),32)

class Experiments(object):
    """ Container for data (keep WOA13 first)
    """
    def __init__(self,exps):
        self.exps = exps
        self.vname = self.exps[0].vname
        if self.vname=='T':
            self.xlabel = "Temperature [$^\circ$C]"
        else:
            self.xlabel = "Salinity [psu]"
        self.ylabel = "depth [m]"
        self.title = "%4.1f$^\circ$N, %5.1f$^\circ$E" % \
                     (exps[0].lat,exps[0].lon)
        modstr = '_'.join([e.dset for e in self.exps])
        self.figfile = "%s_%s_%04d-%04d_%04.1fN_%05.1fE.png" % \
                     (self.vname,modstr,\
                      self.exps[0].syr,self.exps[0].eyr,\
                      exps[0].lat,exps[0].lon)
        # MMM = multimodel mean
        exmmm = copy.copy(exps[1])
        exmmm.dset = 'MMM'
        exmmm.data = np.ma.average([e.data for e in exps[1:]],axis=0)
        self.exps.append(exmmm)
        # Products per 3 panels:
        self.ProductPanel = [{'T':['CGLORS','GECCO2','GSOP_GLORYS2V4','GloSea5_GO5'],\
                              'S':['CGLORS','GECCO2','GSOP_GLORYS2V4','GloSea5_GO5']},\
                             {'T':['ORAP5','TP4','UoR'],\
                              'S':['ORAP5','TP4','UoR']},\
                             {'T':['ECDA','MOVEG2','EN3'],\
                              'S':['ECDA','MOVEG2','EN3v2a']}]

    def plotProfiles(self):
        fig = plt.figure(figsize=(8*2,10))
        ax1 = plt.axes([0.10, 0.1, .2, .8])
        ax2 = plt.axes([0.40, 0.1, .2, .8])
        ax3 = plt.axes([0.70, 0.1, .2, .8])
        for ia, ax in enumerate([ax1,ax2,ax3]):
            lnes, lgnds = [],[]
            # WOA13 climatology
            y = np.ma.hstack((self.exps[0].data,self.exps[0].data[-1]))
            lnes.append(ax.plot(y,self.exps[0].depth,lw=3,\
                                drawstyle='steps-pre',color='black')[0])
            lgnds.append(self.exps[0].dset)
            # multi-model mean
            y = np.ma.hstack((self.exps[-1].data,self.exps[-1].data[-1]))
            lnes.append(ax.plot(y,self.exps[-1].depth,lw=3,linestyle='--',\
                                drawstyle='steps-pre',color='darkgrey')[0])
            lgnds.append(self.exps[-1].dset)
            # then individual models
            for ename in self.ProductPanel[ia][self.vname]:
                exp = [e for e in self.exps if e.dset==ename][0]
                y = np.ma.hstack((exp.data,exp.data[-1]))
                lnes.append(ax.plot(y,exp.depth,lw=2,\
                                    drawstyle='steps-pre',color=ModelLineColors[exp.dset])[0])
                if exp.dset=='GSOP_GLORYS2V4':
                    lgnds.append('GLORYS2V4')
                elif exp.dset=='GloSea5_GO5':
                    lgnds.append('GloSea5')
                elif exp.dset=='EN3v2a':
                    lgnds.append('EN3')
                else:
                    lgnds.append(exp.dset)
            ax.invert_yaxis()
            ax.set_ylim(4000,0)
            if self.vname=='S':
                ax.set_xlim(32,35)
            ax.set_ylabel(self.ylabel)
            ax.set_title("%s) %s" % (Alphabets[ia],self.title))
            ax.set_xlabel(self.xlabel)
            if self.vname=='T':
                ax.legend(lnes,tuple(lgnds),ncol=1,bbox_to_anchor=(1.2, 0.5))
            else:
                ax.legend(lnes,tuple(lgnds),ncol=1,bbox_to_anchor=(0.6, 0.5))
        #plt.show()
        plt.savefig(self.figfile)

    def plotTSdiagram(self,sxps):
        fig = plt.figure(figsize=(8*2.5,10))
        ax1 = plt.axes([0.10, 0.1, .2, .8])
        ax2 = plt.axes([0.40, 0.1, .2, .8])
        ax3 = plt.axes([0.70, 0.1, .2, .8])
        texps, sexps = self.exps, sxps.exps
        # Calculate how many gridcells we need in the x and y dimensions
        tmin = np.ma.min([e.data for e in texps])
        tmax = np.ma.max([e.data for e in texps])
        smin = np.ma.min([e.data for e in sexps])
        smax = np.ma.max([e.data for e in sexps])
        xdim = round((smax-smin)/0.1+1,0)
        ydim = round((tmax-tmin)/0.1+1,0)
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
        dens -= 1000
        for ia, ax in enumerate([ax1,ax2,ax3]):
            CS = ax.contour(si,ti,dens, linestyles='dashed', colors='k')
            ax.clabel(CS, fontsize=12, inline=1, fmt='%2.1f') # Label every second level
            pnts, lgnds = [],[]
            # WOA13 climatology
            t = np.ma.hstack((texps[0].data,texps[0].data[-1]))
            s = np.ma.hstack((sexps[0].data[1:],sexps[0].data[-1]))
            pnts.append(ax.scatter(s,t,lw=3,color='black'))
            lgnds.append(texps[0].dset)
            # multi-model mean
            t = np.ma.hstack((texps[-1].data,texps[-1].data[-1]))
            s = np.ma.hstack((sexps[-1].data[1:],sexps[-1].data[-1]))
            pnts.append(ax.scatter(s,t,lw=3,color='darkgrey'))
            lgnds.append(texps[-1].dset)
            # then individual models
            for ename in self.ProductPanel[ia][self.vname]:
                texp = [e for e in texps if e.dset==ename][0]
                if ename=='EN3':
                    sexp = [e for e in sexps if e.dset=='EN3v2a'][0]
                else:
                    sexp = [e for e in sexps if e.dset==ename][0]
                t = np.ma.hstack((texp.data,texp.data[-1]))
                s = np.ma.hstack((sexp.data[1:],sexp.data[-1]))
                pnts.append(ax.scatter(s,t,lw=3,color=ModelLineColors[texp.dset]))
                if texp.dset=='GSOP_GLORYS2V4':
                    lgnds.append('GLORYS2V4')
                elif texp.dset=='GloSea5_GO5':
                    lgnds.append('GloSea5')
                else:
                    lgnds.append(texp.dset)
            if ia==0:
                ax.set_ylabel("temperature ($^\circ$C)")
            ax.set_title("%s) %s" % (Alphabets[ia],self.title))
            ax.set_xlabel('salinity (ppm)')
            ax.legend(pnts,tuple(lgnds),ncol=1,bbox_to_anchor=(0.65, 1.0))
        #plt.show()
        plt.savefig("S"+self.figfile)

if __name__ == "__main__":
    # lons should be E and lats should be N
    lon, lat = 10., 88.
    #lon, lat = 100., 83.
    #lon, lat = 7., 80.
    #lon, lat = 220., 80.
    vname = 'S' # 'T' or 'S'
    #models = ['CGLORS','ECDA','GloSea5_GO5',\
    #          'MOVEG2','UoR','EN3v2a','GECCO2']
    models = ['CGLORS', 'ECDA','GloSea5_GO5',\
              'MOVEG2', 'UoR','EN3','GECCO2']
    Sxperiments = Experiments([WOA13profile(vname,lon,lat)]+ \
                              [ORAIPprofile(vname,lon,lat,dset=model) \
                               for model in models]+\
                              [GSOP_GLORYS2V4profile(vname,lon,lat),\
                               ORAP5profile(vname,lon,lat),\
                               TOPAZprofile(vname,lon,lat)])
    vname = 'T' # 'T' or 'S'
    Txperiments = Experiments([WOA13profile(vname,lon,lat)]+ \
                              [ORAIPprofile(vname,lon,lat,dset=model) \
                               for model in models]+\
                              [GSOP_GLORYS2V4profile(vname,lon,lat),\
                               ORAP5profile(vname,lon,lat),\
                               TOPAZprofile(vname,lon,lat)])
    for experiments in [Txperiments,Sxperiments]:
        experiments.plotProfiles()
    Txperiments.plotTSdiagram(Sxperiments)
    print "Finnished!"
