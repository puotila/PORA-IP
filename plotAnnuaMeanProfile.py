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
                   "MOVEG2":"cyan",\
                  #"PEODAS":"cyan",\
                   "CGLORS":"lightgreen",\
                   "GloSea":"blue",\
                   "GloSea5":"blue",\
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
        self.vname = vname # T or S
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds[vname]
        self.depth = np.hstack((self.level_bounds[:,0],4000))
        if self.dset in ['GECCO2'] and vname=='S':
            self.data = np.ma.masked_less_equal(self.readGECCO2salinity(),32)
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
                fn = glob.glob(os.path.join(path,"%s_int%s_annmean_????to????_%d-%dm_r360x180.nc" % \
                                                 (dset,vname,0,lb[1])))[0]
                ldatal = self.readOneFile(fn,[0,lb[1]])
                ldata = ldatal - ldatau
            data[li].append(ldata/np.diff(self.level_bounds)[li])
        if vname=='T':
            self.data = np.ma.masked_equal(np.ma.squeeze(data),0)
        else: #S
            self.data = np.ma.masked_less_equal(np.ma.squeeze(data),32)
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
        ncnames = ['S_0_10','S_0_100','S_0_300','S_0_700','S_0_1500','S_0_3000','S_0_4000']
        data = []
        for i,t in enumerate(time[:]):
            year = years[i]
            if year in range(self.syr,self.eyr+1):
                zdata = [0]
                zdata.append(np.ma.array(fp.variables[ncnames[1]][i,iy,ix])/self.level_bounds[1][1])
                for li,ncname in enumerate(ncnames[1:-1]):
                    udata = np.ma.array(fp.variables[ncnames[li+1]][i,iy,ix])
                    ldata = np.ma.array(fp.variables[ncnames[li+2]][i,iy,ix])
                    zdata.append((ldata-udata)/np.diff(self.level_bounds)[li+2])
                data.append(zdata)
        fp.close()
        return np.ma.mean(np.ma.array(data),axis=0)

    def readOneFile(self,fn,lb):
        fpat = ".+/%s_int%s_annmean_(\d+)to(\d+)_%d-%dm_r360x180.nc" % \
               (self.dset,self.vname,lb[0],lb[1])
        m = re.match(fpat,fn)
        dsyr, deyr = [int(i) for i in m.groups()]
        fp = nc.Dataset(fn)
        print "Reading %s." % fn
        if fp.variables.has_key('latitude'):
            lat = np.array(fp.variables['latitude'][:])
        else:
            lat = np.array(fp.variables['lat'][:])
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        if fp.variables.has_key('longitude'):
            lon = np.array(fp.variables['longitude'][:])
        else:
            lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        if fp.variables.has_key('TIME'):
            time = fp.variables['TIME']
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
        self.ProductPanel = [{'T':['CGLORS','GECCO2','GLORYS2V1','GloSea'],\
                              'S':['CGLORS','GECCO2','G2V3','GloSea5']},\
                             {'T':['ORAP5','TP4','UoR'],\
                              'S':['ORAP5','TP4','UoR']},\
                             {'T':['ECDA','MOVEG2','EN3'],\
                              'S':['ECDA','MOVEG2','EN3v2a']}]

    def plotProfiles(self):
        fig = plt.figure(figsize=(8*2,10))
        ax1 = plt.axes([0.10, 0.1, .2, .8])
        ax2 = plt.axes([0.40, 0.1, .2, .8])
        ax3 = plt.axes([0.70, 0.1, .2, .8])
        exps = self.exps[1:-1]
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

if __name__ == "__main__":
    # lons should be E and lats should be N
    lon, lat = 10., 88.
    #lon, lat = 100., 83.
    #lon, lat = 7., 80.
    #lon, lat = 220., 80.
    vname = 'S' # 'T' or 'S'
    for vname in ['T','S']:
        # add TP4, ORAP5, GECCO2
        if vname=='T':
            models = ['CGLORS', 'ECDA','GLORYS2V1', 'GloSea',\
                      'MOVEG2', 'UoR','EN3','GECCO2']
        else: #S
            models = ['CGLORS','ECDA','G2V3','GloSea5',\
                      'MOVEG2','UoR','EN3v2a','GECCO2']
        experiments = Experiments([WOA13profile(vname,lon,lat)]+ \
                                  [ORAIPprofile(vname,lon,lat,dset=model) \
                                   for model in models]+\
                                  [ORAP5profile(vname,lon,lat),TOPAZprofile(vname,lon,lat)])
        experiments.plotProfiles()
    print "Finnished!"
