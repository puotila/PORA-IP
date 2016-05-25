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
from netcdftime import utime

# Global variables
Alphabets = list(string.ascii_lowercase)
LevelBounds = np.array([[   0, 100],[ 100,  300],\
                        [ 300, 700],[ 700, 1500],\
                        [1500,3000],[3000, 6000]])
ModelLineColors = {"GECCO":"blue",\
                   "MOVECORE":"red",\
                   "CFSR":"green",\
                   "GLORYS2V1":"orange",\
                   "MOVEG2":"cyan",\
                   "CGLORS":"violet",\
                   "GloSea":"pink",\
                   "ORAS4":"lightgreen",\
                   "SODA":"brown",\
                   "ECDA":"yellow",\
                   "UoR":"lightblue",\
                   "EN3":"indigo",\
                  }
# Models without Arctic: "ECCOJPL","GODAS","MOVEC","K7"

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
        self.level_bounds = LevelBounds
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
        self.depth = np.array(fp.variables['depth'][:])
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
            iz = np.where((self.depth>=lb[0])&(self.depth<lb[1]))
            self.data.append(tavg_data[iz].mean())

class ORAIPprofile(object):
    """ ORAIP annual means in 1 deg grid.
    """
    def __init__(self,vname,plon,plat,dset='oraip',\
                 syr=1995,eyr=2004,\
                 path='/lustre/tmp/uotilap/ORA-IP/annual_mean/'):
        print "Reading %s" % dset
        self.dset = dset
        self.vname = vname # T or S
        self.path = path
        self.lat = plat
        self.lon = plon
        self.syr = syr
        self.eyr = eyr
        self.level_bounds = LevelBounds
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
        self.data = np.ma.array([np.ma.mean(d) for d in data])

    def readOneFile(self,fn,lb):
        fpat = ".+/%s_int%s_annmean_(\d+)to(\d+)_%d-%dm_r360x180.nc" % \
               (self.dset,self.vname,lb[0],lb[1])
        m = re.match(fpat,fn)
        dsyr, deyr = [int(i) for i in m.groups()]
        fp = nc.Dataset(fn)
        print "Reading %s." % fn
        lat = np.array(fp.variables['lat'][:])
        iy = np.where(np.abs(lat-self.lat)==np.min(np.abs(lat-self.lat)))[0][0]
        lon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        lon[np.where(lon<0.)] += 360.
        ix = np.where(np.abs(lon-self.lon)==np.min(np.abs(lon-self.lon)))[0][0]
        if self.dset in ['ECDA']:
            time = fp.variables['TIME']
        else:
            time = fp.variables['time']
        cdftime = utime(time.units,calendar=time.calendar)
        for i,t in enumerate(time[:]):
            date = cdftime.num2date(t)
            if date.year in range(self.syr,self.eyr+1):
                ldata = np.ma.array(fp.variables[self.ncname][i,iy,ix])
        fp.close()
        return ldata

class Experiments(object):
    """ Container for data (keep WOA13 first)
    """
    def __init__(self,exps):
        self.exps = exps
        if self.exps[0].vname=='T':
            self.xlabel = "Temperature [$^\circ$C]"
        else:
            self.xlabel = "Salinity [psu]"
        self.ylabel = "depth [m]"
        self.title = "%4.1f$^\circ$N, %5.1f$^\circ$E" % \
                     (exps[0].lat,exps[0].lon)
        modstr = '_'.join([e.dset for e in self.exps])
        self.figfile = "%s_%s_%04d-%04d_%04.1fN_%05.1fE.png" % \
                     (self.exps[0].vname,modstr,\
                      self.exps[0].syr,self.exps[0].eyr,\
                      exps[0].lat,exps[0].lon)
        self.depth = np.hstack((exps[0].level_bounds[:,0],4000))
        # MMM = multimodel mean
        exmmm = copy.copy(exps[1])
        exmmm.dset = 'MMM'
        exmmm.data = np.ma.average([e.data for e in exps[1:]],axis=0)
        self.exps.append(exmmm)

    def plotProfiles(self):
        fig = plt.figure(figsize=(8*2,10))
        ax1 = plt.axes([0.10, 0.1, .2, .8])
        ax2 = plt.axes([0.40, 0.1, .2, .8])
        ax3 = plt.axes([0.70, 0.1, .2, .8])
        lns_pplot = [4,4,4] # lines per plot
        exps = self.exps[1:-1]
        li = 0 # line index
        for ia, ax in enumerate([ax1,ax2,ax3]):
            lnes, lgnds = [],[]
            # WOA13 climatology
            y = np.ma.hstack((self.exps[0].data,self.exps[0].data[-1]))
            lnes.append(ax.plot(y,self.depth,lw=3,\
                                drawstyle='steps-pre',color='black')[0])
            lgnds.append(self.exps[0].dset)
            # multi-model mean
            y = np.ma.hstack((self.exps[-1].data,self.exps[-1].data[-1]))
            lnes.append(ax.plot(y,self.depth,lw=3,linestyle='--',\
                                drawstyle='steps-pre',color='darkgrey')[0])
            lgnds.append(self.exps[-1].dset)
            # then individual models
            for exp in exps[li:li+lns_pplot[ia]]:
                y = np.ma.hstack((exp.data,exp.data[-1]))
                lnes.append(ax.plot(y,self.depth,lw=2,\
                                    drawstyle='steps-pre',color=ModelLineColors[exp.dset])[0])
                lgnds.append(exp.dset)
            li += lns_pplot[ia]
            ax.invert_yaxis()
            ax.set_ylim(4000,0)
            ax.set_ylabel(self.ylabel)
            ax.set_title("%s) %s" % (Alphabets[ia],self.title))
            ax.set_xlabel(self.xlabel)
            ax.legend(lnes,tuple(lgnds),ncol=1,bbox_to_anchor=(1.2, 0.5))
        #plt.show()
        plt.savefig(self.figfile)

if __name__ == "__main__":
    # lons should be E and lats should be N
    #lon, lat = 7., 80.
    lon, lat = 220., 80.
    vname = 'T' # or 'S'
    models = ModelLineColors.keys()
    models.sort()
    experiments = Experiments([WOA13profile(vname,lon,lat)]+ \
                              [ORAIPprofile(vname,lon,lat,dset=model) \
                               for model in models])
    experiments.plotProfiles()
    print "Finnished!"
