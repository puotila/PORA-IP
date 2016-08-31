#!/usr/bin/env python
"""
Read Hiroshis (AWI) Arctic T, S climatology
based on recent cruises. Essentially replaces PHC3.
As it is sparse, where a value is missing it is replaced
by WOA13 value (which may be nonsense).
"""

import sys
sys.path.append('/lustre/tmp/uotilap/ORA-IP/annual_mean')
import numpy as np
import netCDF4 as nc
from plotAnnuaMeanProfile import LevelBounds

class Hiroshis(object):
    def __init__(self,fn='/lustre/tmp/uotilap/ORA-IP/annual_mean/ts-clim/hiroshis-clim/archive_v12_QC2_3_DPL_checked_2d_season_int-remapbil-oraip.nc'):
        self.dset, self.syr, self.eyr = 'Sumata', 1980, 2015
        fp = nc.Dataset(fn)
        self.olon = np.array(fp.variables['lon'][:])
        # transfer negative lons to positive
        self.olon[np.where(self.olon<0.)] += 360.
        self.olat = np.array(fp.variables['lat'][:])
        self.odepth = np.array(fp.variables['depth'][:])
        self.temperature = np.ma.array(fp.variables['temperature'][:])
        self.salinity    = np.ma.array(fp.variables['salinity'][:])
        fp.close()

    def getPoint(self,plon,plat,vname):
        """ get the closest point of (plon,plat)
        """
        self.vname, self.lon, self.lat = vname, plon, plat
        iy = np.where(np.abs(plat-self.olat)==np.min(np.abs(plat-self.olat)))[0][0]
        ix = np.where(np.abs(plon-self.olon)==np.min(np.abs(plon-self.olon)))[0][0]
        print "Looking for (%f,%f), closest at (%f,%f)" % \
              (plon,plat,self.olon[ix],self.olat[iy])
        depth = self.odepth
        if vname=='S':
            data, level_bounds = self.salinity, LevelBounds['S']
        else:
            data, level_bounds = self.temperature, LevelBounds['T']
        tavg_data = np.ma.average(data,axis=0)
        # vertical layer averaging
        self.data = []
        for lb in level_bounds:
            iz = np.where((depth>=lb[0])&(depth<lb[1]))
            self.data.append(tavg_data[iz].mean())
            print "Averaged layer %d-%d" % (lb[0],lb[1])
        self.depth = np.hstack((level_bounds[:,0],4000))
        print "depth:", self.depth

if __name__ == "__main__":
    lon, lat = 98., 83.
    hrs = Hiroshis()
    hrs.getPoint(lon,lat,'T')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5,10))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.ma.hstack((hrs.data,hrs.data[-1])),\
            -hrs.depth,\
            lw=2,drawstyle='steps-pre')
    ax.set_ylabel('depth [m]')
    plt.show()
    print "Finnished!"
