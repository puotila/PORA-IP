#!/usr/bin/env python
"""
Plot ORA-IP profile points on a map. Just for a location illustration.
"""
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from netCDF4 import Dataset

# read in etopo5 topography/bathymetry.
url = 'http://iridl.ldeo.columbia.edu/SOURCES/.WORLDBATH/data.nc'
etopodata = Dataset('etopo5.nc')

topoin = etopodata.variables['bath'][:]
lons = etopodata.variables['X'][:]
lats = etopodata.variables['Y'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

# plot topography/bathymetry as an image.

# create the figure and axes instances.
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(projection='npstere',boundinglat=70,lon_0=0,resolution='h',ax=ax)
# transform to nx x ny regularly spaced 5km native projection grid
nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map with imshow.
im = m.imshow(topodat,cm.Blues_r,vmax=0)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.fillcontinents(color='grey')
# draw parallels and meridians.
# label on left and bottom of map.
parallels = np.arange(50.,90,5.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = np.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[0,0,0,1])
#plot locations
xp,yp = m([10,100,220],[88,83,80])
m.scatter(xp,yp,color='red',s=60,edgecolors='k')
# add colorbar
cb = m.colorbar(im,"right", size="5%", pad='2%')
#ax.set_title('ETOPO5 Topography - Lambert Conformal Conic')
plt.savefig('PORA-IP-Arctic-points.png')
#plt.show()

