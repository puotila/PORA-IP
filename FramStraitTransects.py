#!/usr/bin/env python
"""
Plot Polar ORA-IP annual mean transects (T and S)
across the Fram Strait:
    79N, -10E,20E
"""

import sys, os, cPickle, gzip

import numpy as np
import matplotlib.pyplot as plt
import cmocean

from PORAIPHydrography import Products, UoR, MultiModelMean, WOA13, Sumata
from PORAIPHydrography import GloSea5, MOVEG2i, GECCO2, EN4
from PORAIPHydrography import ECDA, ORAP5, SODA331, TOPAZ, GLORYS2V4, CGLORS

class Transect(object):
    def __init__(self,products):
        self.prset = prset
        self.prset.readTransects()
        for vname in ['T','S']:
            self.prset.getMultiModelMean(vname)

    def plotOnePanel(self,ax,product,vname,\
                     lons=np.arange(-19,22),\
                     xidx=np.hstack((np.arange(-19,0),np.arange(22)))):
        vari    = getattr(product,vname)
        z = getattr(vari,'mz')
        if vname=='T':
            ccmap= cmocean.cm.thermal
            vmin, vmax = -2., 3.5
        #    mmin, mmax = -1.9, 10.
        else:
            ccmap= cmocean.cm.haline
            vmin, vmax = 32., 35.
        #    mmin, mmax = 32, 36.
        dtrans = getattr(vari,'data')
        cnt = ax.pcolormesh(lons,z,dtrans[:,xidx],vmin=vmin,vmax=vmax,cmap=ccmap)
        ax.invert_yaxis()
        ax.set_title(product.dset)
        return cnt

    def plotTransects(self,vname):
        products = self.prset.products
        nx, ny = 4, 4
        fig = plt.figure(figsize=(15,20))
        for pidx, objname in enumerate(['sumata','woa13','mmm']):
            product = getattr(self.prset,objname)
            ax = fig.add_subplot(nx,ny,pidx+1)
            if objname in ['mmm']:
                cnt = self.plotOnePanel(ax,product,vname)
            else:
                cnt = self.plotOnePanel(ax,product,vname,xidx=np.arange(160,201))
        cax = fig.add_axes([0.92, 0.72, 0.02, 0.18])
        cb  = fig.colorbar(cnt,cax=cax)
        if vname=='T':
            cb.set_label("Temperature [$^\circ$C]")
        else:
            cb.set_label("Salinity [ppm]")
        for pidx, product in enumerate(products):
            ax = fig.add_subplot(nx,ny,pidx+4)
            cnt = self.plotOnePanel(ax,product,vname)
            if pidx==6:
                ax.set_ylabel('depth [m]')
                ax.set_xlabel('longitude [$^\circ$E]')
        cax = fig.add_axes([0.16, 0.04, 0.7, 0.02])
        cb  = fig.colorbar(cnt,cax=cax,orientation='horizontal')
        if vname=='T':
            cb.set_label("Temperature [$^\circ$C]")
        else:
            cb.set_label("Salinity [ppm]")
        #plt.show()
        plt.savefig("fstransect_%s.png" % vname)

if __name__ == "__main__":
    prset = Products([UoR,GloSea5,MOVEG2i,GECCO2,EN4,\
                      ECDA,ORAP5,GLORYS2V4,CGLORS,SODA331,TOPAZ],'Fram Strait')
    #prset = Products([UoR],'Fram Strait')
    fileout = "fstransect"
    cpzfile = fileout+'.cpickle.gz'
    if os.path.exists(cpzfile):
        fp = gzip.open(cpzfile)
        fstransect = cPickle.load(fp)
        fp.close()
    else:
        fstransect = Transect(prset)
        fp = gzip.open(cpzfile,'w')
        cPickle.dump(fstransect,fp)
        fp.close()
    # plotting
    for vname in ['T','S']:
        fstransect.plotTransects(vname)
    print "Finnished!"
