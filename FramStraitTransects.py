#!/usr/bin/env python
"""
Plot Polar ORA-IP annual mean transects (T and S)
across the Fram Strait:
    79N, -10E,20E
"""

import os, cPickle, gzip

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

    def findProduct(self,product):
        prodpoint = self.prset[self.lons[0]][self.lats[0]]
        prodidx = np.where(np.array([p.__class__==product \
                                     for p in prodpoint.products]))[0][0]
        return prodidx

    def plotOnePanel(self,product,ax):
        prodidx = self.findProduct(product)
        z = getattr(getattr(fstransect.prset[fstransect.lons[0]]\
                           [fstransect.lats[0]].products[prodidx],\
                            self.vname),'mz')
        if self.vname=='T':
            vmin, vmax = -2., 3.5
            mmin, mmax = -1.9, 10.
        else:
            vmin, vmax = 32., 35.2
            mmin, mmax = 32, 36.
        data = []
        for li, lon in enumerate(fstransect.lons):
            lat = self.lats[li]
            product = self.prset[lon][lat].products[prodidx]
            data.append(getattr(getattr(product,self.vname),'data'))
        dtrans = np.ma.masked_outside(np.ma.masked_equal(np.ma.array(data,mask=np.isnan(data)),0.),mmin,mmax).T
        cnt = ax.pcolormesh(self.lons,z,dtrans,cmap=self.ccmap,\
                            vmin=vmin,vmax=vmax)
        ax.invert_yaxis()
        ax.set_title(product.dset)
        return cnt

    def plotObjPanel(self,ax,objname):
        z = getattr(getattr(getattr(fstransect.prset[fstransect.lons[0]]\
                           [fstransect.lats[0]],objname),\
                            self.vname),'mz')
        if self.vname=='T':
            vmin, vmax = -2., 3.5
            mmin, mmax = -1.9, 10.
        else:
            vmin, vmax = 32., 35.2
            mmin, mmax = 32, 36.
        data = []
        for li, lon in enumerate(fstransect.lons):
            lat = self.lats[li]
            product = getattr(self.prset[lon][lat],objname)
            data.append(getattr(getattr(product,self.vname),'data'))
        dtrans = np.ma.masked_outside(np.ma.masked_equal(np.ma.array(data),0.),mmin,mmax).T
        cnt = ax.pcolormesh(self.lons,z,dtrans,cmap=self.ccmap,\
                            vmin=vmin,vmax=vmax)
        ax.invert_yaxis()
        ax.set_title(product.dset)
        return cnt

    def plotTransects(self,products):
        nx, ny = 4, 3
        fig = plt.figure(figsize=(15,20))
        for pidx, objname in enumerate(['sumata','woa13','mmm']):
            ax = fig.add_subplot(nx,ny,pidx+1)
            cnt = self.plotObjPanel(ax,objname)
        cax = fig.add_axes([0.92, 0.72, 0.02, 0.18])
        cb  = fig.colorbar(cnt,cax=cax)
        if self.vname=='T':
            cb.set_label("Temperature [$^\circ$C]")
        else:
            cb.set_label("Salinity [ppm]")
        for pidx, product in enumerate(products):
            ax = fig.add_subplot(nx,ny,pidx+4)
            cnt = self.plotOnePanel(product,ax)
            if pidx==6:
                ax.set_ylabel('depth [m]')
                ax.set_xlabel('longitude [$^\circ$E]')
        cax = fig.add_axes([0.16, 0.04, 0.7, 0.02])
        cb  = fig.colorbar(cnt,cax=cax,orientation='horizontal')
        if self.vname=='T':
            cb.set_label("Temperature [$^\circ$C]")
        else:
            cb.set_label("Salinity [ppm]")
        #plt.show()
        plt.savefig("fstransect_%s.png" % self.vname)

if __name__ == "__main__":
    prset = Products([UoR,GloSea5,MOVEG2i,GECCO2,EN4,\
                      ECDA,ORAP5,GLORYS2V4,CGLORS,SODA331,TOPAZ],'Fram Strait')
    #prset = Products([MOVEG2i,GECCO2,ORAP5,CGLORS,SODA331,TOPAZ],'Fram Strait')
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
    fstransect.plotTransects(prset)
    print "Finnished!"
