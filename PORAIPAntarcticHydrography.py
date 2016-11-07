#!/usr/bin/env python
"""
Plot Polar ORA-IP annual mean profiles (T and S)
and WOA13 1995-2004 profiles for the Southern Ocean.
Basically just excludes Sumata's Arctic climatology.
"""

from PORAIPHydrography import Products, UoR
from PORAIPHydrography import GloSea5, MOVEG2, GECCO2, EN4
from PORAIPHydrography import ECDA, ORAP5, TOPAZ, GLORYS2V4, CGLORS

if __name__ == "__main__":
    # Weddell Sea (deep)
    lon, lat = 330., -65.
    # Ross Sea (deep)
    lon, lat = 200., -70.
    # Amundsen Sea (deep)
    lon, lat = 240., -67.
    #prset = Products([UoR],lon,lat)
    prset = Products([UoR,GloSea5,MOVEG2,GECCO2,EN4,\
                      ECDA,ORAP5,GLORYS2V4,CGLORS],lon,lat)
    for vname in ['T','S']:
        prset.readProfiles(vname)
        prset.getMultiModelMean(vname)
        prset.plotDepthProfile(vname)
    prset.plotTSProfile()
    print "Finnished!"
