#!/usr/bin/env python
"""
Plot Polar ORA-IP annual mean transects (T and S)
across the Fram Strait:
    79N, -10E,20E
"""

import numpy as np

from PORAIPHydrography import Products, UoR
from PORAIPHydrography import GloSea5, MOVEG2, GECCO2, EN4
from PORAIPHydrography import ECDA, ORAP5, TOPAZ, GLORYS2V4, CGLORS

if __name__ == "__main__":
    # Fram Strait collection of points
    lons = np.arange(-10.5,21,1)
    lats = np.array([79]*len(lons))
    for li, lon in enumerate(lons):
        lat = lats[li]
        prset = Products([UoR],lon,lat)
        #prset = Products([UoR,GloSea5,MOVEG2,GECCO2,EN4,\
        #                  ECDA,ORAP5,GLORYS2V4,CGLORS],lon,lat)
        for vname in ['T','S']:
            prset.readProfiles(vname)
    print "Finnished!"
