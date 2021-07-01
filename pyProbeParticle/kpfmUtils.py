#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit


def correct_background(V,Vbias,lvec,Vac_lvl):
    z = np.linspace(0,lvec[3,2],V.shape[0])
    wf = V.mean(-1).mean(-1)

    def correcttion_function(z, a, b, c, d, z0):
        return a*z + b/(np.exp((z-z0)/d)+1) + c

    #initial parameters
    a0 = Vbias
    b0 = (lvec[3,2]-lvec[0,2])*Vbias
    c0 = 0.0
    d0 = 0.1
    z00 = Vac_lvl
    param0 = [a0,b0,c0,d0,z00]

    #select fitting area
    zmin = Vac_lvl - 3.0
    zmax = Vac_lvl + 3.0

    iz_min = int(round(V.shape[0]*(zmin-lvec[0,2])/lvec[3,2]))
    iz_max = int(round(V.shape[0]*(zmax-lvec[0,2])/lvec[3,2]))   

    params = curve_fit(correcttion_function,z[iz_min:iz_max],wf[iz_min:iz_max],p0=param0)

    [a,b,c,d,z0] = params[0]
    print params[0]

    y = correcttion_function(z, a, b, c, d, z0)

    np.savetxt('correction_function', y)

    for i in range(z.shape[0]):
        V[i,:,:] = V[i,:,:] - y[i]
        if (i > iz_min) and (i < iz_max):
              V[i,:,:] = 0.0

    return V