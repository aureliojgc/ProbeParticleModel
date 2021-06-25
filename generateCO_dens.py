#!/usr/bin/python
import sys
import numpy as np
import os
import __main__ as main


import pyProbeParticle                as PPU
#from   pyProbeParticle            import elements   
import pyProbeParticle.GridUtils      as GU
import pyProbeParticle.fieldFFT       as fFFT
import pyProbeParticle.HighLevel      as PPH
import pyProbeParticle.cpp_utils      as cpp_utils
import pyProbeParticle.basUtils       as BU
import pyProbeParticle.kpfmUtils      as kpfmU
from scipy.optimize import curve_fit

bohrRadius2angstroem = 0.5291772109217
rho=None 
rhox=None 
rhoy=None 
rhoz=None 

def generate_CO_density(lvec, nDim, head, borh=False, modulo=8.0, sigma=0.4):
    sampleSize = fFFT.getSampleDimensions( lvec )
    dims = (nDim[2], nDim[1], nDim[0])
    xsize, dx = fFFT.getSize('x', dims, sampleSize)
    ysize, dy = fFFT.getSize('y', dims, sampleSize)
    zsize, dz = fFFT.getSize('z', dims, sampleSize)
    dd = (dx, dy, dz)
    X, Y, Z = fFFT.getMGrid(dims, dd)
    fFFT.fieldInfo( Z, label="fieldInfo Z " )

    rho = fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=sigma, multipole_dict={'s':modulo})
    if (borh):
        rho = rho*(bohrRadius2angstroem**3)

    z = np.linspace(lvec[0,2],lvec[0,2]+lvec[3,2],nDim[0])
    rho_profile = rho[:,0,0]

    datos = np.transpose([z,rho_profile])

    np.savetxt('profile_sigma_'+str(sigma)+'_module_'+str(modulo), datos)


    return rho, lvec, nDim, head

