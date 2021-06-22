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

from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-i", "--input", action="store", type="string", default="CHGCAR.xsf", help="sample 3D data-file (.xsf)")
parser.add_option( "--borh", action="store_true", help="the input is in a.u.", default=False )
parser.add_option( "--module", action="store", type="float", default="1.0" )
parser.add_option( "--sigma", action="store", type="float", default="0.7" )

(options, args) = parser.parse_args()

bohrRadius2angstroem = 0.5291772109217
sigma = options.sigma
rho=None 
rhox=None 
rhoy=None 
rhoz=None 
tilt=0.0

V=None
if(options.input.lower().endswith(".xsf") ):
    print ">>> loading Hartree potential from  ",options.input,"..."
    print "Use loadXSF"
    V, lvec, nDim, head = GU.loadXSF(options.input)
elif(options.input.lower().endswith(".cube") ):
    print " loading Hartree potential from ",options.input,"..."
    print "Use loadCUBE"
    V, lvec, nDim, head = GU.loadCUBE(options.input)

sampleSize = fFFT.getSampleDimensions( lvec )
dims = (nDim[2], nDim[1], nDim[0])
xsize, dx = fFFT.getSize('x', dims, sampleSize)
ysize, dy = fFFT.getSize('y', dims, sampleSize)
zsize, dz = fFFT.getSize('z', dims, sampleSize)
dd = (dx, dy, dz)
X, Y, Z = fFFT.getMGrid(dims, dd)
fFFT.fieldInfo( Z, label="fieldInfo Z " )

rho = fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=options.sigma, multipole_dict={'s':options.module}, tilt=tilt)
if (options.borh):
    rho = rho*(bohrRadius2angstroem**3)

if(options.input.lower().endswith(".xsf") ):
    GU.saveXSF("s_density.xsf", rho, lvec )
elif(options.input.lower().endswith(".cube") ):
    print "soubrutine to write cube files to be developed"

z = np.linspace(lvec[0,2]-options.xyz_shift[2],lvec[0,2]+lvec[3,2]-options.xyz_shift[2],nDim[0])
rho_profile = rho[:,0,0]

datos = [z,rho_profile]

np.savetxt('profile_sigma_'+str(options.sigma)+'_module_'+str(options.module))
