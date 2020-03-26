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


from optparse import OptionParser

## === Functions

def ShiftZpot(V,dd,shift):
    ishift = int(round(shift/dd[2]))
    print ishift
    V = np.roll(V,ishift,axis=0)

    return V

    

## ==============


parser = OptionParser()
parser.add_option( "-s", "--sample_type", action="store", type="string", help="tip model (multipole) {s,pz,dz2,..}", default=None)
parser.add_option( "--charge",       action="store", type="float", help="tip charge [e]", default=1.0 )
parser.add_option(  "--sig", action="store", type="float", default=0.7)
parser.add_option( "--shift", action="store", type="float", help="Shift of the potential in the z direction", default=0.0)
parser.add_option( "--core", action="store", type="float", help="substracts a stiffer s orbital",  nargs=2, default=None)
parser.add_option( "--KPFM", action="store", type="float")
parser.add_option( "--Vbias", action="store", type="float")

(options, args) = parser.parse_args()

sample_type = options.sample_type

if sample_type in {'s','px','py','pz','dx2','dy2','dz2','dxy','dxz','dyz'}:
    rho = None
    multipole={sample_type:1.0}

if sample_type==".py":
    #import tip
    execfile("tip.py")
    print tipMultipole
    sample_type = tipMultipole
    print " PPU.params['tip'] ", sample_type
    rho = None
    multipole = sample_type


lvec = np.array(([0.0,0.0,0.0],[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,20.0]))

sigma = options.sig
tilt = 0.0

nDim = np.array((int(round(lvec[1,0]/0.1)),int(round(lvec[2,1]/0.1)),int(round(lvec[3,2]/0.1))))

sampleSize = fFFT.getSampleDimensions( lvec )
dims = (nDim[0], nDim[1], nDim[2])
    
xsize, dx = fFFT.getSize('x', dims, sampleSize)
ysize, dy = fFFT.getSize('y', dims, sampleSize)
zsize, dz = fFFT.getSize('z', dims, sampleSize)
dd = (dx, dy, dz)

X, Y, Z = fFFT.getMGrid(dims, dd)
if rho is None:
    rho = fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=sigma, multipole_dict=multipole, tilt=tilt )


    #rho[:,:,:] = rho[::-1,:,:] 
    rho =  rho*options.charge

    shift = options.shift
    if shift > 0.0:
        print shift
        #ishift = 10
        rho = ShiftZpot(rho,dd,shift)
        #np.roll(rho,ishift,axis=0)

    if options.core is not None:
        core_sigma = options.core[0]
        core_charge = options.core[1]

        core_rho = core_charge*fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=core_sigma, multipole_dict={"s":1.0}, tilt=tilt )
        rho = rho + core_rho    

    if options.KPFM is not None:

        Vbias = options.Vbias
        kpfm_charge = options.KPFM
        
        KPFM_tip = -1*Vbias*kpfm_charge*fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=0.7, multipole_dict={"pz":1.0}, tilt=tilt )
        KPFM_sample = -1*Vbias*kpfm_charge*fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=0.7, multipole_dict={"pz":1.0}, tilt=tilt )

        tip_check = KPFM_tip
        sample_check = KPFM_sample

        KPFM_sample = KPFM_sample*abs(options.charge)

        tip = fFFT.getProbeDensity(sampleSize, X, Y, Z, dd, sigma=0.7, multipole_dict={"dz2":1.0}, tilt=tilt )

        tip_check = KPFM_tip
        sample_check = KPFM_sample

        KPFM_sample = KPFM_sample + rho
        KPFM_tip = KPFM_tip + tip

    XSF_HEAD = '''CRYSTAL
PRIMVEC ''' + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[1,0],lvec[1,1],lvec[1,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[2,0],lvec[2,1],lvec[2,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[3,0],lvec[3,1],lvec[3,2]) + '''    
CONVVEC''' + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[1,0],lvec[1,1],lvec[1,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[2,0],lvec[2,1],lvec[2,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[3,0],lvec[3,1],lvec[3,2]) + '''    
PRIMCOORD
1  1
9   0.0   0.0   0.0

BEGIN_BLOCK_DATAGRID_3D                        
some_datagrid      
    BEGIN_DATAGRID_3D_whatever 
'''

    GU.saveXSF( "sample.xsf", rho, lvec, head=XSF_HEAD )

    if options.KPFM is not None:

        GU.saveXSF( "sample_kpfm.xsf", KPFM_sample, lvec, head=XSF_HEAD )
        GU.saveXSF( "sample_kpfm_check.xsf", sample_check, lvec, head=XSF_HEAD )

        XSF_HEAD = '''CRYSTAL
PRIMVEC ''' + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[1,0],lvec[1,1],lvec[1,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[2,0],lvec[2,1],lvec[2,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[3,0],lvec[3,1],lvec[3,2]) + '''    
CONVVEC''' + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[1,0],lvec[1,1],lvec[1,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[2,0],lvec[2,1],lvec[2,2]) + '''
    %2.5f  %2.5f  %2.5f  ''' %(lvec[3,0],lvec[3,1],lvec[3,2]) + '''    
PRIMCOORD
2  1
6   0.0   0.0   1.15
8   0.0   0.0   0.0

BEGIN_BLOCK_DATAGRID_3D                        
some_datagrid      
    BEGIN_DATAGRID_3D_whatever 
'''

        GU.saveXSF( "tip_kpfm.xsf", KPFM_tip, lvec, head=XSF_HEAD )
        GU.saveXSF( "tip.xsf", tip, lvec, head=XSF_HEAD )

        GU.saveXSF( "tip_kpfm_check.xsf", tip_check, lvec, head=XSF_HEAD )





