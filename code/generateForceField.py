#!/usr/bin/python -u
import sys
import numpy as np
import basUtils
import elements
import GridUtils     as GU
import ProbeParticle as PP;    PPU = PP.PPU;
import libFFTfin     as LFF
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-i", "--input",       action="store",        type="string", help="format of input file", default='vasp.locpot.xsf')
parser.add_option(       "--lj",          action="store_true",                  help="calculate Lennard-Jones force-field ",  default=False)
parser.add_option(       "--el",          action="store_true",                  help="calculate electrostatic force-field",   default=False)
parser.add_option(       "--all",         action="store_true",                  help="calculate all force-field",             default=True)
(options, args) = parser.parse_args()

num = len(sys.argv)
if (num < 2):
    sys.exit("Number of arguments = "+str(num-1)+". This script should have at least one argument. I am terminating...")
finput = sys.argv[num-1]

if(options.all and not options.el and not options.lj):
    options.el = True
    options.lj = True

sigma  = 1.0 # [ Angstroem ] 

print " >> OVEWRITING SETTINGS by params.ini  "
PPU.loadParams( 'params.ini' )



print " ========= get electrostatic forcefiled from hartree "

# TODO with time implement reading a hartree potential generated by different software
print "   + loading Hartree potential from disk "
if(options.input == 'vasp.locpot.xsf'):
    V, lvec, nDim, head = GU.loadXSF(finput)
elif(options.input == 'aims.cube'):
    V, lvec, nDim, head = GU.loadCUBE(finput)

print "     + update super-cell "
PPU.params['gridA'] = lvec[ 1,:  ].copy()
PPU.params['gridB'] = lvec[ 2,:  ].copy()
PPU.params['gridC'] = lvec[ 3,:  ].copy()
PPU.params['gridN'] = nDim.copy()


if(options.el):
    print " --- computing electrostatic forcefiled from hartree ---"
    print "   + computing convolution with tip by FFT "
    Fel_x,Fel_y,Fel_z = LFF. potential2forces( V, lvec, nDim, sigma = 1.0 )
    
    print "   + saving electrostatic forcefiled into *.xsf files"
    GU.saveXSF('FFel_x.xsf', Fel_x, lvec, head)
    GU.saveXSF('FFel_y.xsf', Fel_y, lvec, head)
    GU.saveXSF('FFel_z.xsf', Fel_z, lvec, head)

    del Fel_x,Fel_y,Fel_z,V
    print ""


if(options.lj):
    print "--- computing Lennard-Jones force-filed ---"
    atoms     = basUtils.loadAtoms('input.xyz', elements.ELEMENT_DICT )
    iZs,Rs,Qs = PP.parseAtoms( atoms, autogeom = False, PBC = True )
    FFLJ      = PP.computeLJ( Rs, iZs, FFLJ=None, FFparams=None)

    GU.limit_vec_field( FFLJ, Fmax=100.0 ) # remove too large values; keeps the same direction; good for visualization 

    print "   + saving  LJ Force-filed into *.xsf files ---"
    GU.saveVecFieldXsf( 'FFLJ', FFLJ, lvec, head)

