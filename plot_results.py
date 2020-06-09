#!/usr/bin/python -u

import os
import numpy as np
import matplotlib as mpl;  mpl.use('Agg'); print "plot WITHOUT Xserver"; # this makes it run without Xserver (e.g. on supercomputer) # see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt
import sys


'''
import basUtils
import elements
import GridUtils as GU
import ProbeParticleUtils as PPU
import PPPlot
'''

import pyProbeParticle                as PPU     
import pyProbeParticle.GridUtils      as GU
import pyProbeParticle.PPPlot         as PPPlot
from   pyProbeParticle            import basUtils
from   pyProbeParticle            import elements 
#import pyProbeParticle.core           as PPC
import pyProbeParticle.HighLevel      as PPH
import pyProbeParticle.cpp_utils      as cpp_utils


#import scipy                        # don't introduce unnecessary dependence on scipy
#from scipy.ndimage import laplace   # don't introduce unnecessary dependence on scipy

# =============== arguments definition

from optparse import OptionParser
parser = OptionParser()
parser.add_option( "-k",       action="store", type="float", help="tip stiffenss [N/m]" )
parser.add_option( "--krange", action="store", type="float", help="tip stiffenss range (min,max,n) [N/m]", nargs=3)
parser.add_option( "-q",       action="store", type="float", help="tip charge [e]" )
parser.add_option( "--qrange", action="store", type="float", help="tip charge range (min,max,n) [e]", nargs=3)
parser.add_option( "-a",       action="store", type="float", help="oscilation amplitude [A]" )
parser.add_option( "--arange", action="store", type="float", help="oscilation amplitude range (min,max,n) [A]", nargs=3)
parser.add_option( "--iets",   action="store", type="float", help="mass [a.u.]; bias offset [eV]; peak width [eV] ", nargs=3 )
parser.add_option( "-V","--Vbias",       action="store", type="float", help="Aplied field [eV/Ang]" )
parser.add_option( "--Vrange",  action="store", type="float", help="set of bias to perform the scan under", nargs=3)
parser.add_option( "--easy_KPFM_b",  action="store_true", default=False, help="calculate the b map in a fast way with the polariz. part of the El force" )
parser.add_option( "--LCPD_maps", action="store_true", default=False, help="print LCPD maps")
parser.add_option("--z0", action="store",type="float", default=0.0 ,help="heigth of the topmost layer of metallic substrate for E to V conversion (Ang)")
parser.add_option("--V2", action="store_false", default=True)

parser.add_option( "--df",       action="store_true", default=False,  help="plot images for dfz " )
parser.add_option( "--save_df" , action="store_true", default=False, help="save frequency shift as df.xsf " )
parser.add_option( "--Laplace",  action="store_true", default=False,  help="plot Laplace-filtered images and save them " )
parser.add_option( "--pos",      action="store_true", default=False, help="save probe particle positions" )
parser.add_option( "--atoms",    action="store_true", default=False, help="plot atoms to images" )
parser.add_option( "--bonds",    action="store_true", default=False, help="plot bonds to images" )
parser.add_option( "--cbar",     action="store_true", default=False, help="plot bonds to images" )
parser.add_option( "--WSxM",     action="store_true", default=False, help="save frequency shift into WsXM *.dat files" )
parser.add_option( "--bI",       action="store_true", default=False, help="plot images for Boltzmann current" )
parser.add_option("-f","--data_format" , action="store" , type="string",help="Specify the input/output format of the vector and scalar field. Supported formats are: xsf,npy", default="xsf")

parser.add_option( "--noPBC", action="store_false",  help="pbc False", dest="PBC",default=None)
(options, args) = parser.parse_args()
opt_dict = vars(options)

PPU.loadParams( 'params.ini' )
PPU.apply_options(opt_dict)

# =============== Setup

if opt_dict['Laplace']:
    from scipy.ndimage import laplace

#PPPlot.params = PPU.params

print " >> OVEWRITING SETTINGS by command line arguments  "
# Ks
if opt_dict['krange'] is not None:
    Ks = np.linspace( opt_dict['krange'][0], opt_dict['krange'][1], opt_dict['krange'][2] )
elif opt_dict['k'] is not None:
    Ks = [ opt_dict['k'] ]
else:
    Ks = [ PPU.params['klat'] ]
# Qs
if opt_dict['qrange'] is not None:
    Qs = np.linspace( opt_dict['qrange'][0], opt_dict['qrange'][1], opt_dict['qrange'][2] )
elif opt_dict['q'] is not None:
    Qs = [ opt_dict['q'] ]
else:
    Qs = [ PPU.params['charge'] ]
# Amps
if opt_dict['arange'] is not None:
    Amps = np.linspace( opt_dict['arange'][0], opt_dict['arange'][1], opt_dict['arange'][2] )
elif opt_dict['a'] is not None:
    Amps = [ opt_dict['a'] ]
else:
    Amps = [ PPU.params['Amplitude'] ]

    #activate the aplied bias
aplied_bias=False
if opt_dict['Vrange'] is not None:
    Vs = np.linspace( opt_dict['Vrange'][0], opt_dict['Vrange'][1], opt_dict['Vrange'][2] )
elif opt_dict['Vbias'] is not None:
    Vs = [ opt_dict['Vbias'] ]
else:
    Vs = 0.0
for iV,Vx in enumerate(Vs):
    if ( abs(Vx) > 1e-7):
        aplied_bias=True     

if (aplied_bias == True):
    print "Vs   =", Vs
print "Ks   =", Ks 
print "Qs   =", Qs 
print "Amps =", Amps 

#sys.exit("  STOPPED ")

print " ============= RUN  "

dz  = PPU.params['scanStep'][2]
xTips,yTips,zTips,lvecScan = PPU.prepareScanGrids( )
extent = ( xTips[0], xTips[-1], yTips[0], yTips[-1] )

atoms_str=""
atoms = None
bonds = None
FFparams = None
if opt_dict['atoms'] or opt_dict['bonds']:
    speciesFile=None
    if os.path.isfile( 'atomtypes.ini' ):
        speciesFile='atomtypes.ini'
    FFparams=PPU.loadSpecies( speciesFile ) 
    atoms_str="_atoms"
    atoms = basUtils.loadAtoms( 'input_plot.xyz' )
    #print "atoms ", atoms
    FFparams            = PPU.loadSpecies( ) 
    elem_dict           = PPU.getFFdict(FFparams);  #print "elem_dict ", elem_dict
    iZs,Rs,Qs_tmp=PPU.parseAtoms(atoms, elem_dict, autogeom = False, PBC = PPU.params['PBC'] )
    atom_colors = basUtils.getAtomColors(iZs,FFparams=FFparams)
    Rs=Rs.transpose().copy()
    atoms= [iZs,Rs[0],Rs[1],Rs[2],atom_colors]
    #print "atom_colors: ", atom_colors
if opt_dict['bonds']:
    bonds = basUtils.findBonds(atoms,iZs,1.0,FFparams=FFparams)
    #print "bonds ", bonds
atomSize = 0.15

cbar_str =""
if opt_dict['cbar']:
	cbar_str="_cbar"

if opt_dict['easy_KPFM_b']:
    try :
        f_kpfm, lvec, nDim = GU.load_scal_field('./OutKPFM_b' , data_format=options.data_format)
    except:
        print "error, canot load kpfm force map"
    df_kpfm = PPU.Fz2df( f_kpfm, dz = dz, k0 = PPU.params['kCantilever'], f0=PPU.params['f0Cantilever'], n=1.0/dz )
    PPPlot.plotImages(
        "./kpfm"+atoms_str+cbar_str, df_kpfm,  slices = range( 0, len(df_kpfm) ), zs=zTips+PPU.params['Amplitude']/2.0,
        extent=extent,cmap=PPU.params['colorscale_kpfm'], atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] 
    )

for iq,Q in enumerate( Qs ):
    for ik,K in enumerate( Ks ):
        dirname = "Q%1.2fK%1.2f" %(Q,K)
        for iv,Vx in enumerate( Vs ):
            if aplied_bias:
                dirname = "Q%1.2fK%1.2fV%1.2f" %(Q,K,Vx)
            if opt_dict['pos']:
                try:
                    PPpos, lvec, nDim = GU.load_vec_field(dirname+'/PPpos' ,data_format=options.data_format)
                    print " plotting PPpos : "
                    PPPlot.plotDistortions( 
                        dirname+"/xy"+atoms_str+cbar_str, PPpos[:,:,:,0], PPpos[:,:,:,1], slices = range( 0, len(PPpos) ), BG=PPpos[:,:,:,2], 
                        extent=extent, atoms=atoms, bonds=bonds, atomSize=atomSize, markersize=2.0, cbar=opt_dict['cbar'] 
                    )
                    del PPpos
                except:
                    print "error: ", sys.exc_info()
                    print "cannot load : " + ( dirname+'/PPpos_?.' + options.data_format ) 
            if opt_dict['iets'] is not None:
                try :
                    eigvalK, lvec, nDim = GU.load_vec_field( dirname+'/eigvalKs' ,data_format=options.data_format)
                    M  = opt_dict['iets'][0]
                    E0 = opt_dict['iets'][1]
                    w  = opt_dict['iets'][2]
                    print " plotting IETS M=%f V=%f w=%f " %(M,E0,w)	
                    hbar       = 6.58211951440e-16 # [eV.s]
                    aumass     = 1.66053904020e-27 # [kg] 
                    eVA2_to_Nm = 16.0217662        # [eV/A^2] / [N/m] 
                    Evib = hbar * np.sqrt( ( eVA2_to_Nm * eigvalK )/( M * aumass ) )
                    IETS = PPH.symGauss(Evib[:,:,:,0], E0, w) + PPH.symGauss(Evib[:,:,:,1], E0, w) + PPH.symGauss(Evib[:,:,:,2], E0, w)
                    PPPlot.plotImages( dirname+"/IETS"+atoms_str+cbar_str, IETS, slices = range(0,len(IETS)), zs=zTips, extent=extent, atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] )
                    PPPlot.plotImages( dirname+"/Evib"+atoms_str+cbar_str, Evib[:,:,:,0], slices = range(0,len(IETS)), zs=zTips, extent=extent, atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] )
                    PPPlot.plotImages( dirname+"/Kvib"+atoms_str+cbar_str, 16.0217662 * eigvalK[:,:,:,0], slices = range(0,len(IETS)), zs=zTips, extent=extent, atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] )
                    del eigvalK; del Evib; del IETS
                except:
                    print "error: ", sys.exc_info()
                    print "cannot load : ", dirname+'/PPpos_?.' + options.data_format
            if ( ( opt_dict['df'] or opt_dict['save_df'] or opt_dict['WSxM'] ) ):
                try :
                    fzs, lvec, nDim = GU.load_scal_field(dirname+'/OutFz' , data_format=options.data_format)

                    if (aplied_bias and opt_dict["V2"] ):
                        permit = PPU.params['permit']
                        #permit = 8.8541878176E-12 * (1.0/1.602176565E-19) * (1E-10)
                        #print "permit = ", permit
                        Rtip = PPU.params['Rtip']
                        for iz,z in enumerate( zTips ):
                            #print iz, z, Vx ,np.pi*permit*((Rtip*Rtip)/(z*(z+Rtip)))*Vx*Vx, fzs[iz,100,100], "bf"
                            fzs[iz,:,:] = fzs[iz,:,:] - np.pi*permit*((Rtip*Rtip)/((z-options.z0)*(z+Rtip)))*Vx*Vx
                            print np.pi*permit*((Rtip*Rtip)/((z-options.z0)*(z+Rtip)))*Vx*Vx, z
                            #print iz, z, Vx ,np.pi*permit*((Rtip*Rtip)/(z*(z+Rtip)))*Vx*Vx, fzs[iz,100,100], "af"

                    for iA,Amp in enumerate( Amps ):
                        AmpStr = "/Amp%2.2f" %Amp
                        print "Amp= ",AmpStr
                        dirNameAmp = dirname+AmpStr
                        if not os.path.exists( dirNameAmp ):
                            os.makedirs( dirNameAmp )
                        dfs = PPU.Fz2df( fzs, dz = dz, k0 = PPU.params['kCantilever'], f0=PPU.params['f0Cantilever'], n=Amp/dz )
                        if opt_dict['save_df']:
                            GU.save_scal_field(dirNameAmp+'/df', dfs, lvec,data_format=options.data_format )
                        if opt_dict['df']:
                            print " plotting df : "
                            PPPlot.plotImages(
                                dirNameAmp+"/df"+atoms_str+cbar_str, dfs,  slices = range( 0, len(dfs) ), zs=zTips+PPU.params['Amplitude']/2.0,
                                extent=extent,cmap=PPU.params['colorscale'], atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] 
                            )
                        if opt_dict['Laplace']:
                            print "plotting Laplace-filtered df : "
                            df_LaplaceFiltered = dfs.copy()
                            laplace( dfs, output = df_LaplaceFiltered )
                            GU.save_scal_field(dirNameAmp+'/df_laplace', df_LaplaceFiltered, lvec,data_format=options.data_format )
                            PPPlot.plotImages(
                                dirNameAmp+"/df_laplace"+atoms_str+cbar_str, df_LaplaceFiltered, slices = range( 0, len(dfs) ), zs=zTips+PPU.params['Amplitude']/2.0,
                                extent=extent,cmap=PPU.params['colorscale'], atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] 
                            )
                        if opt_dict['WSxM']:
                            print " printing df into WSxM files :"
                            GU.saveWSxM_3D( dirNameAmp+"/df" , dfs , extent , slices=None)

                        if opt_dict['LCPD_maps']:
                            if (iv == 0):
                                LCPD_b = - dfs
                            if (iv == (Vs.shape[0]-1)):
                                LCPD_b = (LCPD_b + dfs)/(2*Vx)
                            PPPlot.plotImages(
                                "./b_HzperV"+atoms_str+cbar_str, LCPD_b,  slices = range( 0, len(LCPD_b) ), zs=zTips+PPU.params['Amplitude']/2.0,
                                extent=extent,cmap=PPU.params['colorscale_kpfm'], atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] 
                            )

                            
                        del dfs
                    del fzs
                except:
                    print "error: ", sys.exc_info()
                    print "cannot load : ",dirname+'/OutFz.'+options.data_format
            if opt_dict['bI']:
                try:
                    I, lvec, nDim = GU.load_scal_field(dirname+'/OutI_boltzmann', data_format=options.data_format )
                    print " plotting Boltzmann current: "
                    PPPlot.plotImages( dirname+"/OutI"+atoms_str+cbar_str, I,  slices = range( 0,len(I) ), zs=zTips, extent=extent, atoms=atoms, bonds=bonds, atomSize=atomSize, cbar=opt_dict['cbar'] )
                    del I
                except:
                    print "error: ", sys.exc_info()
                    print "cannot load : " + (dirname+'/OutI_boltzmann.'+options.data_format ) 

print " ***** ALL DONE ***** "

#plt.show()  # for interactive plotting you have to comment "import matplotlib as mpl;  mpl.use('Agg');" at the end
