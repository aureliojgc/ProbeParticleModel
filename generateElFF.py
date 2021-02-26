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
from scipy.optimize import curve_fit


## def functions

def correct_background_old(V,Vbias,lvec,Vac_lvl):
    z = np.linspace(0,lvec[3,2],V.shape[0])
    hev = 0.0*z
    for i in range(z.shape[0]):

        if z[i] > Vac_lvl:
            hev[i] = 1.0

    y = -Vbias*z+Vbias*z[-1]*hev
    for i in range(z.shape[0]):
        V[i,:,:] = V[i,:,:] - y[i]
        if ((z[i] > Vac_lvl-0.3) and (z[i] < Vac_lvl+0.3)):
            V[i,:,:] = 0.0

    return V

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

if __name__=="__main__":
    HELP_MSG="""Use this program in the following way:
    """+os.path.basename(main.__file__) +""" -i <filename> [ --sigma <value> ]
    Supported file fromats are:
       * cube
       * xsf """
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "-i", "--input", action="store", type="string", help="format of input file")
    parser.add_option( "--tip_dens", action="store", type="string", default=None, help="tip enisty file (.xsf)" )
    #parser.add_option( "--sub_core",  action="store_true",  help="subtract core density", default=False )
    parser.add_option( "--Rcore",   default=PPU.params["Rcore"],    action="store", type="float", help="Width of nuclear charge density blob to achieve charge neutrality [Angstroem]" )
    parser.add_option( "-t", "--tip", action="store", type="string", help="tip model (multipole) {s,pz,dz2,..}", default=None)
    parser.add_option( "--tilt", action="store", type="float", help="tilt of tip electrostatic field (radians)", default=0 )
    parser.add_option( "-E", "--energy", action="store_true",  help="pbc False", default=False)
    parser.add_option("--noPBC", action="store_false",  help="pbc False",dest="PBC", default=None)
    parser.add_option( "-w", "--sigma", action="store", type="float",help="gaussian width for convolution in Electrostatics [Angstroem]", default=None)
    parser.add_option("-f","--data_format" , action="store" , type="string", help="Specify the output format of the vector and scalar field. Supported formats are: xsf,npy", default="xsf")
    parser.add_option("--KPFM_tip", action="store",type="string", help="read tip density under bias")
    parser.add_option("--KPFM_sample", action="store",type="string", help="read sample hartree under bias")
    parser.add_option("--Vref", action="store",type="float", help="Field under the KPFM dens. and Vh was calculated in V/Ang")
    parser.add_option("--z0", action="store",type="float", default=0.0 ,help="heigth of the topmost layer of metallic substrate for E to V conversion (Ang)")
    parser.add_option("--Vac", action="store", type="float", help="include a ramp function and a shift in the sample hartree under bias to compensate the introduced by the dft code", default=None)
    (options, args) = parser.parse_args()

    #print "options.tip_dens ", options.tip_dens;  exit() 

    if options.input is None:
        sys.exit("ERROR!!! Please, specify the input file with the '-i' option \n\n"+HELP_MSG)
    opt_dict = vars(options)

    if os.path.isfile( 'params.ini' ):
        FFparams=PPU.loadParams( 'params.ini' ) 
    else:
        print ">> LOADING default params.ini >> 's' ="  
        FFparams = PPU.loadParams( cpp_utils.PACKAGE_PATH+'/defaults/params.ini' )
    #PPU.loadParams( 'params.ini' )
    PPU.apply_options(opt_dict)    

    if os.path.isfile( 'atomtypes.ini' ):
        print ">> LOADING LOCAL atomtypes.ini"  
        FFparams=PPU.loadSpecies( 'atomtypes.ini' ) 
    else:
        FFparams = PPU.loadSpecies( cpp_utils.PACKAGE_PATH+'/defaults/atomtypes.ini' )

    if options.Rcore > 0.0:  # We do it here, in case it crash we don't want to wait for all the huge density files to load
        if options.tip_dens is None: raise Exception( " Rcore>0 but no tip density provided ! " )
        valElDict        = PPH.loadValenceElectronDict()
        Rs_tip,elems_tip = PPH.getAtomsWhichTouchPBCcell( options.tip_dens, Rcut=options.Rcore )

    atoms_samp,nDim_samp,lvec_samp = BU.loadGeometry( options.input, params=PPU.params )
    head_samp                      = BU.primcoords2Xsf( atoms_samp[0], [atoms_samp[1],atoms_samp[2],atoms_samp[3]], lvec_samp )

    V=None
    if(options.input.lower().endswith(".xsf") ):
        print ">>> loading Hartree potential from  ",options.input,"..."
        print "Use loadXSF"
        V, lvec, nDim, head = GU.loadXSF(options.input)
    elif(options.input.lower().endswith(".cube") ):
        print " loading Hartree potential from ",options.input,"..."
        print "Use loadCUBE"
        V, lvec, nDim, head = GU.loadCUBE(options.input)
    
    if PPU.params['tip']==".py":
        #import tip
        execfile("tip.py")
        print tipMultipole
        PPU.params['tip'] = tipMultipole
        print " PPU.params['tip'] ", PPU.params['tip']

    if options.tip_dens is not None:
        '''
        ###  NO NEED TO RENORMALIZE : fieldFFT already works with density
        rho_tip, lvec_tip, nDim_tip, head_tip = GU.loadXSF( options.tip_dens )
        rho_tip *= GU.dens2Q_CHGCARxsf(rho_tip, lvec_tip)
        PPU.params['tip'] = rho_tip
        print " dens_tip check_sum Q =  ", np.sum( rho_tip )
        '''
        print ">>> loading tip density from ",options.tip_dens,"..."

        if (options.tip_dens.lower().endswith(".xsf")):
            rho_tip, lvec_tip, nDim_tip, head_tip = GU.loadXSF( options.tip_dens )
        elif(options.tip_dens.lower().endswith(".cube")):
            rho_tip, lvec_tip, nDim_tip, head_tip = GU.loadCUBE( options.tip_dens, hartree=False )

        if options.Rcore > 0.0:
            print ">>> subtracting core densities from rho_tip ... "
            #subtractCoreDensities( rho_tip, lvec_tip, fname=options.tip_dens, valElDict=valElDict, Rcore=options.Rcore )
            PPH.subtractCoreDensities( rho_tip, lvec_tip, elems=elems_tip, Rs=Rs_tip, valElDict=valElDict, Rcore=options.Rcore, head=head_tip )

        PPU.params['tip'] = rho_tip

    if (options.KPFM_sample is not None):

        V_v0_aux = V.copy()
        V_v0_aux2 = V.copy()

        V_kpfm=None
        if(options.KPFM_sample.lower().endswith(".xsf") ):
            print ">>> loading Hartree potential  under bias from  ",options.KPFM_sample,"..."
            print "Use loadXSF"
            V_kpfm, lvec, nDim, head = GU.loadXSF(options.KPFM_sample)

            if options.Vac is not None:
                print "correnting bias on workfunction with a ramp function"
                V_kpfm = correct_background(V_kpfm,options.Vref,lvec,options.Vac)
                print "printing the hartree for test porpouses"
                GU.saveXSF('corrected_hartree.xsf',V_kpfm,lvec , head=head)
        elif(options.KPFM_sample.lower().endswith(".cube") ):
            print " loading Hartree potential under bias from ",options.KPFM_sample,"..."
            print "Use loadCUBE"
            V_kpfm, lvec, nDim, head = GU.loadCUBE(options.KPFM_sample)

            
        dV_kpfm = (V_kpfm - V_v0_aux)#/(options.Vref)

    if (options.KPFM_tip is not None):

        print ">>> loading tip density under bias from ",options.KPFM_tip,"..."

        if (options.KPFM_tip.lower().endswith(".xsf")):
            rho_tip_v0_aux = rho_tip.copy()
            rho_tip_kpfm, lvec_tip, nDim_tip, head_tip = GU.loadXSF( options.KPFM_tip )
            drho_kpfm = (rho_tip_kpfm - rho_tip_v0_aux)#/(options.Vref)
        elif(options.KPFM_tip.lower().endswith(".cube")):
            rho_tip_v0_aux = rho_tip.copy()
            rho_tip_kpfm, lvec_tip, nDim_tip, head_tip = GU.loadCUBE( options.KPFM_tip, hartree=False )
            drho_kpfm = (rho_tip_kpfm - rho_tip_v0_aux)#/(options.Vref)
        elif options.KPFM_tip in {'fit', 'dipole', 'pz'}:
            if ( PPU.params['probeType'] == 8 ):
                drho_kpfm={'pz':0.025} # compared with DFT VASP 0.015. As VASP goes with q=-1.0 and dz2 goes with -0.2 -> 5*0.015. COAg tip is over 0.025, so 0.125
                sigma = 0.7
            if ( PPU.params['probeType'] == 24 ):
                drho_kpfm={'pz':0.016}
                sigma = 0.8
                        #Calculate ~V terms
        FFkpfm_t0sV,Eel_t0sV=PPH.computeElFF(dV_kpfm,lvec,nDim,PPU.params['tip'],computeVpot=options.energy , tilt=opt_dict['tilt'] ,)
        FFkpfm_tVs0,Eel_tVs0=PPH.computeElFF(V_v0_aux2,lvec,nDim,drho_kpfm,computeVpot=options.energy , tilt=opt_dict['tilt'] )

        #debug save tippol
        if options.KPFM_tip in {'fit', 'dipole', 'pz'}:
            rho, lvec_tip, nDim_tip, tiphead = GU.loadXSF("rhoTip.xsf")
            GU.saveXSF( "Tip_bias_pol.xsf", rho, lvec )
        #debug save tippol

        zpos = np.linspace(lvec[0,2]-options.z0,lvec[3,2]-options.z0,nDim[0])
        for i in range(nDim[0]):
                FFkpfm_t0sV[i,:,:]=FFkpfm_t0sV[i,:,:]/((options.Vref)*(zpos[i]+0.1))
                FFkpfm_tVs0[i,:,:]=FFkpfm_tVs0[i,:,:]/((options.Vref)*(zpos[i]+0.1))

        print ">>> saving electrostatic forcefiled ... "
        GU.save_vec_field('FFkpfm_t0sV',FFkpfm_t0sV,lvec_samp ,data_format=options.data_format, head=head_samp)
        GU.save_vec_field('FFkpfm_tVs0',FFkpfm_tVs0,lvec_samp ,data_format=options.data_format, head=head_samp)


    print ">>> calculating electrostatic forcefiled with FFT convolution as Eel(R) = Integral( rho_tip(r-R) V_sample(r) ) ... "
    #FFel,Eel=PPH.computeElFF(V,lvec,nDim,PPU.params['tip'],Fmax=10.0,computeVpot=options.energy,Vmax=10, tilt=opt_dict['tilt'] )
    FFel,Eel=PPH.computeElFF(V,lvec,nDim,PPU.params['tip'],computeVpot=options.energy , tilt=opt_dict['tilt'] )

    print ">>> saving electrostatic forcefiled ... "
    
    GU.save_vec_field('FFel',FFel,lvec_samp ,data_format=options.data_format, head=head_samp)
    if options.energy:
        GU.save_scal_field( 'Eel', Eel, lvec_samp, data_format=options.data_format)
    del FFel,V;
    
    
    
    
