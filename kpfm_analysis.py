import os

import numpy as np
import scipy.integrate
import pyProbeParticle.GridUtils as GU
import matplotlib.pyplot as plt
import matplotlib.cm as cm 

import pyProbeParticle as PPU     


from scipy.optimize import curve_fit


def get_LCPD(image_point,V,Q,K,Amp):
    df_V = np.zeros([ image_point.shape[0], V.shape[0] ] )

    
    for iv,Vx in enumerate( V ):
        
        dirname = "Q%1.2fK%1.2fV%1.2f/Amp%1.2f/" %(Q,K,Vx,Amp)
        df, image_lvec, image_nDim, image_head = GU.loadXSF(dirname+'df.xsf')
        for i in range(image_point.shape[0]):
            image_relative_point = np.array([image_point[i,0]/image_lvec[1,0],image_point[i,1]/image_lvec[2,1],(image_point[i,2]-image_lvec[0,2])/image_lvec[3,2]])
        
            image_mtx_elementx = int(round(image_relative_point[0]*df.shape[2]))
            image_mtx_elementy = int(round(image_relative_point[1]*df.shape[1]))
            image_mtx_elementz = int(round(image_relative_point[2]*df.shape[0]))
        
            df_V[i,iv] =  df[image_mtx_elementz,image_mtx_elementy,image_mtx_elementx]

            df_fit=np.zeros([image_point.shape[0],4])
            x_axis = np.linspace(-1,1,200)

            acd = np.zeros((image_point.shape[0],3))
            fit_data = np.zeros((image_point.shape[0],x_axis.shape[0]))
            x = V
            y = df_V[i,:]

            def fit_func(x, a, d, c):
                return a*x**2 + d*x+c

            params = curve_fit(fit_func, x, y)

            [a, d, c] = params[0]

            acd[i,:] = params[0]

            fit_data[i,:] = a*(x_axis**2) + d*x_axis + c

            df_fit[i,0] = image_point[i,2]
            df_fit[i,1] = acd[i,0]
            df_fit[i,2] = acd[i,1]
            df_fit[i,3] = acd[i,2]


            
    return df_V, df_fit, fit_data, x_axis



from optparse import OptionParser
parser = OptionParser()

parser.add_option( "-k",       action="store", type="float", help="tip stiffenss [N/m]" )
parser.add_option( "--krange", action="store", type="float", help="tip stiffenss range (min,max,n) [N/m]", nargs=3)
parser.add_option( "-q",       action="store", type="float", help="tip charge [e]" )
parser.add_option( "--qrange", action="store", type="float", help="tip charge range (min,max,n) [e]", nargs=3)
parser.add_option( "-a",       action="store", type="float", help="oscilation amplitude [A]" )
parser.add_option( "--arange", action="store", type="float", help="oscilation amplitude range (min,max,n) [A]", nargs=3)

parser.add_option( "--Vrange",  action="store", type="float", help="set of bias to perform the scan under", nargs=3)
parser.add_option( "--image_point", action="store", type="float", help="coordinates on the image for LCPD calculation", nargs=3)
parser.add_option( "--save_LCPD" , action="store_true", default=False, help="save LCPD as data file " )

(options, args) = parser.parse_args()
opt_dict = vars(options)

PPU.loadParams( 'params.ini' )
PPU.apply_options(opt_dict)

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


if opt_dict['image_point']:
    image_point = np.zeros((1,3))
    image_point[0,0] = opt_dict['image_point'][0]
    image_point[0,1] = opt_dict['image_point'][1]
    image_point[0,2] = opt_dict['image_point'][2]
#else if multiple image points:
    #here deal with a list of points
else:
    print 'Scan point unspecified'

Vs = np.linspace( opt_dict['Vrange'][0], opt_dict['Vrange'][1], opt_dict['Vrange'][2] )


for iq,Q in enumerate( Qs ):
    for ik,K in enumerate( Ks ):
        for iA,Amp in enumerate( Amps ):
            df_V, df_fit, df_fitted, V_fitted = get_LCPD(image_point,Vs,Q,K,Amp)

            for i in range(image_point.shape[0]):

                x = opt_dict['image_point'][0]
                y = opt_dict['image_point'][1]
                z = opt_dict['image_point'][2]

                plt.plot(V_fitted[:], df_fitted[i,:], linewidth=3)
                plt.plot(Vs[:],df_V[i,:], label="x%1.2fy%1.2fz%1.2f" %(x,y,z), marker="o")
                plt.plot(Vs[:],df_fit[i,2]*Vs[:]+df_fit[i,3])
                plt.legend(["x%1.2fy%1.2fz%1.2f"  %(x,y,z)])


            plt.xlabel('bias (V)', fontsize=10)
            plt.ylabel('df (Hz)',fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)

            filename = "Q%1.2fK%1.2fAmp%1.2fx%1.2fy%1.2fz%1.2f" %(Q,K,Amp,x,y,z)
            plt.savefig('LCPD'+filename+'.png') 


            if opt_dict['save_LCPD']:
                outdata = np.append([Vs], df_V, axis=0)
                outdata = np.transpose(outdata)
                np.savetxt('LCPD'+filename+'.dat', outdata)
                np.savetxt('LCPD_fit'+filename+'.dat', df_fit)