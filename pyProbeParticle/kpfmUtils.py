#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit
import GridUtils      as GU


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

def pos2index(r,nDim,lvec):

    delta = np.array([
        [lvec[1,0]/nDim[2], lvec[1,1]/nDim[1], lvec[1,2]/nDim[0]],
        [lvec[2,0]/nDim[2], lvec[2,1]/nDim[1], lvec[2,2]/nDim[0]],
        [lvec[3,0]/nDim[2], lvec[3,1]/nDim[1], lvec[3,2]/nDim[0]]
    ])

    #ix = int(round(((r[0]-lvec[0,0])/(lvec[1,0]))*nDim[2]))
    #iy = int(round(((r[1]-lvec[0,1])/(lvec[2,1]))*nDim[1]))
    denominador = delta[0,0]*delta[1,1] - delta[0,1]*delta[1,0]

    ix = int(round( ( (r[0]-lvec[0,0])*delta[1,1] - (r[1]-lvec[0,1])*delta[1,0] )/( denominador ) ))
    iy = int(round( ( (r[1]-lvec[0,1])*delta[0,0] - (r[0]-lvec[0,1])*delta[1,0] )/( denominador ) ))

    iz = int(round(((r[2]-lvec[0,2])/(lvec[3,2]))*nDim[0]))
    ir = np.array([ix,iy,iz])

    return ir

def index2pos(ir,nDim,lvec):

    delta = np.array([
        [lvec[1,0]/nDim[2], lvec[1,1]/nDim[1], lvec[1,2]/nDim[0]],
        [lvec[2,0]/nDim[2], lvec[2,1]/nDim[1], lvec[2,2]/nDim[0]],
        [lvec[3,0]/nDim[2], lvec[3,1]/nDim[1], lvec[3,2]/nDim[0]]
    ])

    x = lvec[0,0] + delta[0,0]*ir[0] + delta[1,0]*ir[1]
    y = lvec[0,1] + delta[0,1]*ir[0] + delta[1,1]*ir[1]
    z = lvec[0,2] + delta[2,2]*ir[2]

    return np.array([x,y,z])

def core_correction(V,lvec,nDim,Rs_tip,elems_tip,tol):
    
    V_corrected = V
    V_corrected_zero = V
    atoms_tip = np.array(Rs_tip)

    delta = np.array([lvec[1,0]/nDim[2], lvec[2,1]/nDim[1], lvec[3,2]/nDim[0]])

    pbc = 0 # int(2*tol/delta)

    for iatm in range(atoms_tip.shape[0]):
        print iatm

        iatm_pos = pos2index(np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]]), nDim, lvec)
        atom_pos = np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]])
        list_points = None
        list_dist = None

        print iatm_pos, atom_pos

        for i in range(0, nDim[2], 1):
            for j in range(0, nDim[1], 1):
                for k in range(0, nDim[0], 1):

                    r = index2pos([i,j,k],nDim,lvec)
                    vector = r - atom_pos
                    dist = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
                    #print i,j,k,dist
                    if (dist < tol):
                        if list_points is None:
                            print "initiate" , dist, i, j ,k
                            list_points = np.zeros([1,3], dtype=int)
                            list_dist = np.zeros([1], dtype=float)
                            list_points[0,:] = [i,j,k]
                            list_dist[0] = dist
                        else:
                            print "append" , dist, i, j ,k
                            list_points_aux = np.zeros([1,3], dtype=int)
                            list_dist_aux = np.zeros([1], dtype=float)
                            list_points_aux[0,:] = [i,j,k]
                            list_dist_aux[0] = dist
                            list_points=np.append(list_points, list_points_aux, axis=0)
                            list_dist=np.append(list_dist, list_dist_aux, axis=0)
        print "list_points"
        print list_points
        print "list_dist"
        print list_dist

        media = 0.0
        Norm = 0.0

        #i_min = iatm_pos[0]
        #i_max = iatm_pos[0]
        #j_min = iatm_pos[1]
        #j_max = iatm_pos[1]     
        #k_min = iatm_pos[2]
        #k_max = iatm_pos[2]

        for i in range(list_dist.shape[0]):

            media = media + list_dist[i]*V[ list_points[i,2], list_points[i,1], list_points[i,0] ]
            Norm = Norm + list_dist[i]

        #    if (list_points[0]>i_max): i_max = list_points[0]
        #    if (list_points[0]<i_min): i_min = list_points[0]
        #    if (list_points[1]>j_max): j_max = list_points[1]
        #    if (list_points[1]<j_min): j_min = list_points[1]
        #    if (list_points[2]>k_max): k_max = list_points[2]
        #    if (list_points[2]<k_min): k_min = list_points[2]           

        media = media/Norm
        #media_ring = ( V[i_min,iatm_pos[1],iatm_pos[2]] + V[i_max,iatm_pos[1],iatm_pos[2]] + V[iatm_pos[0],j_min,iatm_pos[2]] + V[iatm_pos[0],j_max,iatm_pos[2]] + V[iatm_pos[0],iatm_pos[1],k_min] + V[iatm_pos[0],iatm_pos[2],k_max] )/6.0

        for i in range(list_dist.shape[0]):
            print i
            r = index2pos([list_points[i,0], list_points[i,1], list_points[i,2]],nDim,lvec)
        #    V_corrected = media + (media_ring - media)*((r/tol)**2)
            #V_corrected[list_points[i,2], list_points[i,1], list_points[i,0]] = media + media*((tol - dist)/tol)*0.1
            V_corrected[list_points[i,2], list_points[i,1], list_points[i,0]] = media
            print "making point" , list_points[i,2], list_points[i,1], list_points[i,0] , "zero"
            V_corrected_zero[list_points[i,2], list_points[i,1], list_points[i,0]] = 0.0

        print iatm_pos, atom_pos

        print media

        print "list point and distances"
        print list_points
        print list_dist
    print "saving Corrected tip"
    np.savetxt('profile', np.array(V[:,0,0]) )
    GU.saveXSF( "Corrected_tip.xsf", V_corrected, lvec )
    np.savetxt('profile_corrected', np.array(V_corrected[:,0,0]) )
    GU.saveXSF( "Corrected_tip_zero.xsf", V_corrected_zero, lvec )
    np.savetxt('profile_corrected_zero', np.array(V_corrected_zero[:,0,0]) )

    return V_corrected_zero


