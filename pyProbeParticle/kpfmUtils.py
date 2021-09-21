#!/usr/bin/python

import numpy as np
from numpy.core.numeric import Inf
from numpy.lib.function_base import piecewise
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

def pbc(i,imax):
    if (i > imax-1): i0 = i - imax
    elif (i < 0): i0 = i + imax
    else: i0 = i

    return i0

def TD_gaussian_fit(atom_pos, points_positions, points_value, A0, point2plot ):

    def gauss(R, x0, y0, z0, A, sx, sy, sz):
        x, y, z = R
        fx = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( (x - x0)**2 )/(2*sx**2) )
        fy = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( (y - y0)**2 )/(2*sy**2) )
        fz = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( (z - z0)**2 )/(2*sz**2) )
        G = A*fx*fy*fz
        return G
    def gauss0(R, A, sx, sy, sz):
        x, y, z = R
        fx = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( x**2 )/(2*sx**2) )
        fy = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( y**2 )/(2*sy**2) )
        fz = (1.0/(np.sqrt(2*np.pi)))*np.exp( -( z**2 )/(2*sz**2) )
        G = A*fx*fy*fz
        return G

    #initialize parameters 
    s0 = [A0, 0.5,0.5,0.5]
    smax = [Inf, 10.0, 10.0, 10.0]
    smin = [0.0, 0.2, 0.2, 0.2]

    #R = points_positions
    #R[:,0], R[:,1], R[:,2] = [points_positions[:,0] - atom_pos[0], points_positions[:,1] - atom_pos[1], points_positions[:,2] - atom_pos[2]]
    #R = points_positions
    x_norm = points_positions[:,0] - atom_pos[0]
    y_norm = points_positions[:,1] - atom_pos[1]
    z_norm = points_positions[:,2] - atom_pos[2]
    params, junk = curve_fit(gauss0, (x_norm, y_norm, z_norm), points_value, p0=s0, bounds= (smin, smax))

    #print "A, sx, sy, sz", params
    [A, sx, sy, sz] = params
    #print "center of the gaussian is", gauss0([0.0,0.0,0.0], A, sx, sy, sz )
    #print    point2plot - atom_pos
    V = gauss0(point2plot - atom_pos, A, sx, sy, sz )

    return V


def core_correction(V,lvec,nDim,Rs_tip,elems_tip, tol_in ,tol_out):
    print "V with shape ", V.shape, "and nDim", nDim 
    V_corrected = V.copy()
    atoms_tip = np.array(Rs_tip)

    delta = np.array([lvec[1,0]/nDim[2], lvec[2,1]/nDim[1], lvec[3,2]/nDim[0]])

    for iatm in range(atoms_tip.shape[0]):
        print iatm

        iatm_pos = pos2index(np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]]), nDim, lvec)
        atom_pos = np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]])
        list_points = None
        list_dist = None
        list_pos = None
        list_V = None
        list_points_inner = None
        list_dist_inner = None

        check_area = np.zeros([3], dtype=int)
        check_area[0] = int((3*tol_out)/delta[0])
        check_area[1] = int((3*tol_out)/delta[1])
        check_area[2] = int((3*tol_out)/delta[2])

        limits = np.array( [
            [iatm_pos[0]-check_area[0], iatm_pos[0]+check_area[0]],
            [iatm_pos[1]-check_area[1], iatm_pos[1]+check_area[1]],
            [iatm_pos[2]-check_area[2], iatm_pos[2]+check_area[2]]
            ] )

        for i in range(limits[0,0], limits[0,1], 1):
            for j in range(limits[1,0], limits[1,1], 1):
                for k in range(limits[2,0], limits[2,1], 1):   

                    r = index2pos([i,j,k],nDim,lvec)
                    vector = r - atom_pos
                    dist = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

                    if ((dist > tol_in) and (dist < tol_out)):
                        if list_points is None:
                            #print "initiate" , dist, i, j ,k
                            list_points = np.zeros([1,3], dtype=int)
                            list_dist = np.zeros([1], dtype=float)
                            list_pos = np.zeros([1,3], dtype=float)
                            list_V = np.zeros([1], dtype=float)

                            list_points[0,:] = [i,j,k]
                            list_dist[0] = dist
                            list_pos[0,:] = index2pos( [i,j,k] ,nDim,lvec)
                            list_V[0] = V[pbc(k,nDim[0]), pbc(j,nDim[1]), pbc(i,nDim[2])]

                        else:
                            #print "append" , dist, i, j ,k
                            list_points_aux = np.zeros([1,3], dtype=int)
                            list_dist_aux = np.zeros([1], dtype=float)
                            list_pos_aux = np.zeros([1,3], dtype=float)
                            list_V_aux = np.zeros([1], dtype=float)

                            list_points_aux[0,:] = [i,j,k]
                            list_dist_aux[0] = dist
                            list_pos_aux[0,:] = index2pos( [i,j,k] ,nDim,lvec)
                            list_V_aux[0] = V[pbc(k,nDim[0]), pbc(j,nDim[1]), pbc(i,nDim[2])]

                            list_points=np.append(list_points, list_points_aux, axis=0)
                            list_dist=np.append(list_dist, list_dist_aux, axis=0)    
                            list_pos=np.append(list_pos, list_pos_aux, axis=0)
                            list_V=np.append(list_V, list_V_aux, axis=0)
                    
                    if (dist < tol_in):
                        if list_points_inner is None:
                            #print "initiate" , dist, i, j ,k
                            list_points_inner = np.zeros([1,3], dtype=int)
                            list_dist_inner = np.zeros([1], dtype=float)
                            list_pos_inner = np.zeros([1,3], dtype=float)

                            list_points_inner[0,:] = [i,j,k]
                            list_dist_inner[0] = dist
                            list_pos_inner[0,:] = index2pos( [i,j,k] ,nDim,lvec)

                        else:
                            #print "append" , dist, i, j ,k
                            list_points_aux_inner = np.zeros([1,3], dtype=int)
                            list_dist_aux_inner = np.zeros([1], dtype=float)
                            list_pos_aux_inner = np.zeros([1,3], dtype=float)

                            list_points_aux_inner[0,:] = [i,j,k]
                            list_dist_aux_inner[0] = dist
                            list_pos_aux_inner[0,:] = index2pos( [i,j,k] ,nDim,lvec)

                            list_points_inner=np.append(list_points_inner, list_points_aux_inner, axis=0)
                            list_dist_inner=np.append(list_dist_inner, list_dist_aux_inner, axis=0)    
                            list_pos_inner=np.append(list_pos_inner, list_pos_aux_inner, axis=0)



        Centro = 50.0
        #print "data2fit_x, data2fit_y",   data2fit_x, data2fit_y
        for i in range(list_dist_inner.shape[0]):
            V_corrected[ 
                pbc(list_points_inner[i,2],nDim[0]), 
                pbc(list_points_inner[i,1],nDim[1]), 
                pbc(list_points_inner[i,0],nDim[2]) ] = TD_gaussian_fit(atom_pos, list_pos, list_V, Centro, [list_pos_inner[i,0], list_pos_inner[i,1], list_pos_inner[i,2]])
            #print "fitted value" , V_corrected[ pbc(list_points_inner[i,2],nDim[0]), pbc(list_points_inner[i,1],nDim[1]), pbc(list_points_inner[i,0],nDim[2]) ]

    GU.saveXSF( "Corrected_tip.xsf", V_corrected, lvec )
    return V_corrected



def core_correction_oldnew(V,lvec,nDim,Rs_tip,elems_tip,tol):
    print pbc(301,nDim[0])
    print "V with shape ", V.shape, "and nDim", nDim 
    V_corrected = V.copy()
    atoms_tip = np.array(Rs_tip)

    delta = np.array([lvec[1,0]/nDim[2], lvec[2,1]/nDim[1], lvec[3,2]/nDim[0]])

    for iatm in range(atoms_tip.shape[0]):
        print iatm

        iatm_pos = pos2index(np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]]), nDim, lvec)
        atom_pos = np.array([atoms_tip[iatm,0],atoms_tip[iatm,1],atoms_tip[iatm,2]])
        list_points = None
        list_dist = None

        check_area = np.zeros([3], dtype=int)
        check_area[0] = int((3*tol)/delta[0])
        check_area[1] = int((3*tol)/delta[1])
        check_area[2] = int((3*tol)/delta[2])

        limits = np.array( [
            [iatm_pos[0]-check_area[0], iatm_pos[0]+check_area[0]],
            [iatm_pos[1]-check_area[1], iatm_pos[1]+check_area[1]],
            [iatm_pos[2]-check_area[2], iatm_pos[2]+check_area[2]]
            ] )

        for i in range(limits[0,0], limits[0,1], 1):
            for j in range(limits[1,0], limits[1,1], 1):
                for k in range(limits[2,0], limits[2,1], 1):   

                    r = index2pos([i,j,k],nDim,lvec)
                    vector = r - atom_pos
                    dist = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

                    if (dist < tol):
                        if list_points is None:
                            #print "initiate" , dist, i, j ,k
                            list_points = np.zeros([1,3], dtype=int)
                            list_dist = np.zeros([1], dtype=float)
                            list_points[0,:] = [i,j,k]
                            list_dist[0] = dist
                        else:
                            #print "append" , dist, i, j ,k
                            list_points_aux = np.zeros([1,3], dtype=int)
                            list_dist_aux = np.zeros([1], dtype=float)
                            list_points_aux[0,:] = [i,j,k]
                            list_dist_aux[0] = dist
                            list_points=np.append(list_points, list_points_aux, axis=0)
                            list_dist=np.append(list_dist, list_dist_aux, axis=0)    
        i_min = iatm_pos[0]
        i_max = iatm_pos[0]
        j_min = iatm_pos[1]
        j_max = iatm_pos[1]     
        k_min = iatm_pos[2]
        k_max = iatm_pos[2]

        for i in range(list_dist.shape[0]):
            
            #find the ring edges
            if (list_points[i,0]>i_max): i_max = list_points[i,0]
            if (list_points[i,0]<i_min): i_min = list_points[i,0]
            if (list_points[i,1]>j_max): j_max = list_points[i,1]
            if (list_points[i,1]<j_min): j_min = list_points[i,1]
            if (list_points[i,2]>k_max): k_max = list_points[i,2]
            if (list_points[i,2]<k_min): k_min = list_points[i,2]      

        i_max_pbc = pbc(i_max,nDim[2])
        i_min_pbc = pbc(i_min,nDim[2])
        j_max_pbc = pbc(j_max,nDim[1])
        j_min_pbc = pbc(j_min,nDim[1])
        k_max_pbc = pbc(k_max,nDim[0])
        k_min_pbc = pbc(k_min,nDim[0])

        iatm_pos_pbc = iatm_pos
        iatm_pos_pbc[0] = pbc(iatm_pos_pbc[0],nDim[2])
        iatm_pos_pbc[1] = pbc(iatm_pos_pbc[1],nDim[1])
        iatm_pos_pbc[2] = pbc(iatm_pos_pbc[2],nDim[0])

        #media_ring_x = ( V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_min_pbc] + V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_max_pbc] )/2.0
        #media_ring_y = ( V[iatm_pos_pbc[2],j_min_pbc,iatm_pos_pbc[0]] + V[iatm_pos_pbc[2],j_max_pbc,iatm_pos_pbc[0]] )/2.0
        #media_ring_z = ( V[k_min_pbc,iatm_pos_pbc[1],iatm_pos_pbc[0]] + V[k_max_pbc,iatm_pos_pbc[2],iatm_pos_pbc[0]] )/2.0
        Centro = (media_ring_x + media_ring_y + media_ring_z)/3.0
        print "media and values x", media_ring_x , V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_min_pbc] , V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_max_pbc] 
        print "media and values y", media_ring_y , V[iatm_pos_pbc[2],j_min_pbc,iatm_pos_pbc[0]] , V[iatm_pos_pbc[2],j_max_pbc,iatm_pos_pbc[0]]
        print "media and values z", media_ring_z , V[k_min_pbc,iatm_pos_pbc[1],iatm_pos_pbc[0]] , V[k_max_pbc,iatm_pos_pbc[2],iatm_pos_pbc[0]]

        data2fit_x = np.array([
            index2pos( [iatm_pos[2],iatm_pos[1],i_min] ,nDim,lvec),
            index2pos( [iatm_pos[2],iatm_pos[1],i_max] ,nDim,lvec),
            index2pos( [iatm_pos[2],j_min,iatm_pos[0]] ,nDim,lvec),
            index2pos( [iatm_pos[2],j_max,iatm_pos[0]] ,nDim,lvec),
            index2pos( [k_min,iatm_pos[1],iatm_pos[0]] ,nDim,lvec),
            index2pos( [k_max,iatm_pos[2],iatm_pos[0]] ,nDim,lvec),
            ] )
        
        #print "data2fit_x, data2fit_y",   data2fit_x


        data2fit_y = np.array([
            V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_min_pbc],
            V[iatm_pos_pbc[2],iatm_pos_pbc[1],i_max_pbc],
            V[iatm_pos_pbc[2],j_min_pbc,iatm_pos_pbc[0]],
            V[iatm_pos_pbc[2],j_max_pbc,iatm_pos_pbc[0]],
            V[k_min_pbc,iatm_pos_pbc[1],iatm_pos_pbc[0]],
            V[k_max_pbc,iatm_pos_pbc[2],iatm_pos_pbc[0]],
            ])

        #print "data2fit_x, data2fit_y",   data2fit_x, data2fit_y
        for i in range(list_dist.shape[0]):
            V_corrected[ 
                pbc(list_points[i,2],nDim[0]), 
                pbc(list_points[i,1],nDim[1]), 
                pbc(list_points[i,0],nDim[2]) ] = TD_gaussian_fit(atom_pos, data2fit_x, data2fit_y, Centro, [list_points[i,0], list_points[i,1], list_points[i,2]])
            print "fitted value" , TD_gaussian_fit(atom_pos, data2fit_x, data2fit_y, Centro, [list_points[i,0], list_points[i,1], list_points[i,2]])

    GU.saveXSF( "Corrected_tip.xsf", V_corrected, lvec )
    return V_corrected


def core_correction_old(V,lvec,nDim,Rs_tip,elems_tip,tol):
    
    V_corrected = V.copy()
    V_corrected_zero = V.copy()
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
        check_area = np.zeros([3], dtype=int)
        check_area[0] = int((3*tol)/delta[0])
        check_area[1] = int((3*tol)/delta[1])
        check_area[2] = int((3*tol)/delta[2])
        print "check area" , check_area
#        for i in range(0, nDim[2], 1):
#            for j in range(0, nDim[1], 1):
#                for k in range(0, nDim[0], 1):

        limits = np.array( [
            [iatm_pos[0]-check_area[0], iatm_pos[0]+check_area[0]],
            [iatm_pos[1]-check_area[1], iatm_pos[1]+check_area[1]],
            [iatm_pos[2]-check_area[2], iatm_pos[2]+check_area[2]]
            ] )
        print 'limist!!!!!!!!!!!!!!!111', limits
        
        if (limits[0,1] > nDim[2]): limits[0,1] = nDim[2]
        if (limits[1,1] > nDim[1]): limits[1,1] = nDim[1]
        if (limits[2,1] > nDim[0]): limits[2,1] = nDim[0]

        

        for i in range(limits[0,0], limits[0,1], 1):
            for j in range(limits[1,0], limits[1,1], 1):
                for k in range(limits[2,0], limits[2,1], 1):                    

                    r = index2pos([i,j,k],nDim,lvec)
                    vector = r - atom_pos
                    dist = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
                    #print i,j,k,dist
                    if (dist < tol):
                        if list_points is None:
                            #print "initiate" , dist, i, j ,k
                            list_points = np.zeros([1,3], dtype=int)
                            list_dist = np.zeros([1], dtype=float)
                            list_points[0,:] = [i,j,k]
                            list_dist[0] = dist
                        else:
                            #print "append" , dist, i, j ,k
                            list_points_aux = np.zeros([1,3], dtype=int)
                            list_dist_aux = np.zeros([1], dtype=float)
                            list_points_aux[0,:] = [i,j,k]
                            list_dist_aux[0] = dist
                            list_points=np.append(list_points, list_points_aux, axis=0)
                            list_dist=np.append(list_dist, list_dist_aux, axis=0)
        #print "list_points"
        #print list_points
        #print "list_dist"
        #print list_dist

        media = 0.0
        Norm = 0.0

        #i_min = iatm_pos[0]
        #i_max = iatm_pos[0]
        #j_min = iatm_pos[1]
        #j_max = iatm_pos[1]     
        #k_min = iatm_pos[2]
        #k_max = iatm_pos[2]

        for i in range(list_dist.shape[0]):
            if (list_dist[i] > 0.5*tol):
            #if ( abs(V[ list_points[i,2], list_points[i,1], list_points[i,0] ]) < 2.0 ):
                #print "at k j i", list_points[i,2], list_points[i,1], list_points[i,0] , "V is", V[ list_points[i,2], list_points[i,1], list_points[i,0] ]
                media = media + (list_dist[i]**4)*V[ list_points[i,2], list_points[i,1], list_points[i,0] ]
                Norm = Norm + list_dist[i]**4
                #print "i and partial media and norm",  i, media, Norm, media/Norm

            #if (list_points[i,0]>i_max): i_max = list_points[i,0]
            #if (list_points[i,0]<i_min): i_min = list_points[i,0]
            #if (list_points[i,1]>j_max): j_max = list_points[i,1]
            #if (list_points[i,1]<j_min): j_min = list_points[i,1]
            #if (list_points[i,2]>k_max): k_max = list_points[i,2]
            #if (list_points[i,2]<k_min): k_min = list_points[i,2]           

        media = media/Norm
        #media_ring = ( V[i_min,iatm_pos[1],iatm_pos[2]] + V[i_max,iatm_pos[1],iatm_pos[2]] + V[iatm_pos[0],j_min,iatm_pos[2]] + V[iatm_pos[0],j_max,iatm_pos[2]] + V[iatm_pos[0],iatm_pos[1],k_min] + V[iatm_pos[0],iatm_pos[2],k_max] )/6.0
        #media_ring = ( 
        #    V[iatm_pos[2],iatm_pos[1],i_min] + V[iatm_pos[2],iatm_pos[1],i_max] + 
        #    V[iatm_pos[2],j_min,iatm_pos[0]] + V[iatm_pos[2],j_max,iatm_pos[0]] + 
        #    V[k_min,iatm_pos[1],iatm_pos[0]] + V[k_max,iatm_pos[2],iatm_pos[0]] )/6.0

        for i in range(list_dist.shape[0]):
            
            r = index2pos([list_points[i,0], list_points[i,1], list_points[i,2]],nDim,lvec)
        #    V_corrected = media + (media_ring - media)*((r/tol)**2)
            #V_corrected[list_points[i,2], list_points[i,1], list_points[i,0]] = media + media*((tol - dist)/tol)*0.1
            V_corrected[list_points[i,2], list_points[i,1], list_points[i,0]] = media
            #print "making point" , list_points[i,2], list_points[i,1], list_points[i,0] , "zero"
            V_corrected_zero[list_points[i,2], list_points[i,1], list_points[i,0]] = 0.0

        print iatm_pos, atom_pos

        print "MEDIA", media

        #print "list point and distances"
        #print list_points
        #print list_dist
    #V_corrected[V_corrected < -0.5] = -0.5
    #V_corrected_zero[V_corrected_zero < -0.5] = -0.5

    print "saving Corrected tip"
    np.savetxt('profile', np.array(V[:,0,0]) )
    GU.saveXSF( "Corrected_tip.xsf", V_corrected, lvec )
    np.savetxt('profile_corrected', np.array(V_corrected[:,0,0]) )
    GU.saveXSF( "Corrected_tip_zero.xsf", V_corrected_zero, lvec )
    np.savetxt('profile_corrected_zero', np.array(V_corrected_zero[:,0,0]) )

    return V_corrected


