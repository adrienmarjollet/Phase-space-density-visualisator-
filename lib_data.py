'''
This module gather the selected trajectories (replicas, classical and centroids).
'''
from units import kcal2cm, femtoSec2au,kcal2cm,au2cm,ang2au,eV2cm,mass_n_au,mass_p_au,mass_mu_au,mass_e_au,mH,mMu

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

import lib_coordinates


def store_coordinates(R_range,r_range,number_trajs,N,D,n,m,Energy_col,excitation,number_frames,start_frame):
    PLOT_DYNAMIC_2D_HISTO_CLASSICAL_COMBINED = 1
    PLOT_DYNAMIC_2D_HISTO_REPLICAS_COMBINED = 1

    PLOT_STATIC_2D_HISTO_CENTROID = 0
    PLOT_STATIC_2D_HISTO_REPLICAS = 0

    PLOT_DYNAMIC_2D_HISTO_REPLICAS_SINGLES = 0

    PLOT_STATIC_2D_HISTO_CENTROID = 0
    PLOT_STATIC_2D_HISTO_REPLICAS = 0

    PLOT_INSTANTON_TRAJECTORY = False
    PLOT_LINES_BETWEEN_POINTS = False
    '''                  
    PLOT OPTIONS
    '''

    def define_FramePerSec(R, R_array, fps_array):
        '''
        Used later to enhance the number of frame per seconds around the collision moment.
        '''
        ind_array = len(R_array)
        for i in range(len(R_array)-1):
            if R > R_array[i] and R <= R_array[i+1]:
                ind_array = i
        return fps_array[ind_array]

    '''
    Storing the coordinates for the subsequent density plot.
    '''
    R_add_c, r_add_c, R_add_repli, r_add_repli  = np.array([]), np.array([]), np.array([]), np.array([])
    lst_replicas, P  = [int(num) for num in range(n)], len(Energy_col)

    current_directory = os.getcwd()

    for p in range(len(Energy_col)):

        loc_dir = current_directory+"/R_Ecol_"+str(Energy_col[p])
        count_xyz = 1
        '''
        COMBINED HIST2D
        '''
        list_frames_trajs_R_cl, list_frames_trajs_r_cl  = [[] for j in range(number_frames)], [[] for j in range(number_frames)]
        list_frames_trajs_R, list_frames_trajs_r = [[] for j in range(number_frames)], [[] for j in range(number_frames)]

        while os.path.exists(loc_dir+"/num"+str(float(count_xyz))+"_MuH2v"+str(int(excitation))+".xyz") and (count_xyz<number_trajs):

            print("traj number:", count_xyz)
            R_centro = np.array([])
            r_centro = np.array([])
            with open(loc_dir+"/num"+str(float(count_xyz))+"_MuH2v"+str(int(excitation))+".xyz", 'r') as xyz_file:
                count_xyz += 1
                lines = xyz_file.readlines()
                xyz_file.close()
                Ntot, num_start_cl, delta_start = 0, 0, 0
                num_start_replicas = np.array([(i+2) for i in range(N*n)])
                Ntot = int(lines[0])
                for num_line, line in enumerate(lines):
                    pp = line.split()
                    pp = np.asarray(pp)
                    if len(pp) > 2 and pp[0] == "M":
                        num_start_cl = num_line
                    if len(pp) == 1 and num_line > 0 and pp[0] == str(Ntot):
                        delta_start = num_line
                    if num_line > Ntot+3:
                        break
                number_lines = len(lines)
                R_succ_c, r_succ_c = np.array([None]), np.array([None])
                R_centro, r_centro = np.array([]), np.array([])
                for ii in range(0, number_lines, delta_start):
                    if ii >= (start_frame*delta_start) and ii < ((number_frames)*delta_start):
                        ii_centro = ii + 3 # locate to the centroid coordinate 
                        Zq = np.zeros(N*D*n)
                        i = 1
                        for k, line in enumerate(lines[ii_centro:ii_centro+n]):
                            pp = line.split()
                            pp = np.asarray(pp)
                            if pp[0] == 'H':
                                for j in range(D):
                                    Zq[i*D*n+j*n+k] = float(pp[j+1])*ang2au
                            else:
                                print('centro problem with H1')
                                sys.exit()
                        i = 2
                        for k, line in enumerate(lines[ii_centro+n-1:ii_centro+2*n-1]):
                            pp = line.split()
                            pp = np.asarray(pp)
                            if pp[0] == 'H':
                                for j in range(D):
                                    Zq[i*D*n+j*n+k] = float(pp[j+1])*ang2au
                            else:
                                print('centro problem with H2')
                                print(pp[0])
                                sys.exit()
                        i = 0
                        for k, line in enumerate(lines[ii_centro+2*n-1:ii_centro+3*n-1]):
                            pp = line.split()
                            pp = np.asarray(pp)
                            if pp[0] == 'Mu' or pp[0] == 'D':
                                for j in range(D):
                                    Zq[i*D*n+j*n+k] = float(pp[j+1])*ang2au
                            else:
                                print('centro problem with Mu')
                                print(pp[0])
                                sys.exit()
                        ZKq_temp = lib_coordinates.centroids_q(Zq)
                        R_next_centro, r_next_centro = lib_coordinates.xyzcoords_to_Rr(m, ZKq_temp)
                        if PLOT_STATIC_2D_HISTO_CENTROID == 1:
                            ZKq_temp = lib_coordinates.centroids_q(Zq)
                            R_next_centro, r_next_centro = lib_coordinates.xyzcoords_to_Rr(
                                m, ZKq_temp)
                            R_add_c = np.append(R_add_c, R_next_centro)
                            r_add_c = np.append(r_add_c, r_next_centro)

                        if PLOT_STATIC_2D_HISTO_REPLICAS == 1:
                            for k in lst_replicas:
                                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(
                                    N, D, n, Zq, k, m)
                                R_add_repli = np.append(R_add_repli, R_next_repli)
                                r_add_repli = np.append(r_add_repli, r_next_repli)

                        if PLOT_DYNAMIC_2D_HISTO_REPLICAS_SINGLES == 1:
                            '''
                            #TODO: THIS ONE IS VERY UNOPTIMIZE
                            '''
                            R_dyn_repli = np.array([])
                            r_dyn_repli = np.array([])
                            for k in lst_replicas:
                                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(
                                    N, D, n, Zq, k, m)
                                R_dyn_repli = np.append(R_dyn_repli, R_next_repli)
                                r_dyn_repli = np.append(r_dyn_repli, r_next_repli)
                            plt.hist2d(R_dyn_repli, r_dyn_repli, bins=[90, 70], cmap='Blues', range=[
                                    [min(R_range), max(R_range)], [min(r_range), max(r_range)]])  # ,cmin=1,normed=0)
                            # plt.pause(0.02)
                            # WE ADAPT THE FPS TO EMPHASIZE AROUND THE TRANSITION STATE
                            plt.pause(1.0/define_FramePerSec(R_next_centro, R_range, R_range*3))

                        '''
                        COMBINED HIST2D
                        '''
                        if PLOT_DYNAMIC_2D_HISTO_REPLICAS_COMBINED == 1:
                            R_dyn_repli = np.array([])
                            r_dyn_repli = np.array([])
                            for k in lst_replicas:
                                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(
                                    N, D, n, Zq, k, m)
                                R_dyn_repli = np.append(R_dyn_repli, R_next_repli)
                                r_dyn_repli = np.append(r_dyn_repli, r_next_repli)
                            #print('frame :',ii//delta_start)
                            list_frames_trajs_R[ii //
                                                delta_start].extend(list(R_dyn_repli))
                            list_frames_trajs_r[ii //
                                                delta_start].extend(list(r_dyn_repli))

                        '''
                        Classical part
                        '''

                        ii_cl = ii + num_start_cl

                        q_cl = np.zeros(N*D)

                        for line in lines[ii_cl:ii_cl+3]:
                            pp = line.split()
                            pp = np.asarray(pp)
                            if pp[0] == 'P':
                                i = 0
                                pp = line.split()
                                pp = np.asarray(pp)
                                for j in range(D):
                                    q_cl[i*D+j] = float(pp[j+1])*ang2au

                            elif pp[0] == 'M':
                                i = 1
                                pp = line.split()
                                pp = np.asarray(pp)
                                for j in range(D):
                                    q_cl[i*D+j] = float(pp[j+1])*ang2au

                            elif pp[0] == 'L':
                                i = 2
                                pp = line.split()
                                pp = np.asarray(pp)
                                for j in range(D):
                                    q_cl[i*D+j] = float(pp[j+1])*ang2au

                        R_next_cl, r_next_cl = lib_coordinates.xyzcoords_to_Rr(m, q_cl)

                        if PLOT_DYNAMIC_2D_HISTO_CLASSICAL_COMBINED == 1:
                            R_dyn_cl = np.array([])
                            r_dyn_cl = np.array([])

                            R_dyn_cl = np.append(R_dyn_cl, R_next_cl)
                            r_dyn_cl = np.append(r_dyn_cl, r_next_cl)

                            list_frames_trajs_R_cl[ii //
                                                delta_start].extend(list(R_dyn_cl))
                            list_frames_trajs_r_cl[ii //
                                                delta_start].extend(list(r_dyn_cl))

    return R_add_c, r_add_c, R_add_repli, r_add_repli, list_frames_trajs_R_cl, list_frames_trajs_r_cl,list_frames_trajs_R, list_frames_trajs_r
