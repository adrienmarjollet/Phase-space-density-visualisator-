'''
This module gathers the selected trajectories (replicas, classical and centroids).
'''

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from units import kcal2cm, femtoSec2au, kcal2cm, au2cm, ang2au, eV2cm, mass_n_au, mass_p_au, mass_mu_au, mass_e_au, mH, mMu
import lib_coordinates

class TrajectoryData:
    
    def __init__(self, R_range, r_range, number_trajs, N, D, n, m, Energy_col, excitation, number_frames, start_frame):
        self.R_range = R_range
        self.r_range = r_range
        self.number_trajs = number_trajs
        self.N = N
        self.D = D
        self.n = n
        self.m = m
        self.Energy_col = Energy_col
        self.excitation = excitation
        self.number_frames = number_frames
        self.start_frame = start_frame
        self.current_directory = os.getcwd()
        self.R_add_c = np.array([])
        self.r_add_c = np.array([])
        self.R_add_repli = np.array([])
        self.r_add_repli = np.array([])
        self.list_frames_trajs_R_cl = [[] for _ in range(number_frames)]
        self.list_frames_trajs_r_cl = [[] for _ in range(number_frames)]
        self.list_frames_trajs_R = [[] for _ in range(number_frames)]
        self.list_frames_trajs_r = [[] for _ in range(number_frames)]
        self.lst_replicas = [int(num) for num in range(n)]
        self.P = len(Energy_col)

    def define_FramePerSec(self, R, R_array, fps_array):
        '''
        Used later to enhance the number of frame per seconds around the collision moment.
        '''
        ind_array = len(R_array)
        for i in range(len(R_array) - 1):
            if R > R_array[i] and R <= R_array[i + 1]:
                ind_array = i
        return fps_array[ind_array]

    def process_trajectory(self, loc_dir, count_xyz):
        '''
        Process a single trajectory file.
        '''
        R_centro = np.array([])
        r_centro = np.array([])
        with open(loc_dir + f"/num{float(count_xyz)}_MuH2v{int(self.excitation)}.xyz", 'r') as xyz_file:
            lines = xyz_file.readlines()
        Ntot, num_start_cl, delta_start = 0, 0, 0
        num_start_replicas = np.array([(i + 2) for i in range(self.N * self.n)])
        Ntot = int(lines[0])
        for num_line, line in enumerate(lines):
            pp = line.split()
            pp = np.asarray(pp)
            if len(pp) > 2 and pp[0] == "M":
                num_start_cl = num_line
            if len(pp) == 1 and num_line > 0 and pp[0] == str(Ntot):
                delta_start = num_line
            if num_line > Ntot + 3:
                break
        number_lines = len(lines)
        for ii in range(0, number_lines, delta_start):
            if ii >= (self.start_frame * delta_start) and ii < ((self.number_frames) * delta_start):
                self.process_frame(ii, lines, delta_start, num_start_cl)

    def process_frame(self, ii, lines, delta_start, num_start_cl):
        '''
        Process a single frame of a trajectory.
        '''
        ii_centro = ii + 3  # locate to the centroid coordinate
        Zq = np.zeros(self.N * self.D * self.n)
        self.extract_centroid_coordinates(ii_centro, lines, Zq)
        ZKq_temp = lib_coordinates.centroids(Zq, self.n)
        R_next_centro, r_next_centro = lib_coordinates.xyzcoords_to_Rr(self.N, self.D, self.m, ZKq_temp)
        self.update_histograms(R_next_centro, r_next_centro, Zq)
        self.process_classical_part(ii, lines, delta_start, num_start_cl)

    def extract_centroid_coordinates(self, ii_centro, lines, Zq):
        '''
        Extract centroid coordinates from the lines.
        '''
        i = 1
        for k, line in enumerate(lines[ii_centro:ii_centro + self.n]):
            pp = line.split()
            pp = np.asarray(pp)
            if pp[0] == 'H':
                for j in range(self.D):
                    Zq[i * self.D * self.n + j * self.n + k] = float(pp[j + 1]) * ang2au
            else:
                print('centro problem with H1')
                sys.exit()
        i = 2
        for k, line in enumerate(lines[ii_centro + self.n - 1:ii_centro + 2 * self.n - 1]):
            pp = line.split()
            pp = np.asarray(pp)
            if pp[0] == 'H':
                for j in range(self.D):
                    Zq[i * self.D * self.n + j * self.n + k] = float(pp[j + 1]) * ang2au
            else:
                print('centro problem with H2')
                print(pp[0])
                sys.exit()
        i = 0
        for k, line in enumerate(lines[ii_centro + 2 * self.n - 1:ii_centro + 3 * self.n - 1]):
            pp = line.split()
            pp = np.asarray(pp)
            if pp[0] == 'Mu' or pp[0] == 'D':
                for j in range(self.D):
                    Zq[i * self.D * self.n + j * self.n + k] = float(pp[j + 1]) * ang2au
            else:
                print('centro problem with Mu')
                print(pp[0])
                sys.exit()

    def update_histograms(self, R_next_centro, r_next_centro, Zq):
        '''
        Update histograms with the new data.
        '''
        if PLOT_STATIC_2D_HISTO_CENTROID == 1:
            self.R_add_c = np.append(self.R_add_c, R_next_centro)
            self.r_add_c = np.append(self.r_add_c, r_next_centro)

        if PLOT_STATIC_2D_HISTO_REPLICAS == 1:
            for k in self.lst_replicas:
                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(self.N, self.D, self.n, Zq, k, self.m)
                self.R_add_repli = np.append(self.R_add_repli, R_next_repli)
                self.r_add_repli = np.append(self.r_add_repli, r_next_repli)

        if PLOT_DYNAMIC_2D_HISTO_REPLICAS_SINGLES == 1:
            R_dyn_repli = np.array([])
            r_dyn_repli = np.array([])
            for k in self.lst_replicas:
                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(self.N, self.D, self.n, Zq, k, self.m)
                R_dyn_repli = np.append(R_dyn_repli, R_next_repli)
                r_dyn_repli = np.append(r_dyn_repli, r_next_repli)
            plt.hist2d(R_dyn_repli, r_dyn_repli, bins=[90, 70], cmap='Blues', range=[
                [min(self.R_range), max(self.R_range)], [min(self.r_range), max(self.r_range)]])
            plt.pause(1.0 / self.define_FramePerSec(R_next_centro, self.R_range, self.R_range * 3))

        if PLOT_DYNAMIC_2D_HISTO_REPLICAS_COMBINED == 1:
            R_dyn_repli = np.array([])
            r_dyn_repli = np.array([])
            for k in self.lst_replicas:
                R_next_repli, r_next_repli = lib_coordinates.replica_xyzcoords_to_Rr(self.N, self.D, self.n, Zq, k, self.m)
                R_dyn_repli = np.append(R_dyn_repli, R_next_repli)
                r_dyn_repli = np.append(r_dyn_repli, r_next_repli)
            self.list_frames_trajs_R[ii // delta_start].extend(list(R_dyn_repli))
            self.list_frames_trajs_r[ii // delta_start].extend(list(r_dyn_repli))

    def process_classical_part(self, ii, lines, delta_start, num_start_cl):
        '''
        Process the classical part of the trajectory.
        '''
        ii_cl = ii + num_start_cl
        q_cl = np.zeros(self.N * self.D)
        for line in lines[ii_cl:ii_cl + 3]:
            pp = line.split()
            pp = np.asarray(pp)
            if pp[0] == 'P':
                i = 0
                for j in range(self.D):
                    q_cl[i * self.D + j] = float(pp[j + 1]) * ang2au
            elif pp[0] == 'M':
                i = 1
                for j in range(self.D):
                    q_cl[i * self.D + j] = float(pp[j + 1]) * ang2au
            elif pp[0] == 'L':
                i = 2
                for j in range(self.D):
                    q_cl[i * self.D + j] = float(pp[j + 1]) * ang2au
        R_next_cl, r_next_cl = lib_coordinates.xyzcoords_to_Rr(self.N, self.D, self.m, q_cl)
        if PLOT_DYNAMIC_2D_HISTO_CLASSICAL_COMBINED == 1:
            R_dyn_cl = np.array([])
            r_dyn_cl = np.array([])
            R_dyn_cl = np.append(R_dyn_cl, R_next_cl)
            r_dyn_cl = np.append(r_dyn_cl, r_next_cl)
            self.list_frames_trajs_R_cl[ii // delta_start].extend(list(R_dyn_cl))
            self.list_frames_trajs_r_cl[ii // delta_start].extend(list(r_dyn_cl))

    def store_coordinates(self):
        '''
        Main function to store the coordinates for the subsequent density plot.
        '''
        for p in range(self.P):
            loc_dir = os.path.join(self.current_directory, f"R_Ecol_{self.Energy_col[p]}")
            count_xyz = 1
            while os.path.exists(loc_dir + f"/num{float(count_xyz)}_MuH2v{int(self.excitation)}.xyz") and (count_xyz < self.number_trajs):
                print("traj number:", count_xyz)
                self.process_trajectory(loc_dir, count_xyz)
                count_xyz += 1
        return self.R_add_c, self.r_add_c, self.R_add_repli, self.r_add_repli, self.list_frames_trajs_R_cl, self.list_frames_trajs_r_cl, self.list_frames_trajs_R, self.list_frames_trajs_r
