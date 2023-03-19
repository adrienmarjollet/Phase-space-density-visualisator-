'''
This module arrange the data so that it can be treated with the function store_coordinates
'''
import os
import shutil

class arange_file:

    def __init__(self,SP,Energy_col,b_max,NUMBER_TRAJ,reaction='MuH2v1'):
        self.SP=SP
        self.Energy_col = Energy_col
        self.b_max = b_max
        self.NUMBER_TRAJ = NUMBER_TRAJ
        self.reaction = reaction


    def select_store_data(self,check_Path=True,stored=False):
        '''
        CREATE FOLDER TO STORE THE XYZ FILES FOR THE REACTIVE TRAJECTORIES 
        '''
        current_directory, P = os.getcwd(), len(self.Energy_col)
        if check_Path and not stored:
            for p in range(P):
                final_directory = os.path.join(
                    current_directory, "R_Ecol_"+str(self.Energy_col[p]))
                if not os.path.exists(final_directory):
                    os.makedirs(final_directory)
                    stored = True
        '''
        We first spot the reactive traj and then put then move them in the created folder (final_directory)
        '''
        if stored == True:
            count_Ecol_R, list_k_R = [float(0) for ii in range(P)], []
            for p in range(P):
                for k in range(1, self.SP+1):
                    p_i_k_dir = current_directory+"/Ecol_" + \
                        str(self.Energy_col[p])+"/b"+str(self.b_max)+"/job"+str(k)
                    if os.path.exists(p_i_k_dir):
                        for tj in range(self.NUMBER_TRAJ):
                            if os.path.exists(p_i_k_dir+"/reactive_"+str(tj)+".xyz"):
                                with open(p_i_k_dir+"/reactive_"+str(tj)+".xyz", 'r') as xyz_file:
                                    lines = xyz_file.readlines()
                                    xyz_file.close()
                                    if len(lines) > 1:
                                        count_Ecol_R[p] += 1.0
                                        shutil.move(p_i_k_dir+"/reactive_"+str(tj)+".xyz", current_directory+"/R_Ecol_"+str(
                                            self.Energy_col[p])+"/num"+str(count_Ecol_R[p])+"_"+self.reaction+".xyz")
