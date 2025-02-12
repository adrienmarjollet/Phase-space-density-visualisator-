'''
This module arranges the data so that it can be treated with the function store_coordinates
'''
import os
import shutil
import logging

class ArangeFile:

    def __init__(self, SP, Energy_col, b_max, NUMBER_TRAJ, reaction='MuH2v1'):
        self.SP = SP
        self.Energy_col = Energy_col
        self.b_max = b_max
        self.NUMBER_TRAJ = NUMBER_TRAJ
        self.reaction = reaction
        self.setup_logging()

    def setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    def create_directories(self, current_directory, P):
        for p in range(P):
            final_directory = os.path.join(current_directory, "R_Ecol_" + str(self.Energy_col[p]))
            if not os.path.exists(final_directory):
                os.makedirs(final_directory)
                logging.info(f"Created directory: {final_directory}")

    def move_reactive_trajectories(self, current_directory, P):
        count_Ecol_R = [0.0 for _ in range(P)]
        for p in range(P):
            for k in range(1, self.SP + 1):
                p_i_k_dir = os.path.join(current_directory, f"Ecol_{self.Energy_col[p]}", f"b{self.b_max}", f"job{k}")
                if os.path.exists(p_i_k_dir):
                    for tj in range(self.NUMBER_TRAJ):
                        reactive_file = os.path.join(p_i_k_dir, f"reactive_{tj}.xyz")
                        if os.path.exists(reactive_file):
                            with open(reactive_file, 'r') as xyz_file:
                                lines = xyz_file.readlines()
                            if len(lines) > 1:
                                count_Ecol_R[p] += 1.0
                                destination = os.path.join(current_directory, f"R_Ecol_{self.Energy_col[p]}", f"num{count_Ecol_R[p]}_{self.reaction}.xyz")
                                shutil.move(reactive_file, destination)
                                logging.info(f"Moved {reactive_file} to {destination}")

    def select_store_data(self, check_Path=True, stored=False):
        '''
        CREATE FOLDER TO STORE THE XYZ FILES FOR THE REACTIVE TRAJECTORIES 
        '''
        try:
            current_directory = os.getcwd()
            P = len(self.Energy_col)
            if check_Path and not stored:
                self.create_directories(current_directory, P)
                stored = True

            if stored:
                self.move_reactive_trajectories(current_directory, P)
        except Exception as e:
            logging.error(f"An error occurred: {e}")
