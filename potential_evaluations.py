'''
Library to evaluate potential energy functions for collisional triatomic systems 
with certain insightful constraints.
'''
import numpy as np
import os
import pot
import sys

class pot_eval:
    '''
    Default potential: .so file which is an accurate compiled potential program for a system consisting of 3 hydrogen atoms ('H3')
    '''
    def __init__(self,N,D,which_pot='H3'): 
        self.N = N
        self.D = D
        self.which_pot = which_pot
        if self.which_pot=='H3' and os.path.isfile('./pot.cpython-36m-x86_64-linux-gnu.so')
            pot.pes_init()
            print('H3 Potential initialized.')
        elif os.path.exists('./pot.py'):
            print('Molecular potential imported.')
        else: 
            print('No molecular potential to work with.')
            sys.exit()
            
    @staticmethod
    def POT_H3(position_matrix):
        '''
        Potential energy as function of the Cartesian coordinates
        Indices [i,j] corresponds to ith atom and jth dimension
        Return the potential in atomic units (Hartree)
        '''
        v, dummy = pot.pot(position_matrix)
        v += 0.174496  #To set the zero of the potential (does not influence the dynamics)
        return v
    
    @staticmethod
    def POT_H3_rRa(r, R, angle):
        '''
        We use cylindric symmetry of the target when its cm is on the main axis 
        '''
        qq_madeup = np.zeros((self.N, self.D))
        qq_madeup[1, 0] = -R + 0.5*r*np.cos(angle)
        qq_madeup[1, 1] = 0.5*r*np.sin(angle)
        qq_madeup[2, 0] = -R - 0.5*r*np.cos(angle)
        qq_madeup[2, 1] = -0.5*r*np.sin(angle)
        v, dummy = pot.pot(qq_madeup)
        v += 0.174496
        return v

    def POT_H3_rR_min_BruteForced(self,R_arr, r_arr):
        '''
        BF algo (used once per compilation) to find the global minimums for a set of distances.
        R_arr : distances for R
        r_arr : distance for r

        Returns a matrix V with rows ~ r and columns ~ R
        '''
        dang = (1.0/50)*(np.pi/2) #angular step of 0.02*pi/2 radian
        V = np.zeros((len(r_arr), len(R_arr)))
        for i in range(len(r_arr)):
            for j in range(len(R_arr)):
                min_v = 1.0
                for k in range(int(round(np.pi/2/dang))):
                    x = np.array([float(0) for i in range(self.D)])
                    x[0], x[1], x[2] = R_arr[j], r_arr[i], dang*k
                    v_angle = POT_H3_rRa(x)
                    if v_angle < min_v:
                        min_v = v_angle
                V[i, j] = min_v
        return V
