import numpy as np


class TriatomicCongifuration:

    def __init__(self, N, D, m):
        self.N = N
        self.D = D
        self.m = m

    def Q_coords(self, q):
        '''
        Notations are cumbersome but help to see the physical meaning of Q.
        q: position array of size N*D
        Returns
        Q: generalized coordinates (handy coordinates for triatomic systems)
        '''
        try:
            Q = np.zeros(self.N * self.D)
            Q[(self.N-3)*self.D:(self.N-2)*self.D] = q[(self.N-1)*self.D:self.N*self.D] - q[(self.N-2)*self.D:(self.N-1)*self.D]
            Q[(self.N-2)*self.D:(self.N-1)*self.D] = q[(self.N-3)*self.D:self.D] - (self.m[(self.N-2)*self.D:(self.N-1)*self.D]*q[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D]*q[(self.N-1)*self.D:self.N*self.D]) / (self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D])
            Q[(self.N-1)*self.D:self.N*self.D] = (self.m[(self.N-3)*self.D:(self.N-2)*self.D]*q[(self.N-3)*self.D:(self.N-2)*self.D] + self.m[(self.N-2)*self.D:(self.N-1)*self.D]*q[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D]*q[2*self.D:self.N*self.D]) / (self.m[(self.N-3)*self.D:(self.N-2)*self.D] + self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D])
            return Q
        except Exception as e:
            raise ValueError(f"Error in Q_coords: {e}")

    def q_coords(self, Q):
        '''
        Reverse transformation back to displacement (Cartesian) coordinates.
        '''
        try:
            qb = np.zeros(self.N * self.D)
            qb[(self.N-3)*self.D:(self.N-2)*self.D] = ((self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:3*self.D]) / (self.m[(self.N-3)*self.D:self.D] + self.m[self.D:2*self.D] + self.m[2*self.D:3*self.D])) * Q[(self.N-2)*self.D:(self.N-1)*self.D] + Q[(self.N-1)*self.D:self.N*self.D]
            qb[(self.N-2)*self.D:(self.N-1)*self.D] = -((self.m[(self.N-1)*self.D:self.N*self.D]) / (self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D])) * Q[(self.N-3)*self.D:(self.N-2)*self.D] - \
                                  ((self.m[(self.N-3)*self.D:(self.N-2)*self.D]) / (self.m[(self.N-3)*self.D:(self.N-2)*self.D] + self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D])) * Q[(self.N-2)*self.D:(self.N-1)*self.D] + Q[(self.N-2)*self.D:self.N*self.D]
            qb[(self.N-1)*self.D:self.N*self.D] = (self.m[self.D:2*self.D] / (self.m[self.D:2*self.D] + self.m[2*self.D:3*self.D])) * Q[:self.D] - self.m[(self.N-3)*self.D:(self.N-2)*self.D] / \
                              (self.m[(self.N-3)*self.D:(self.N-2)*self.D] + self.m[(self.N-2)*self.D:(self.N-1)*self.D] + self.m[(self.N-1)*self.D:self.N*self.D]) * Q[(self.N-2)*self.D:(self.N-1)*self.D] + Q[(self.N-1)*self.D:self.N*self.D]
            return qb
        except Exception as e:
            raise ValueError(f"Error in q_coords: {e}")

    def xyzcoords_to_Rr(self, q):
        '''
        Input:  Cartesian coordinates q (np.array)
        Output: 
                R: distance between the incoming atom A and BC's center of mass 
        '''
        try:
            Q = self.Q_coords(q)
            return np.sqrt(Q[3]**2 + Q[4]**2 + Q[5]**2), np.sqrt(Q[0]**2 + Q[1]**2 + Q[2]**2)
        except Exception as e:
            raise ValueError(f"Error in xyzcoords_to_Rr: {e}")

    def replica_xyzcoords_to_Rr(self, n, Zq, k):
        '''
        Same as xyzcoords_to_Rr but for each slice of imaginary time, that is to say each ensemble of replica 
        corresponding to the kth slice (among the n total slices).
        '''
        try:
            return self.xyzcoords_to_Rr(Zq[k:self.N*self.D*n+k:n])
        except Exception as e:
            raise ValueError(f"Error in replica_xyzcoords_to_Rr: {e}")

    def centroids(self, Zq, n):
        '''Returns centroid coordinates
           Note: n must be even and greater or equal with 4.     
        '''
        try:
            if n < 4 or n % 2 != 0:
                raise ValueError("n must be even and greater or equal to 4.")
            return np.array([np.mean(Zq[i*n:i*n+n]) for i in range(len(Zq)//n)])
        except Exception as e:
            raise ValueError(f"Error in centroids: {e}")
