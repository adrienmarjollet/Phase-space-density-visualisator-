import numpy as np 

'''
ONLY USEFUL FOR TRIATOMIC SYSTEMS (N=3) BUT STRAIGHTFORWARDLY GENERALIZABLE TO BIGGER SYSTEMS.

Q : handy generalized coordinates for the reaction between an atom A and dimatomic molecule BC.

q : Cartesian coordinates

i*D+j index refers to the ith atom and the jth spatial dimension.

'''

def Q_coords(N, D, m, q):
    '''
    Notations are cumbersome but help to see the physical meaning of Q.
    q: position array of size N*D
    m: mass array of size N*D
    Returns
    Q: generalized coordinates (handy coordinates for triatomic systems)
    '''
    Q  = np.zeros(N*D)
    Q[(N-3)*D:(N-2)*D] = q[(N-1)*D:N*D] - q[(N-2)*D:(N-1)*D]

    Q[(N-2)*D:(N-1)*D] = q[(N-3)*D:D] - (m[(N-2)*D:(N-1)*D]*q[(N-2)*D:(N-1)*D] + m[(N-1)*D:N*D]*q[(N-1)*D:N*D] ) / ( m[(N-2)*D:(N-1)*D] + m[(N-1)*D:N*D] )

    Q[(N-1)*D:N*D] = (m[(N-3)*D:(N-2)*D]*q[(N-3)*D:(N-2)*D] + m[(N-2)*D:(N-1)*D]*q[(N-2)*D:(N-1)*D] + 
                    m[(N-1)*D:N*D]*q[2*D:N*D])/(m[(N-3)*D:(N-2)*D] + m[(N-2)*D:(N-1)*D] + m[(N-1)*D:N*D])
    return Q


def q_coords(N,D,m,Q):
    '''
    Reverse transformation back to displacement (Cartesian) coordinates.
    '''
    qb = np.zeros(N*D)
    qb[(N-3)*D:(N-2)*D] = ((m[(N-2)*D:(N-1)*D]+m[(N-1)*D:3*D])/(m[(N-3)*D:D]+m[D:2*D]+m[2*D:3*D]))*Q[(N-2)*D:(N-1)*D] + Q[(N-1)*D:N*D]

    qb[(N-2)*D:(N-1)*D] =  - ((m[(N-1)*D:N*D])/(m[(N-2)*D:(N-1)*D]+m[(N-1)*D:N*D]))*Q[(N-3)*D:(N-2)*D] - \
                            ((m[(N-3)*D:(N-2)*D])/(m[(N-3)*D:(N-2)*D]+m[(N-2)*D:(N-1)*D]+m[(N-1)*D:N*D]))*Q[(N-2)*D:(N-1)*D] + Q[(N-2)*D:N*D]
    
    qb[(N-1)*D:N*D]   =    (m[D:2*D]/(m[D:2*D]+m[2*D:3*D]))*Q[:D]   -  m[(N-3)*D:(N-2)*D]   / \
                          (m[(N-3)*D:(N-2)*D]+m[(N-2)*D:(N-1)*D]+m[(N-1)*D:N*D])*Q[(N-2)*D:(N-1)*D] + Q[(N-1)*D:N*D]
    return qb

def xyzcoords_to_Rr(N, D, m, q):
    Q = Q_coords(N,D,m,q)
    return np.sqrt(Q[3]**2 + Q[4]**2 + Q[5]**2), np.sqrt(Q[0]**2 + Q[1]**2 + Q[2]**2)

def replica_xyzcoords_to_Rr(N, D, n, Zq, k, m):
    return xyzcoords_to_Rr(m, Zq[k:N*D*n+k:n])

def centroids(Zq,n):
    '''Returns centroid coordinates
       Note: n must be even and greater or equal with 4.     
    '''
    return np.array([np.mean(Zq[i*n:i*n+n]) for i in range(len(Zq)//n)])
