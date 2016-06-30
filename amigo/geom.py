import numpy as np
from scipy.interpolate import interp1d

def length(R,Z,norm=True):
    L = np.append(0,np.cumsum(np.sqrt(np.diff(R)**2+np.diff(Z)**2)))
    if norm: L = L/L[-1]
    return L
    
def vector_length(X,norm=True):
    L = np.append(0,np.cumsum(np.sqrt(np.diff(X[:,0])**2+
                                      np.diff(X[:,1])**2+
                                      np.diff(X[:,2])**2)))
    if norm: L = L/L[-1]
    return L
    
def space(R,Z,npoints):
    L = length(R,Z)
    l = np.linspace(0,1,npoints)
    R,Z = interp1d(L,R)(l),interp1d(L,Z)(l)
    return R,Z
    