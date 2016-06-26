import numpy as np
import pylab as pl
from scipy.sparse import lil_matrix,csr_matrix
from scipy.sparse.linalg import spsolve
from collections import OrderedDict

import seaborn as sns
rc = {'figure.figsize':[10,10*12/16],'savefig.dpi':100, # 
      'savefig.jpeg_quality':100,'savefig.pad_inches':0.1,
      'lines.linewidth':2}
sns.set(context='talk',style='white',font='sans-serif',palette='Set2',
        font_scale=7/8,rc=rc)

def delete_row_csr(mat, i):
    if not isinstance(mat, csr_matrix):
        raise ValueError('works only for CSR format -- use .tocsr() first')
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

# element stiffness matrix
class beam(object):
    
    def __init__(self,X,frame='1D',E=1,A=1,Iy=1,Iz=1,G=1,J=1):
        self.E,self.G,self.J = E,G,J
        self.Iy,self.Iz = Iy,Iz
        self.A = A
        self.X = X
        self.frame = frame 
        self.grid()
        self.set_stiffness()
        self.initalize_BC()  # boundary conditions
        self.initalize_L()  # load
        self.initalize_D()  # displacment
        
    def grid(self):
        self.nel = np.shape(self.X)[0]-1  # element number
        dx = np.array(np.diff(self.X,axis=0))  # element aligned coordinate
        self.check_frame(dx)
        dl = np.linalg.norm(dx,axis=1)  # element length
        norm = np.dot(np.matrix(dl).T,np.ones((1,3)))
        dx /= norm
        Rz90 = np.array([[0,-1, 0],[1,0,0],[0,0,1]])  # pi/2 rotation about z
        #Ry90 = np.array([[0, 0,-1],[0,1,0],[1,0,0]])  # -pi/2 rotation about y
        dy = np.dot(Rz90,dx.T).T/norm  # rotate about z
        dz = np.cross(dx,dy)
        #dz = np.dot(Ry90,dx.T).T  # rotate about y
        self.dX,self.dY,self.dZ,self.dL = dx,dy,dz,dl
        
    def initalize_D(self):  # initalize displacment vector [u,v,w,rx,ry,rz]
        self.D = {}
        for disp in ['x','y','z','tx','ty','tz']:
            self.D[disp] = np.zeros(self.nel+1)
            
    def initalize_L(self):
        self.L = np.zeros(self.nK)  # loading vector

    def initalize_BC(self):
        self.BC = OrderedDict()
        self.BC['fix'] = np.array([])# type:[node,...] 
        for key in self.dof:
            self.BC[key] = np.array([])# type:[node,...] 
        
    def set_stiffness(self):  # set element stiffness matrix
        if self.frame=='1D':
            dof = ['v','rz']
            disp = ['y','tz']
            stiffness = self.stiffness_1D
            
        elif self.frame=='2D':
            dof = ['u','v','rz']
            disp = ['x','y','tz']
            stiffness = self.stiffness_2D
        else:
            dof = ['u','v','w','rx','ry','rz']
            disp = ['x','y','z','tz','ty','tz']
            stiffness = self.stiffness_3D  
        self.dof = dof
        self.disp = disp
        self.ndof = int(2*len(dof))  # frame element degrees of fredom
        self.node_dof = int(self.ndof/2)  # ndof per node
        self.nK = int((self.nel-1)*self.ndof/2+self.ndof)
        self.stiffness = stiffness
        return dof
        
    def stiffness_1D(self,dx,dy,dz,dl):  # dof [v,rz]
        a = dl/2
        E,Iz = self.E,self.Iz
        k = E*Iz/(2*a**3)*np.matrix([[3,  3*a,   -3,  3*a],
                                     [3*a,4*a**2,-3*a,2*a**2],
                                     [-3,-3*a,    3, -3*a],
                                     [3*a,2*a**2,-3*a,4*a**2]])    
        return k

    def stiffness_2D(self,dx,dy,dz,dl):  # dof [u,v,rz]
        a = dl/2
        E,Iz,A = self.E,self.Iz,self.A
        k = E*np.matrix([[A/(2*a), 0,             0,
                         -A/(2*a), 0,             0],
                         [0,       3*Iz/(2*a**3), 3*Iz/(2*a**2),
                          0,      -3*Iz/(2*a**3), 3*Iz/(2*a**2)],
                         [0,       3*Iz/(2*a**2), 2*Iz/a,
                          0,      -3*Iz/(2*a**2), Iz/a],
                         [-A/(2*a),0,             0,
                          A/(2*a), 0,             0],
                         [0,      -3*Iz/(2*a**3),-3*Iz/(2*a**2),
                          0,       3*Iz/(2*a**3),-3*Iz/(2*a**2)],
                         [0,       3*Iz/(2*a**2), Iz/a,
                          0,      -3*Iz/(2*a**2), 2*Iz/a]])
        k = self.rotate(k,dx,dy,dz)  # transfer to global coordinates
        return k
        
    def stiffness_3D(self,dx,dy,dz,dl):  # dof [u,v,w,rx,ry,rz]
        a = dl/2
        E,A,G,J = self.E,self.A,self.G,self.J
        Iy,Iz = self.Iy,self.Iz
        k = np.matrix([[A*E/(2*a), 0,               0,
                        0,         0,               0,
                        -A*E/(2*a),0,               0,
                        0,         0,               0],
                       [0,         3*E*Iz/(2*a**3), 0,
                        0,         0,               3*E*Iz/(2*a**2),
                        0,        -3*E*Iz/(2*a**3), 0,
                        0,         0,               3*E*Iz/(2*a**2)],
                       [0,         0,               3*E*Iy/(2*a**3),
                        0,        -3*E*Iy/(2*a**2), 0,
                        0,         0,              -3*E*Iy/(2*a**3),
                        0,        -3*E*Iy/(2*a**2), 0],
                       [0,         0,               0,
                        G*J/(2*a), 0,               0,
                        0,         0,               0,
                        -G*J/(2*a),0,               0],
                       [0,         0,              -3*E*Iy/(2*a**2),
                        0,         2*E*Iy/a,        0,
                        0,         0,               3*E*Iy/(2*a**2),
                        0,         E*Iy/a,          0],
                       [0,         3*E*Iz/(2*a**2), 0,
                        0,         0,               2*E*Iz/a,
                        0,        -3*E*Iz/(2*a**2), 0,
                        0,         0,               E*Iz/a],
                       [-A*E/(2*a),0,               0,
                        0,         0,               0,
                        A*E/(2*a), 0,               0,
                        0,         0,               0],
                       [0,        -3*E*Iz/(2*a**3), 0,
                        0,         0,              -3*E*Iz/(2*a**2),
                        0,         3*E*Iz/(2*a**3), 0,
                        0,         0,              -3*E*Iz/(2*a**2)],
                       [0,         0,              -3*E*Iy/(2*a**3),
                        0,         3*E*Iy/(2*a**2), 0,
                        0,         0,               3*E*Iy/(2*a**3),
                        0,         3*E*Iy/(2*a**2), 0],
                       [0,         0,               0,
                       -G*J/(2*a), 0,               0,
                        0,         0,               0,
                        G*J/(2*a), 0,               0],
                       [0,         0,              -3*E*Iy/(2*a**2),
                        0,         E*Iy/a,          0,
                        0,         0,               3*E*Iy/(2*a**2),
                        0,         2*E*Iy/a,        0],
                       [0,         3*E*Iz/(2*a**2), 0,
                        0,         0,               E*Iz/a,
                        0,        -3*E*Iz/(2*a**2), 0,
                        0,         0,               2*E*Iz/a]])
        k = self.rotate(k,dx,dy,dz)  # transfer to global coordinates
        return k

    def check_frame(self,dx):
        if self.frame=='1D':
            if np.sum(abs(np.diff(dx[:,1]))) > 0:
                err_txt = 'error: rotation in z for 1D element'
                err_txt += ' - select frame=2D'
                raise ValueError(err_txt)
        elif self.frame=='2D':
            if np.sum(abs(np.diff(dx[:,2]))) > 0:
                err_txt = 'error: rotation in y for 2D element'
                err_txt += ' - select frame=3D'
                raise ValueError(err_txt)
        
    def rotate(self,k,dx,dy,dz):
        X = np.array([1,0,0]).T
        Y = np.array([0,1,0]).T
        Z = np.array([0,0,1]).T
        T3 = np.matrix([[np.dot(dx,X),np.dot(dx,Y),np.dot(dx,Z)],
                        [np.dot(dy,X),np.dot(dy,Y),np.dot(dy,Z)],
                        [np.dot(dz,X),np.dot(dz,Y),np.dot(dz,Z)]])
        T = np.zeros((self.ndof,self.ndof))                
        for i in range(int(self.ndof/3)):
            T[3*i:3*i+3,3*i:3*i+3] = T3
        k = T.T*k*T
        return k
        
    def assemble(self):  # assemble global stiffness matrix
        self.K = lil_matrix((self.nK,self.nK))  # stiffness matrix
        self.F = np.zeros(self.nK)
        for i,(dx,dy,dz,dl) in enumerate(zip(self.dX,self.dY,self.dZ,self.dL)):
            ke = self.stiffness(dx,dy,dz,dl)
            j = i*int(self.ndof/2)
            self.K[j:j+self.ndof,j:j+self.ndof] += ke
        self.K = self.K.tocsr()
        self.K = self.K.toarray()
        
    def check_dof(self,dof):
        if dof not in self.BC:
            err_txt = 'dof \''+dof+'\' not present in frame='+self.frame
            err_txt +=' ['+', '.join(self.dof)+']'
            raise ValueError(err_txt)
            
    
    def add_nodal_load(self,dof,nodes,values):
        self.check_dof(dof)
        index = self.dof.index(dof)
        self.L[nodes*self.node_dof+index] = values

        '''
        distributed
        point
        gravity
        '''
        
    def addBC(self,dof,nodes):
        self.check_dof(dof)
        self.BC[dof] = np.append(self.BC[dof],np.array(nodes))
        
    def setBC(self):  # colapse constrained dof
        self.extractBC()
        if len(self.BCindex) > 0:
            self.L = np.delete(self.L,self.BCindex)
            for i in range(2):  # rows and cols
                self.K = np.delete(self.K,self.BCindex,i)
        
    def extractBC(self):  # extract matrix indicies for constraints
        self.BCindex = np.array([])
        for j,key in enumerate(self.BC.keys()):
            ui = self.node_dof*self.BC[key]  # node dof
            if len(ui) > 0:
                if key == 'fix':  # fix all dofs
                    for i in range(self.node_dof):
                        self.BCindex = np.append(self.BCindex,ui+i)
                else:
                    self.BCindex = np.append(self.BCindex,ui+j-1)
        self.BCindex = list(map(int,set(self.BCindex)))
        self.Dindex = np.arange(0,self.nK)
        self.Dindex = np.delete(self.Dindex,self.BCindex)
        
    def solve(self):
        bm.assemble()  # assemble stiffness matrix
        bm.setBC()  # remove constrained equations from stiffness + load 
        D = np.zeros(self.nK)
        D[self.Dindex] = np.linalg.solve(bm.K,bm.L)  # global displacment
        for i,disp in enumerate(self.disp):
            self.D[disp] = D[i::self.node_dof]

    def plot(self):
        pl.plot(self.X[:,0],self.X[:,1],'-o')
        pl.plot(self.X[:,0],self.X[:,1]+self.D['y'],'-o')
        pl.axis('equal')
        sns.despine()
        

                 

        
L = 1
Nel = 12
X = np.zeros((Nel+1,3))
X[:,0] = np.linspace(0,L,Nel+1)
X[:,1] = np.linspace(0,L/3,Nel+1)

v = 80

'''
bm = beam(X,frame='1D')
bm.addBC('fix',[0])  # fix left hand node
bm.add_nodal_load('v',bm.nel-1,v)
bm.solve()
bm.plot()
'''
bm = beam(X,frame='2D')
bm.addBC('fix',[0])  # fix left hand node
bm.add_nodal_load('v',bm.nel,v)
bm.solve()
bm.plot()

bm = beam(X,frame='3D')
bm.addBC('fix',[0])  # fix left hand node
bm.add_nodal_load('v',bm.nel,v)
bm.solve()
bm.plot()



           
'''
# [k][v] = [F]       
# [k][v1,t1,v2,t2]' = [P1,M1,P2,M2]'
     
L = 10
     
Nel = 51
dL = L/Nel
n = 4+2*(Nel-1)
k = np.zeros((n,n))

for i in range(Nel):
    k[2*i:2*i+4,2*i:2*i+4] += elem_k(dL)
    
k = k[2:,2:] # fix end


F = np.zeros((n-2,1))  # loading vector
P = [[0.3,-0.1],[1,-0.95],[0.5,1.5]]
for p in P:
    n = int(np.ceil(p[0]*Nel)*2)-1
    F[n] = p[1]

vs = np.linalg.lstsq(k,F)[0]*dL**3/EI
Fo = np.dot(k,np.reshape(vs,(len(vs),1)))

pl.plot(vs[::2],'o-')
'''