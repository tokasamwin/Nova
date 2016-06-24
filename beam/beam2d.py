import numpy as np
import pylab as pl

EI = 1


# element stiffness matrix
def elem_k(dL):
    k = np.matrix([[12,6*dL,-12,6*dL],
                   [6*dL,4*dL**2,-6*dL,2*dL**2],
                   [-12,-6*dL,12,-6*dL],
                   [6*dL,2*dL**2,-6*dL,4*dL**2]])
    return k
               
               
               
# [k][v] = [F]       
# [k][v1,t1,v2,t2]' = [P1,M1,P2,M2]'
     
L = 10
     
Nel = 201
dL = L/Nel
n = 4+2*(Nel-1)
k = np.zeros((n,n))

for i in range(Nel):
    k[2*i:2*i+4,2*i:2*i+4] += elem_k(dL)
    
k = k[2:,2:] # fix end

#print(np.linalg.matrix_rank(k))
#print(np.shape(k))

F = np.zeros((n-2,1))

P = [[0.3,-0.1],[1,-0.95],[0.5,1.5]]
for p in P:
    n = int(np.ceil(p[0]*Nel)*2)-1
    F[n] = p[1]
    print(n)

#F[-1] = 1.5

#F[int(n/2)] = 5
#F[-2] = -0.5

vs = np.linalg.lstsq(k,F)[0]*dL**3
Fo = np.dot(k,np.reshape(vs,(len(vs),1)))

#print(Fo) 

#print(np.sum(vs[::2]))

pl.plot(vs[::2],'o-')

#pl.plot(vs[1::2],'o-')
#pl.plot(0.2*np.cumsum(vs[1::2]),'o-')

'''
ks = np.append(ke,np.zeros((2,np.shape(ke)[1])),axis=0)          
ks[4,0] = 1
ks[5,1] = 1

F = np.array([0,0,1,0,0,0])
vs = np.linalg.lstsq(ks,F)[0]
#print(vs)           

ko = ke[2:,2:]
Fo = np.array([1,0])
vo = np.linalg.solve(ko,Fo)
#print(vo)
               

dL = 1
ke = np.matrix([[12,6*dL,-12,6*dL],
               [6*dL,4*dL**2,-6*dL,2*dL**2],
               [-12,-6*dL,12,-6*dL],
               [6*dL,2*dL**2,-6*dL,4*dL**2]])           

k = np.zeros((6,6))  # two elements
k[:4,:4] += ke
k[2:,2:] += ke

v = np.linalg.solve(k[2:,2:],[10,0,0,0])
print(v)


k = np.append(k,np.zeros((1,6)),axis=0)
k = np.append(k,np.zeros((1,6)),axis=0)
#k = np.append(k,np.zeros((1,6)),axis=0)


k[6,0] = 1
k[7,1] = 1

F = np.zeros(8)
F[4] = 1
v = np.linalg.lstsq(k,F)

print(v)
'''

