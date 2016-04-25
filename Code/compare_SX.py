import pylab as pl
import numpy as np
import matplotlib as mpl
from itertools import cycle
from eqConfig import Config
from radial_build import RB
from streamfunction import SF
from feild_calc import solCalc
#from elliptic import grid
#import cross_coil as cc
import seaborn as sns
rc = {'figure.figsize':[2*3.14*12/16,2*3.14],'savefig.dpi':150, # 
      'savefig.jpeg_quality':100,'savefig.pad_inches':0.1,
      'lines.linewidth':1.3}
sns.set(context='talk',style='white',font='sans-serif',palette='Set2',
        font_scale=7/8,rc=rc)
Color = cycle(sns.color_palette('Set2'))

figheight = 8.0
mpl.rcParams['figure.figsize'] = [figheight,figheight*5/9]
xlim,ylim = [6,12],[-7,-6]

pl.figure()
pl.axis('equal')
pl.axis('off')
pl.xlim(xlim)
pl.ylim(ylim)

conf = Config('SX7')
sf = SF(conf,sample=1)
conf.TF(sf)
rb = RB(conf,sf,Np=1e3)
rb.divertor_outline(False,color=next(Color))
Rt,Zt = rb.targets['outer']['R'],rb.targets['outer']['Z']
Ro,Zo = rb.targets['outer']['Rsol'][-1],rb.targets['outer']['Zsol'][-1]
index = np.argmin((Rt-Ro)**2+(Zt-Zo)**2)
R,Z = Rt[:index+1],Zt[:index+1]
graze = rb.sol.get_graze(R,Z)
print(graze*180/np.pi,rb.targets['outer']['theta']*180/np.pi)

conf = Config('SX8')
sf = SF(conf,sample=1)
conf.TF(sf)
rb = RB(conf,sf,Np=1e3)
rb.divertor_outline(False,color=next(Color))

Ro,Zo = rb.targets['outer']['Rsol'][-1],rb.targets['outer']['Zsol'][-1]
index = np.argmin((Rt-Ro)**2+(Zt-Zo)**2)
R,Z = Rt[:index+1],Zt[:index+1]
graze = rb.sol.get_graze(R,Z)
print(graze*180/np.pi,rb.targets['outer']['theta']*180/np.pi)
