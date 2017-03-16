import pylab as pl
from nova.config import Setup,select
from nova.streamfunction import SF
from nova.radial_build import RB
from nova.elliptic import EQ
from nova.coils import PF,TF
from nova.inverse import INV,SWING
import numpy as np
from time import time
import amigo.geom as geom
from nova import loops
from nova.loops import Profile,plot_oppvar
from nova.shape import Shape
from nova.DEMOxlsx import DEMO
from nova.force import force_feild
from nova.firstwall import divertor,main_chamber
from amigo.IO import trim_dir
from nova.coil_cage import coil_cage

import seaborn as sns
rc = {'figure.figsize':[5,5*16/12],'savefig.dpi':150, # 
      'savefig.jpeg_quality':200,'savefig.pad_inches':0.1,
      'lines.linewidth':1.5}
sns.set(context='talk',style='white',font='sans-serif',palette='Set2',
        font_scale=5/8,rc=rc)


nTF,ripple = 18,True
base = {'TF':'demo_nTF','eq':'DEMO_SN_EOF'}    
config,setup = select(base,nTF=nTF,update=False)     
sf = SF(setup.filename) 
pf = PF(sf.eqdsk) 
  
   
pf.plot(coils=pf.coil,color=0.75*np.ones(3),label=False,plasma=False)


demo = DEMO()
profile = Profile(config['TF'],family='S',part='TF',nTF=nTF,obj='L')
shp = Shape(profile,eqconf=config['eq'],ny=8)
shp.add_vessel(demo.parts['Vessel']['out'])
#shp.minimise(ripple=ripple,verbose=False)

tf = TF(x_in=demo.parts['TF_Coil']['in'],nTF=nTF,sf=sf)  
tf.x['out'] = demo.parts['TF_Coil']['out']
tf.x['cl']['r'],tf.x['cl']['z'] = geom.rzSLine(tf.x['cl']['r'],tf.x['cl']['z'])
config['eqdsk'] += '_baseline'
#tf = shp.tf

cage = coil_cage(nTF=nTF,rc=tf.rc,plasma={'config':config['eq']},ny=shp.ny,
                 smooth=False)
cage.set_TFcoil(tf.x['cl'])
cage.output()  

cage.pattern(plot=True)

s = np.zeros(3)
s[0],s[2] = sf.LFPr,sf.LFPz
ripple = cage.point(s,variable='ripple')
print('ripple',ripple)
'''

inv = INV(pf,tf,dCoil=0.5)
sc = scenario(inv,sf)
sc.flat_top(n=31)

sc.output()

swing = 'SOF'
if swing == 'SOF':
    swing_index = 0
else:
    swing_index = -1
    config['eqdsk'] = config['eqdsk'].replace('SOF','EOF')
inv.solve_slsqp(inv.swing['flux'][swing_index])
inv.eq.run(update=False)

          

#demo.fill_part('Blanket')
demo.fill_part('Vessel')
#sf.contour()
pf.plot(coils=pf.coil,label=False,current=True,plasma=False,alpha=0.75)
#pf.plot(coils=pf.sub_coil)
inv.ff.plot(scale=3)

cage.plot_contours(variable='ripple',n=2e3,loop=demo.fw)
tf.fill()

sf.eqwrite(pf,config=config['eqdsk'])
pl.savefig('../../Figs/{}.pdf'.format(config['eqdsk']))

'''
