import pylab as pl
from config import Setup
from streamfunction import SF
from radial_build import RB
from coils import PF,TF
import geqdsk



setup = Setup('SXex')
sf = SF(setup.filename)
rb = RB(setup,sf)


pf = PF(sf.eqdsk)

tf = TF()

xCoil = [ 5.63566835,  2.5318409 ,  4.4570691 ,  4.52603387, -2.61044312]
r,z,l = tf.drawTF(xCoil)
pl.plot(r,z)
pl.axis('equal')



pf.plot_coils(coils=pf.coil,label=True,plasma=False,current=False) 

#tf.TFcoil(False)