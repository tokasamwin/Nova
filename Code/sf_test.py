import pylab as pl
from streamfunction import SF
from radial_build import RB
from elliptic import EQ

import seaborn as sns
rc = {'figure.figsize':[10,10*12/16],'savefig.dpi':100, # 
      'savefig.jpeg_quality':100,'savefig.pad_inches':0.1,
      'lines.linewidth':2}
sns.set(context='talk',style='white',font='sans-serif',palette='Set2',
        font_scale=7/8,rc=rc)

pl.clf()


filename = 'Equil_AR3d1_16coils_SFminus_v4_2015_09_bt_1d03li_0d8_Ipl_20d25_SOF.eqdsk'
#filename = 'Equil_AR3d1_16coils_v3_bt_1d03li_0d8_Ipl_20d25_XDext_v4_NOEDDY_SOF_FINAL3_v1.eqdsk'
sf = SF(filename,config='X',dSOL=3e-3,Nsol=51)
sf.get_sol_psi()

rb = RB(sf,Np=250)





eq = EQ(sf,sigma=0,limit=[6,10,-7,-3],n=1e4)  # resample



eq.plotj()

levels = sf.contour()


sf.sol(plot=True)
