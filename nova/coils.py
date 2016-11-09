import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from itertools import count
import seaborn as sns
import matplotlib
import collections
import amigo.geom as geom
from nova.loops import Profile
from nova.config import Setup
from nova.streamfunction import SF
import nova.cross_coil as cc
colors = sns.color_palette('Paired',12)
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from nova.TF.DEMOxlsx import DEMO
    
class PF(object):
    def __init__(self,eqdsk):
        self.set_coils(eqdsk)
        self.categorize_coils()
        
    def set_coils(self,eqdsk):
        self.coil = collections.OrderedDict()
        if eqdsk['ncoil'] > 0: 
            nC = count(0)
            CSindex = np.argmin(eqdsk['rc'])  # CS radius and width
            self.rCS,self.drCS = eqdsk['rc'][CSindex],eqdsk['drc'][CSindex]
            for i,(r,z,dr,dz,I) in enumerate(zip(eqdsk['rc'],eqdsk['zc'],
                                                 eqdsk['drc'],eqdsk['dzc'],
                                                 eqdsk['Ic'])):
                name = 'Coil{:1.0f}'.format(next(nC))

                self.coil[name] = {'r':r,'z':z,'dr':dr,'dz':dz,'I':I,
                                   'rc':np.sqrt(dr**2+dz**2)/2}
                if eqdsk['ncoil'] > 100 and i>=eqdsk['ncoil']-101:
                    print('exit set_coil loop - coils')
                    break
        
    def categorize_coils(self):
        self.CS_coils,self.PF_coils = [],[]
        catogory = np.zeros(len(self.coil),dtype=[('r','float'),('z','float'),
                                            ('index','int'),('name','object')])
        for i,name in enumerate(self.coil):
            catogory[i]['r'] = self.coil[name]['r']
            catogory[i]['z'] = self.coil[name]['z']
            catogory[i]['index'] = i
            catogory[i]['name'] = name
        CSsort = np.sort(catogory,order=['r','z'])  # sort CS, r then z
        CSsort = CSsort[CSsort['r'] < self.rCS+self.drCS]
        PFsort = np.sort(catogory,order='z')  # sort PF,  z
        PFsort = PFsort[PFsort['r'] > self.rCS+self.drCS]
        self.index = {'PF':{'n':len(PFsort['index']),
                            'index':PFsort['index'],'name':PFsort['name']},
                      'CS':{'n':len(CSsort['index']),
                            'index':CSsort['index'],'name':CSsort['name']}}
      
    def unpack_coils(self):
        nc = len(self.coil.keys())
        Ic = np.zeros(nc)
        rc,zc,drc,dzc = np.zeros(nc),np.zeros(nc),np.zeros(nc),np.zeros(nc)
        names = []
        for i,name in enumerate(self.coil.keys()):
            rc[i] = self.coil[name]['r']
            zc[i] = self.coil[name]['z']
            drc[i] = self.coil[name]['dr']
            dzc[i] = self.coil[name]['dz']
            Ic[i] = self.coil[name]['I']
            names.append(name)
        return nc,rc,zc,drc,dzc,Ic,names
        
    def plot_coil(self,coils,label=False,current=False,coil_color=None,fs=12,
                  alpha=1):
        if coil_color==None:
            color = colors
        else:
            color = coil_color  # color itterator
        
        for i,name in enumerate(coils.keys()):
            coil = coils[name]
            r,z,dr,dz = coil['r'],coil['z'],coil['dr'],coil['dz']
            Rfill = [r+dr/2,r+dr/2,r-dr/2,r-dr/2]
            Zfill = [z-dz/2,z+dz/2,z+dz/2,z-dz/2]
            if coil['I'] != 0:
                edgecolor = 'k'
            else:
                edgecolor = 'r'  
            coil_color = color[4] 
            if name in self.index['CS']['name']:
                drs = -2.5/3*dr
                ha = 'right'
                coil_color = color[5]
            elif name in self.index['PF']['name']:
                drs = 2.5/3*dr
                ha = 'left'
                coil_color = color[5]
            pl.fill(Rfill,Zfill,facecolor=coil_color,alpha=alpha,
                    edgecolor=edgecolor)
            
            if label and current:
                zshift = max([coil['dz']/4,0.4])
            else:
                zshift = 0
            if label: 
                pl.text(r+drs,z+zshift,name,fontsize=fs*1.1,
                        ha=ha,va='center',color=0.5*np.ones(3))
            if current: 
                pl.text(r+drs,z-zshift,'{:1.1f}MA'.format(coil['I']*1e-6),
                        fontsize=fs*1.1,ha=ha,va='center',
                        color=0.5*np.ones(3))
                                  
    def plot(self,color=None,coils=None,label=False,plasma=False,
                   current=False,alpha=1):
        fs = matplotlib.rcParams['legend.fontsize']
        if coils is None:
            coils = self.coil
        self.plot_coil(coils,label=label,current=current,fs=fs,
                       coil_color=color,alpha=alpha)
        if plasma:
            coils = self.plasma_coil                
            self.plot_coil(coils,coil_color=color,alpha=alpha)
            
    def coil_corners(self,coils):
        R,Z = np.array([]),np.array([])
        Nc = len(coils['id'])
        dR,dZ = np.zeros(Nc),np.zeros(Nc)
        if len(coils['dR'])>0: 
            dR[Nc-len(coils['dR']):]=coils['dR']
        if len(coils['dZ'])>0: 
            dZ[Nc-len(coils['dZ']):]=coils['dZ']
        for Cid,Cdr,Cdz in zip(coils['id'],dR,dZ):
            r = self.coil['Coil'+str(Cid)]['r']
            z = self.coil['Coil'+str(Cid)]['z']
            dr = self.coil['Coil'+str(Cid)]['dr']
            dz = self.coil['Coil'+str(Cid)]['dz']
            if Cdr==0 and Cdz==0:
                R = np.append(R,[r+dr/2,r+dr/2,r-dr/2,r-dr/2])
                Z = np.append(Z,[z+dz/2,z-dz/2,z+dz/2,z-dz/2])
            else:
                R = np.append(R,r+Cdr)
                Z = np.append(Z,z+Cdz)
        return R,Z

    def fit_coils(self,Cmove,dLo=0.1):
        coils = collections.OrderedDict()
        for side in Cmove.keys():
            if isinstance(dLo,(list,tuple)):
                if not 'in' in side:
                    dL = dLo[0]
                else:
                    dL = dLo[-1]
            else:
                dL = dLo
            for index in Cmove[side]:
                coil = 'Coil'+str(index)
                dr,dz = self.Cshift(self.sf.coil[coil],side,dL)
                coils[coil] = self.sf.coil[coil]
                coils[coil]['r'] = coils[coil]['r']+dr
                coils[coil]['z'] += dz
                coils[coil]['shiftr'] = dr
                coils[coil]['shiftz'] = dz
                
        with open('../Data/'+self.conf.config+'_coil_shift.txt','w') as f:
            f.write('Ncoils = {:1.0f}\n\n'.format(len(coils.keys())))
            f.write('Name\tR[m]\t\tZ[m]\t\tshiftR[m]\tshiftZ[m]\n')
            index = sorted(map((lambda s: int(s.strip('Coil'))),coils.keys()))
            for i in index:
                coil = 'Coil'+str(i)
                f.write(coil+'\t{:1.6f}\t{:1.6f}\t{:1.6f}\t{:1.6f}\n'.format(\
                coils[coil]['r'],coils[coil]['z'],\
                coils[coil]['shiftr'],coils[coil]['shiftz']))                                   
        return coils

class TF(object):
    
    def __init__(self,profile,**kwargs):
        self.initalise_loops()  # initalise loop family
        self.profile = profile  # loop profile
        self.cross_section(**kwargs)  # coil cross-sections
        self.get_loops(profile.loop.draw())

    def cross_section(self,J=18.25,twall=0.045,**kwargs):  # MA/m2 TF current density
        self.section = {}
        self.section['case'] = {'side':0.1,'nose':0.51,'inboard':0.04,
                                'outboard':0.19,'external':0.225}
        if 'sf' in kwargs and hasattr(self.profile,'nTF'):
            sf = kwargs.get('sf')
            BR = sf.eqdsk['bcentr']*sf.eqdsk['rcentr']
            Iturn = 1e-6*abs(2*np.pi*BR/(self.profile.nTF*cc.mu_o))
            Acs = Iturn/J
            ro = np.min(self.profile.loop.draw()['r'])
            rwp1 = ro-self.section['case']['inboard']
            theta = np.pi/self.profile.nTF
            rwall = twall/np.sin(theta)
            depth = np.tan(theta)*(rwp1-rwall+\
            np.sqrt((rwall-rwp1)**2-4*Acs/(2*np.tan(theta))))
            width = Acs/depth
            self.section['winding_pack'] = {'width':width,'depth':depth}
        else:
            self.section['winding_pack'] = {'width':0.625,'depth':1.243}
        self.rc = self.section['winding_pack']['width']/2

    def initalise_loops(self):
        self.x = {}
        for loop in ['in','wp_in','cl','wp_out','out','nose','loop']:
            self.x[loop] = {'r':[],'z':[]}

    def transition_index(self,r_in,z_in,eps=1e-6):
        npoints = len(r_in)
        r_cl = r_in[0]+eps
        upper = npoints-next((i for i,r_in_ in enumerate(r_in) if r_in_>r_cl))+1
        lower = next((i for i,r_in_ in enumerate(r_in) if r_in_>r_cl))
        top,bottom = np.argmax(z_in), np.argmin(z_in)
        index = {'upper':upper,'lower':lower,'top':top,'bottom':bottom}
        return index
        
    def loop_dt(self,r,z,dt_in,dt_out,index):
        l = geom.length(r,z)
        L = np.array([0,l[index['lower']],l[index['bottom']],
                      l[index['top']],l[index['upper']],1])
        dR = np.array([dt_in,dt_in,dt_out,dt_out,dt_in,dt_in])
        dt = interp1d(L,dR)(l)
        return dt
 
    def get_loops(self,x):
        r,z = x['r'],x['z']
        wp = self.section['winding_pack']
        case = self.section['case']
        inboard_dt = [case['inboard'],wp['width']/2,wp['width']/2,case['nose']]
        outboard_dt = [case['outboard'],wp['width']/2,wp['width']/2,
                       case['external']]
        loops = ['wp_in','cl','wp_out','out']
        self.x['in']['r'],self.x['in']['z'] = r,z
        index = self.transition_index(self.x['in']['r'],self.x['in']['z'])
        for loop,dt_in,dt_out in zip(loops,inboard_dt,outboard_dt):
            dt = self.loop_dt(r,z,dt_in,dt_out,index)
            r,z = geom.offset(r,z,dt,close_loop=True)
            self.x[loop]['r'],self.x[loop]['z'] = r,z
        return self.x
        
    def split_loop(self):  # split inboard/outboard for fe model
        r,z = self.x['cl']['r'],self.x['cl']['z']
        index = self.transition_index(r,z)
        upper,lower = index['upper'],index['lower']
        self.x['nose']['r'] = np.append(r[upper:],r[1:lower+1])
        self.x['nose']['z'] = np.append(z[upper:],z[1:lower+1])
        self.x['loop']['r'] = r[lower+1:upper]
        self.x['loop']['z'] = z[lower+1:upper]

    def loop_interpolators(self,trim=[0,1],offset=0.75):  # outer loop coordinate interpolators
        r,z = self.x['cl']['r'],self.x['cl']['z']
        index = self.transition_index(r,z)
        self.fun = {'in':{},'out':{}}
        for side,sign in zip(['in','out','cl'],[-1,1,1]):  # inner/outer loop offset
            r,z = self.x[side]['r'],self.x[side]['z']
            r = r[index['lower']+1:index['upper']]
            z = z[index['lower']+1:index['upper']]
            r,z = geom.offset(r,z,sign*offset)
            l = geom.length(r,z)
            l[(l>0.5) & (l<0.7)] = 0.5
            lt = np.linspace(trim[0],trim[1],int(np.diff(trim)*len(l)))
            r,z = interp1d(l,r)(lt),interp1d(l,z)(lt)
            l = np.linspace(0,1,len(r))
            self.fun[side] = {'r':IUS(l,r),'z':IUS(l,z)}
            self.fun[side]['L'] = geom.length(r,z,norm=False)[-1]
            self.fun[side]['dr'] = self.fun[side]['r'].derivative(n=1)
            self.fun[side]['dz'] = self.fun[side]['r'].derivative(n=1)
     
    def Cshift(self,coil,side,dL):  # shift pf coils to tf track
        if 'in' in side:
            R,Z = self.x['in']['r'],self.x['in']['z']
        else:
            R,Z = self.x['out']['r'],self.x['out']['z']
        rc,zc = coil['r'],coil['z']
        drc,dzc = coil['dr'],coil['dz']
        rp = rc+np.array([-drc,drc,drc,-drc])/2
        zp = zc+np.array([-dzc,-dzc,dzc,dzc])/2
        nR,nZ = geom.normal(R,Z)[:2]       
        mag = np.sqrt(nR**2+nZ**2)
        nR /= mag
        nZ /= mag
        i = []
        L = np.empty(len(rp))
        dn = np.empty((2,len(rp)))
        for j,(r,z) in enumerate(zip(rp,zp)):
            i.append(np.argmin((R-r)**2+(Z-z)**2))
            dr = [r-R[i[-1]],z-Z[i[-1]]]  
            dn[:,j] = [nR[i[-1]],nZ[i[-1]]]
            L[j] = np.dot(dr,dn[:,j])
            if 'in' in side: L[j] *= -1
        jc = np.argmin(L)
        fact = dL-L[jc]
        if 'in' in side: fact *= -1
        delta = fact*dn[:,jc]
        return delta[0],delta[1]
        
    def get_loop(self,expand=0):  # generate boundary dict for elliptic
        R,Z = self.x['cl']['r'],self.x['cl']['z']
        boundary = {'R':R,'Z':Z,'expand':expand}
        return boundary

    def fill(self,write=False,plot=True,alpha=1):
        geom.polyparrot(self.x['in'],self.x['wp_in'],
                        color=0.4*np.ones(3),alpha=alpha)
        geom.polyparrot(self.x['wp_in'],self.x['wp_out'],
                        color=0.6*np.ones(3),alpha=alpha)
        geom.polyparrot(self.x['wp_out'],self.x['out'],
                        color=0.4*np.ones(3),alpha=alpha)
        pl.plot(self.x['cl']['r'],self.x['cl']['z'],'-.',color=0.5*np.ones(3))
        pl.axis('equal')
        pl.axis('off')

    def support(self,**kwargs):
        self.rzGet()
        self.fill(**kwargs)
    '''    
    def volume_ratio(self,Rp,Zp):
        self.Pvol = loop_vol(Rp,Zp,plot=False)
        Rtf,Ztf,L = self.drawTF(self.xCoil, npoints=300)
        self.TFvol = loop_vol(Rtf,Ztf,plot=False)
        self.Rvol = self.TFvol/self.Pvol
       
    def length_ratio(self):
        self.Plength = self.length(self.Rp,self.Zp,norm=False)[-1]
        Rtf,Ztf,L = self.drawTF(self.xCoil, npoints=300)
        R,Z = geom.offset(Rtf,Ztf,self.dRsteel+self.dRcoil/2)
        self.TFlength = self.length(R,Z,norm=False)[-1]
        self.Rlength = self.TFlength/self.Plength
    '''    
if __name__ is '__main__':  # test functions

    nPF,nTF = 6,16
    config = {'TF':'SN','eq':'SN_{:d}PF_{:d}TF'.format(nPF,nTF)}
    setup = Setup(config['eq'])
    sf = SF(setup.filename)
    profile = Profile(config['TF'],family='S',part='TF',nTF=nTF,obj='L',
                      load=True)
    #profile.loop.plot()
    tf = TF(profile,sf=sf)
    tf.fill()

    demo = DEMO()
    demo.fill_part('Vessel')
    demo.fill_part('Blanket')
    demo.plot_ports()
  

    rp,zp = demo.port['P0']['right']['r'],demo.port['P0']['right']['z']                  
    pl.plot(rp,zp,'k',lw=3)
    '''
    xc = [rp[0],zp[0]]
    nhat = np.array([rp[1]-rp[0],zp[1]-zp[0]])
    sal.tf.loop_interpolators(offset=0)  # construct TF interpolators
    loop = sal.tf.fun['out']

    Lo = minimize_scalar(SALOME.OIS_placment,method='bounded',
                         args=(loop,(rp[0],zp[0])),bounds=[0,1]).x
    xo = [Lo,0]             
    L = minimize(SALOME.intersect,xo,method='L-BFGS-B', 
                 bounds=([0.1,0.9],[0,15]),args=(xc,nhat,loop)).x
    pl.plot(loop['r'](L[0]),loop['z'](L[0]),'o')
    '''


    
