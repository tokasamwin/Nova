import numpy as np
import pylab as pl
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import scipy.optimize as op
from itertools import cycle,count
import seaborn as sns
import matplotlib
import collections
from amigo.geom import Loop
import amigo.geom as geom
import nova.geqdsk
from nova.inductance import neumann
from nova.loops import Acoil,Dcoil,Scoil
from nova.TF.ripple import ripple
from nova.config import Setup,trim_dir
colors = sns.color_palette('Set2',12)
import pandas as pd
    
class PF(object):
    def __init__(self,eqdsk):
        self.set_coils(eqdsk)
        
    def set_coils(self,eqdsk):
        self.coil = {}
        if eqdsk['ncoil'] > 0: 
            nC = count(0)
            for i,(r,z,dr,dz,I) in enumerate(zip(eqdsk['rc'],eqdsk['zc'],
                                                 eqdsk['drc'],eqdsk['dzc'],
                                                 eqdsk['Ic'])):
                name = 'Coil{:1.0f}'.format(next(nC))
                self.coil[name] = {'r':r,'z':z,'dr':dr,'dz':dz,'I':I,
                                   'rc':np.sqrt(dr**2+dz**2)/2}
                if i>=10:
                    print('exit set_coil loop - coils')
                    break
         
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
        
    def plot_coil(self,coils,label=False,current=False,coil_color=None,fs=12):
        if coil_color==None:
            color = cycle(sns.color_palette('Set2',len(coils.keys())))
        else:
            color = coil_color  # color itterator
        
        color = sns.color_palette('Set3',300)
        r_med = sp.median([coils[coil]['r'] for coil in coils])
        for i,name in enumerate(coils.keys()):
            coil = coils[name]
            if label and current:
                zshift = coil['dz']/4
            else:
                zshift = 0
            r,z,dr,dz = coil['r'],coil['z'],coil['dr'],coil['dz']
            Rfill = [r+dr/2,r+dr/2,r-dr/2,r-dr/2]
            Zfill = [z-dz/2,z+dz/2,z+dz/2,z-dz/2]
            if coil['I'] != 0:
                edgecolor = 'k'
            else:
                edgecolor = 'r'  
            #coil_color = next(color)
            try:
                ic = int(name.split('_')[0].split('Coil')[-1])
            except:
                ic = 0  # plasma coils
            coil_color = color[ic] 
            coil_color = color[8]
            pl.fill(Rfill,Zfill,facecolor=coil_color,alpha=1,
                    edgecolor=edgecolor)
            if label: 
                pl.text(r,z+zshift,name.replace('Coil',''),fontsize=fs*0.6,
                        ha='center',va='center',color=[0.25,0.25,0.25])
            if current: 
                if r < r_med:
                    drs = -2/3*dr
                    ha = 'right'
                else:
                    drs = 2/3*dr
                    ha = 'left'
                pl.text(r+drs,z-zshift,'{:1.1f}MA'.format(coil['I']*1e-6),
                        fontsize=fs*0.6,ha=ha,va='center',
                        color=[0.25,0.25,0.25])
                                  
    def plot(self,color=None,coils=None,label=False,plasma=False,
                   current=False):
        fs = matplotlib.rcParams['legend.fontsize']
        if coils is None:
            coils = self.coil
        self.plot_coil(coils,label=label,current=current,fs=fs,
                       coil_color=color)
        if plasma:
            coils = self.plasma_coil                
            self.plot_coil(coils,coil_color=color)
            
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
              
    def Cshift(self,coil,side,dL):
        if 'in' in side:
            R,Z = self.TFinR,self.TFinZ
        else:
            R,Z = self.TFoutR,self.TFoutZ
        rc,zc = coil['r'],coil['z']
        drc,dzc = coil['dr'],coil['dz']
        rp = rc+np.array([-drc,drc,drc,-drc])/2
        zp = zc+np.array([-dzc,-dzc,dzc,dzc])/2
        nR,nZ = self.normal(R,Z)       
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
    
    def __init__(self,config,coil_type='S',npoints=100):
        self.x = {}
        for loop in ['in','wp_in','cl','wp_out','out','nose','loop']:
            self.x[loop] = {'r':[],'z':[]}
        self.npoints = npoints  # number of loop divisions
        self.cross_section()
        self.config = config  # save file prefix
        self.load_coil(coil_type)  # initalize coil object
        data_dir = trim_dir('../../Data/')
        self.dataname = data_dir+self.config+'_TF.pkl'
        self.read_coil_dict()
        
    def read_coil_dict(self):
        try:
            with open(self.dataname,'rb') as input:
                self.coil_dict = pickle.load(input)
        except:
            print('file '+self.dataname+' not found')
            print('initializing new coil_dict')
            self.coil_dict = {}
        self.frame_data()
            
    def frame_data(self):
        self.data_frame = {}   
        for coil_type in self.coil_dict:
            data = {}
            for nTF in self.coil_dict[coil_type]:
                data[nTF] = {}
                for obj in self.coil_dict[coil_type][nTF]:
                    data[nTF][obj] = True
            self.data_frame[coil_type] = pd.DataFrame(data)
 
    def load(self,nTF='unset',objective='unset'):
        if objective in self.coil_dict.get(self.coil_type,{}).get(nTF,{}):
            coil_dict = self.coil_dict[self.coil_type][nTF][objective]
            for key in coil_dict:
                if hasattr(self.coil,key):
                    setattr(self.coil,key,coil_dict[key])
        else:
            errtxt = '\n'
            errtxt += 'data not found:\n'
            errtxt += 'coil type {}, nTF {}, objective {}\n'.\
            format(self.coil_type,nTF,objective)
            errtxt += self.avalible_data(verbose=False)
            raise ValueError(errtxt)
        self.get_loops(self.coil.draw())
        
    def avalible_data(self,verbose=True):
        if len(self.coil_dict) == 0:
            datatxt = 'no data avalible'
        else:
            datatxt ='\n{}: data avalible [objective,nTF]'.format(self.config)
            for coil_type in self.data_frame:
                datatxt += '\n\ncoil type {}:\n{}'.\
                format(coil_type,self.data_frame[coil_type].fillna(''))
        if verbose:
            print(datatxt)
        else:
            return datatxt
            
    def clear_data(self):
        with open(self.dataname, 'wb') as output:
            self.coil_dict = {}
            pickle.dump(self.coil_dict,output,-1)
        
    def write(self,nTF='unset',objective='unset'):  # write xo and oppvar to fil
        if self.coil_type in self.coil_dict:
            if nTF not in self.coil_dict[self.coil_type]:
                self.coil_dict[self.coil_type][nTF] = {objective:[]}
        else:
            self.coil_dict[self.coil_type] = {nTF:{objective:[]}}
        cdict = {}
        for key in ['xo','oppvar','coil_type','symetric','tension','limits']:
            if hasattr(self.coil,key):
                cdict[key] = getattr(self.coil,key)
        self.coil_dict[self.coil_type][nTF][objective] = cdict
        with open(self.dataname, 'wb') as output:
            pickle.dump(self.coil_dict,output,-1)
        self.frame_data()

    def cross_section(self):
        self.section = {}
        self.section['winding_pack'] = {'width':0.625,'depth':1.243}
        self.section['case'] = {'side':0.1,'nose':0.51,'inboard':0.04,
                                'outboard':0.19,'external':0.225}
        self.rc = self.section['winding_pack']['width']/2
        
    def load_coil(self,coil_type):
        self.coil_type = coil_type  # A==arc, D==Princton-D, S==spline
        if self.coil_type == 'A':  # tripple arc (5-7 parameter)
            self.coil = Acoil(npoints=self.npoints)
        elif self.coil_type == 'D':  # princton D (3 parameter)
            self.coil = Dcoil(npoints=self.npoints)
        elif self.coil_type == 'S':  # polybezier (8-16 parameter)
            self.coil = Scoil(npoints=self.npoints)
        else:
            errtxt = '\n'
            errtxt += 'coil type \''+self.coil_type+'\'\n'
            errtxt += 'select from [A,D,S]\n'
            raise ValueError(errtxt)

    def setargs(self,shape):  # set input data based on shape contents

        if 'R' in shape and 'Z' in shape:  # define coil by R,Z arrays
            self.datatype = 'array'
            self.Rcl,self.Zcl = shape['R'],shape['Z']  # coil centreline
            
        elif 'vessel' in shape and 'pf' in shape:
            self.datatype = 'fit'  # fit coil to internal incoil + pf coils
            self.set_bound(shape)  # TF boundary conditions
            #self.bound['ro_min'] -= 0.5  # offset minimum radius
            if 'plot' in shape:
                if shape['plot']:
                    self.plot_bounds() 
        elif 'config' in shape and 'coil_type' in shape:
            self.datatype = 'file'  # load coil from file
            self.setup = Setup(shape.get('config'))
            self.setup.coils['external']['id']=[]
            self.dataname = self.setup.dataname  # TF coil geometory
            self.filename = self.setup.filename  # eqdsk file
        else:
            errtxt = '\n'
            errtxt += 'unable to load TF coil\n'
            raise ValueError(errtxt)
            
    def transition_index(self,r_in,z_in,eps=1e-6):
        npoints = len(r_in)
        r_cl = r_in[0]+eps
        upper = npoints-next((i for i,r_in_ in enumerate(r_in) if r_in_>r_cl))-1
        lower = next((i for i,r_in_ in enumerate(r_in) if r_in_>r_cl))-1
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
        upper,lower = index['upper']+2,index['lower']+1
        self.x['nose']['r'] = np.append(r[upper:],r[1:lower+1])
        self.x['nose']['z'] = np.append(z[upper:],z[1:lower+1])
        self.x['loop']['r'] = r[lower+1:upper]
        self.x['loop']['z'] = z[lower+1:upper]
        
    def get_loop(self,expand=0):  # generate boundary dict for elliptic
        R,Z = self.x['cl']['r'],self.x['cl']['z']
        boundary = {'R':R,'Z':Z,'expand':expand}
        return boundary

    def fill(self,write=False,plot=True):
        geom.polyloopfill(self.x['in'],self.x['wp_in'],color=0.4*np.ones(3))
        geom.polyloopfill(self.x['wp_in'],self.x['wp_out'],color=0.6*np.ones(3))
        geom.polyloopfill(self.x['wp_out'],self.x['out'],color=0.4*np.ones(3))
        pl.plot(self.x['cl']['r'],self.x['cl']['z'],'-.',color=0.5*np.ones(3))
        pl.axis('equal')
        pl.axis('off')

        '''
        if write:
            with open('../Data/'+self.dataname+'_TFcoil.txt','w') as f:
                f.write('Rin m\t\tZin m\t\tRout m\t\tZout m\n')
                for rin,zin,rout,zout in zip(self.TFinR,self.TFinZ,
                                             self.TFoutR,self.TFoutZ):
                    f.write('{:1.6f}\t{:1.6f}\t{:1.6f}\t{:1.6f}\n'.format(\
                    rin,zin,rout,zout))            
        '''
    def support(self,**kwargs):
        self.rzGet()
        self.fill(**kwargs)
        
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
        
if __name__ is '__main__':  # test functions

    config = 'DEMO_SN'
    tf = TF(config,coil_type='S',npoints=80)
    
    tf.avalible_data()
    
    
    #tf.coil.set_input(inputs={'l':1.5})
    #tf.write(nTF=18,objective='E')
    
    tf.load(nTF=18,objective='L')
    
    tf.split_loop()
    #pl.plot(tf.x['cl']['r'],tf.x['cl']['z'])
    pl.plot(tf.x['nose']['r'],tf.x['nose']['z'])
    pl.plot(tf.x['loop']['r'],tf.x['loop']['z'])
    
    tf.load(nTF=16,objective='L')
    #pl.plot(tf.x['cl']['r'],tf.x['cl']['z'])
    
    
    pl.axis('equal')
    pl.axis('off')