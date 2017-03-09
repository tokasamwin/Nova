import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline as spline
import seaborn as sns
colors = sns.color_palette('Paired',12)
from itertools import cycle
import amigo.geom as geom
from nova.loops import Profile
from nova.shape import Shape
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.optimize import minimize
from collections import OrderedDict
from amigo.IO import trim_dir
import json
from nova.firstwall import divertor,main_chamber

class RB(object):
    def __init__(self,sf,setup,npoints=500):
        self.setup = setup
        self.sf = sf
        self.npoints = npoints
        self.dataname = setup.configuration
        self.datadir = trim_dir('../../../Data/') 
        self.segment = {}  # store section segments (divertor,fw,blanket...)
        
    def generate(self,mc,mode='calc',color=0.6*np.ones(3),
                  plot=True,debug=False,symetric=False,DN=False):
        
        if mode == 'eqdsk':  # read first wall from eqdsk
            self.Rb,self.Zb = self.sf.xlim,self.sf.ylim
        elif mode == 'calc':
            div = divertor(self.sf,self.setup,debug=debug)
            if DN:  # double null
                self.sf.get_Xpsi(select='upper')  # upper X-point
                div.place(debug=debug)
                self.segment = div.join(mc)
                ru = self.segment['first_wall']['r']
                zu = self.segment['first_wall']['z']
                blanket_u = self.place_blanket(select='upper',store=False)[0]
                self.sf.get_Xpsi(select='lower')  # lower X-point
                div.place(debug=debug)
                self.segment = div.join(mc)
                rl = self.segment['first_wall']['r']
                zl = self.segment['first_wall']['z']
                blanket_l = self.place_blanket(select='lower',store=False)[0]
                r,z = self.upper_lower(ru,zu,rl,zl,Zjoin=self.sf.Mpoint[1])
                self.segment['first_wall'] = {'r':r,'z':z} 
                rfw,zfw = self.upper_lower(blanket_u[0]['r'][::-1],
                                           blanket_u[0]['z'][::-1],
                                           blanket_l[0]['r'],
                                           blanket_l[0]['z'],
                                           Zjoin=self.sf.Mpoint[1])
                self.blanket = wrap({'r':r,'z':z},{'r':rfw,'z':zfw})
                for loop in ['inner','outer']:
                    self.blanket.sort_z(loop,select='lower')
                bfill = self.blanket.fill(plot=False)
                self.segment['blanket_fw']  = bfill[0]
                self.segment['blanket']  = bfill[1]
            else:
                self.sf.get_Xpsi(select='lower')  # lower X-point
                div.place(debug=debug)
                self.segment = div.join(mc)   
                self.blanket = self.place_blanket(select='lower')[-1]
            self.vessel = self.vessel_fill()
            self.Rb = self.segment['first_wall']['r']
            self.Zb = self.segment['first_wall']['z']
        else:
            errtxt = 'set input mode \'eqdsk\',\'calc\''
            raise ValueError(errtxt)

        if mode != 'eqdsk':  # update eqdsk
            self.sf.xlim = self.Rb
            self.sf.ylim = self.Zb
            self.sf.nlim = len(self.sf.xlim)
            
        if plot:
            pl.plot(self.segment['first_wall']['r'],
                    self.segment['first_wall']['z'],
                    lw=1.75,color=0.5*np.ones(3))
            self.blanket.fill(plot=True,color=colors[0])
            self.vessel.fill(plot=True,color=colors[1])
            pl.axis('equal')
            pl.axis('off')
            
    def upper_lower(self,Ru,Zu,Rl,Zl,Zjoin=0):
        u_index = np.arange(len(Ru),dtype='int')
        iu_lower = u_index[Zu<Zjoin]
        iu_out,iu_in = iu_lower[0]-1,iu_lower[-1]+1
        l_index = np.arange(len(Rl),dtype='int')
        il_upper = l_index[Zl>Zjoin]
        il_out,il_in = il_upper[0]-1,il_upper[-1]+1
        R = np.append(Rl[:il_out],Ru[:iu_out][::-1])
        R = np.append(R,Ru[iu_in:][::-1])
        R = np.append(R,Rl[il_in:])
        Z = np.append(Zl[:il_out],Zu[:iu_out][::-1])
        Z = np.append(Z,Zu[iu_in:][::-1])
        Z = np.append(Z,Zl[il_in:])
        return R,Z
    
    def place_blanket(self,select='upper',store=True,plot=False):
        blanket = wrap(self.segment['first_wall'],self.segment['inner_loop'])
        blanket.sort_z('inner',select=select)
        blanket.offset('outer',self.setup.build['BB'],
                       ref_o=3/10*np.pi,dref=np.pi/3)
        bfill = blanket.fill(plot=plot,color=colors[0])
        if store:
            self.segment['blanket_fw']  = bfill[0]
            self.segment['blanket']  = bfill[1]
        return bfill,blanket

    def vessel_fill(self,gap=True):
        r,z = self.segment['blanket_fw']['r'],self.segment['blanket_fw']['z']
        loop = geom.Loop(r,z)
        r,z = loop.fill(dt=0.05)
        rb = np.append(r,r[0])
        zb = np.append(z,z[0])
        profile = Profile(self.setup.configuration,family='S',part='vv',
                          npoints=400,read_write=False)
        shp = Shape(profile,objective='L')
        shp.loop.set_l({'value':0.5,'lb':0.6,'ub':1.5})  # 1/tesion)
        shp.loop.xo['upper'] = {'value':0.33,'lb':0.5,'ub':1}  
        shp.loop.xo['lower'] = {'value':0.33,'lb':0.5,'ub':1}
        #shp.loop.oppvar.remove('flat')
        r,z = geom.rzSLine(rb,zb,200)  # sub-sample
        rup,zup = r[z>self.sf.Xpoint[1]],z[z>self.sf.Xpoint[1]]
        shp.add_bound({'r':rup,'z':zup},'internal')  # vessel inner bounds
        rd,zd = r[z<self.sf.Xpoint[1]],z[z<self.sf.Xpoint[1]]
        ro,zo = geom.offset(rd,zd,0.1)  # divertor offset
        shp.add_bound({'r':rd,'z':zd},'internal')  # vessel inner bounds
        shp.add_bound({'r':rd,'z':zd-0.25},'internal')  # gap below divertor
        #shp.plot_bounds()
        shp.minimise()
        #shp.loop.plot()
        x = profile.loop.draw()
        rin,zin = x['r'],x['z']
        loop = geom.Loop(rin,zin)
        r,z = loop.fill(dt=self.setup.build['VV'],ref_o=2/8*np.pi,dref=np.pi/6)
        shp.clear_bound()
        shp.add_bound({'r':r,'z':z},'internal')  # vessel outer bounds
        shp.minimise()
        x = profile.loop.draw()
        r,z = x['r'],x['z']
        if 'SX' in self.setup.configuration or gap == True:
            vv = wrap({'r':rin,'z':zin},{'r':r,'z':z})
        else:
            vv = wrap({'r':rb,'z':zb},{'r':r,'z':z})
        vv.sort_z('inner',select=self.sf.Xloc)
        vv.sort_z('outer',select=self.sf.Xloc)
        self.segment['vessel_inner'] = {'r':rin,'z':zin}
        self.segment['vessel_outer'] = {'r':r,'z':z}
        return vv
      
    def get_sol(self,plot=False):
        self.trim_sol(plot=plot)
        for leg in list(self.sf.legs)[2:]:
            L2D,L3D,Rsol,Zsol = self.sf.connection(leg,0)
            Ro,Zo = Rsol[-1],Zsol[-1]
            L2Dedge,L3Dedge = self.sf.connection(leg,-1)[:2]
            if leg not in self.setup.targets:
                self.setup.targets[leg] = {}
            Xi = self.sf.expansion([Ro],[Zo])
            graze,theta = np.zeros(self.sf.Nsol),np.zeros(self.sf.Nsol)
            pannel = self.sf.legs[leg]['pannel']
            for i in range(self.sf.Nsol):
                ro = self.sf.legs[leg]['R'][i][-1]
                zo = self.sf.legs[leg]['Z'][i][-1]
                graze[i] = self.sf.get_graze((ro,zo),pannel[i])
                theta[i] = self.sf.strike_point(Xi,graze[i])                             
            self.setup.targets[leg]['graze_deg'] = graze*180/np.pi
            self.setup.targets[leg]['theta_deg'] = theta*180/np.pi
            self.setup.targets[leg]['L2Do'] = L2D
            self.setup.targets[leg]['L3Do'] = L3D
            self.setup.targets[leg]['L2Dedge'] = L2Dedge
            self.setup.targets[leg]['L3Dedge'] = L3Dedge
            self.setup.targets[leg]['Ro'] = Ro
            self.setup.targets[leg]['Zo'] = Zo 
            self.setup.targets[leg]['Rsol'] = Rsol
            self.setup.targets[leg]['Zsol'] = Zsol 
        
    def trim_sol(self,color='k',plot=True):
        self.sf.sol()
        color = sns.color_palette('Set2',self.sf.nleg+5)
        #color = 'k'
        for c,leg in enumerate(self.sf.legs.keys()):
            if 'core' not in leg:
                Rsol = self.sf.legs[leg]['R']
                Zsol = self.sf.legs[leg]['Z']
                self.sf.legs[leg]['pannel'] = [[] for i in range(self.sf.Nsol)]
                for i in range(self.sf.Nsol):
                    if len(Rsol[i]) > 0:
                        R,Z = Rsol[i],Zsol[i]
                        for j in range(2):  # predict - correct
                            R,Z,pannel = self.trim(self.Rb,self.Zb,R,Z)
                        self.sf.legs[leg]['R'][i] = R  # update sf
                        self.sf.legs[leg]['Z'][i] = Z
                        self.sf.legs[leg]['pannel'][i] = pannel
                        if plot:
                            if color != 'k' and i > 0:
                                pl.plot(R,Z,'-',color=0.7*np.ones(3))  # color[c+3]
                                #pl.plot(R,Z,'-',color=color[c+3])  
                            elif color == 'k':
                                pl.plot(R,Z,'-',color='k',alpha=0.15)
                            else:
                                #pl.plot(R,Z,color=color[c])
                                pl.plot(R,Z,'--',color=[0.5,0.5,0.5])

    def crossed_lines(self,Ro,Zo,R1,Z1):
        index = np.zeros(2)
        dl = np.zeros(len(Ro))
        for i,(ro,zo) in enumerate(zip(Ro,Zo)):
            dl[i] = np.min((R1-ro)**2+(Z1-zo)**2)
        index[0] = np.argmin(dl)
        index[1] = np.argmin((R1-Ro[index[0]])**2+(Z1-Zo[index[0]])**2)
        return index
    
    def trim(self,Rloop,Zloop,R,Z):
        Rloop,Zloop = geom.order(Rloop,Zloop)
        L = geom.length(R,Z)
        index = np.append(np.diff(L)!=0,True)
        R,Z = R[index],Z[index]  # remove duplicates
        nRloop,nZloop,Rloop,Zloop = geom.normal(Rloop,Zloop)
        Rin,Zin = np.array([]),np.array([])
        for r,z in zip(R,Z):
            i = np.argmin((r-Rloop)**2+(z-Zloop)**2)
            dr = [Rloop[i]-r,Zloop[i]-z]  
            dn = [nRloop[i],nZloop[i]]
            if np.dot(dr,dn) > 0:
                Rin,Zin = np.append(Rin,r),np.append(Zin,z)
        i = np.argmin((Rin[-1]-R)**2+(Zin[-1]-Z)**2)
        Rin,Zin = R[:i+2],Z[:i+2] # extend past target
        i = np.argmin((Rin[-1]-R)**2+(Zin[-1]-Z)**2)  # sol crossing bndry
        jo = np.argmin((R[i]-Rloop)**2+(Z[i]-Zloop)**2)
        j = np.array([jo,jo+1])
        s = np.array([R[i],Z[i]])
        ds = np.array([R[i]-R[i-1],Z[i]-Z[i-1]])
        b = np.array([Rloop[j[0]],Zloop[j[0]]])
        bstep,db = self.get_bstep(s,ds,b,Rloop,Zloop,j)
        if bstep < 0:
            j = np.array([jo-1,jo])  # switch target pannel
            bstep,db = self.get_bstep(s,ds,b,Rloop,Zloop,j)
        step = np.cross(b-s,db)/np.cross(ds,db)
        intersect = s+step*ds  # predict - correct
        if step < 0:  # step back
            Rin,Zin = Rin[:-1],Zin[:-1]
        Rin,Zin = np.append(Rin,intersect[0]),np.append(Zin,intersect[1])
        return Rin,Zin,db
        
    def get_bstep(self,s,ds,b,Rloop,Zloop,j):
        db = np.array([Rloop[j[1]]-Rloop[j[0]],Zloop[j[1]]-Zloop[j[0]]])
        step = np.cross(b-s,ds)/np.cross(ds,db)
        return step,db

    def json(self,**kwargs):
        data = {}
        for loop in ['first_wall','blanket','vessel_inner','vessel_outer']:
            data[loop] = {}
            for var in self.segment[loop]:
                data[loop][var] = self.segment[loop][var].tolist()
        if 'tf' in kwargs:  # add tf profile
            tf = kwargs.get('tf')
            for loop,label in zip(['in','out'],['TF_inner','TF_outer']):
                data[label] = {}
                for var in tf.x[loop]:
                    data[label][var] = tf.x[loop][var].tolist()

        with open(self.datadir+'{}.json'.format(self.dataname),'w') as f:
            json.dump(data,f,indent=4,sort_keys=True)
            

class wrap(object):
    
    def __init__(self,inner_points,outer_points):
        self.loops = OrderedDict()
        self.loops['inner'] = {'points':inner_points}
        self.loops['outer'] = {'points':outer_points}
        
    def get_segment(self,loop):
        segment = self.loops[loop]['points']
        return segment['r'],segment['z']
        
    def set_segment(self,loop,r,z):
        self.loops[loop]['points'] = {'r':r,'z':z}
        
    def sort_z(self,loop,select='lower'):
        r,z = self.get_segment(loop)
        r,z = geom.order(r,z,anti=True)  # order points
        if select == 'lower':
            imin = np.argmin(z)  # locate minimum
            r = np.append(r[imin:],r[:imin])  # sort
            z = np.append(z[imin:],z[:imin])
        else:
            imax = np.argmax(z)  # locate minimum
            r = np.append(r[imax:],r[:imax])  # sort
            z = np.append(z[imax:],z[:imax])
        self.set_segment(loop,r,z)

    def offset(self,loop,dt,**kwargs):
        r,z = self.get_segment(loop)
        gloop = geom.Loop(r,z)
        r,z = gloop.fill(dt=dt,**kwargs)
        self.set_segment(loop,r,z)
           
    def interpolate(self):
        for loop in self.loops:
            r,z = self.get_segment(loop)
            r,z,l = geom.unique(r,z)
            interpolant = {'r':IUS(l,r),'z':IUS(l,z)}
            self.loops[loop]['fun'] = interpolant  

    def interp(self,loop,l):
        interpolant = self.loops[loop]['fun']
        return interpolant['r'](l),interpolant['z'](l)

    def sead(self,dl,N=500): # course search
        l = np.linspace(dl[0],dl[1],N)
        r,z = np.zeros((N,2)),np.zeros((N,2))
        for i,loop in enumerate(self.loops):
            r[:,i],z[:,i] = self.interp(loop,l)
        dr_min,i_in,i_out = np.max(r[:,1]),0,0
        for i,(rin,zin) in enumerate(zip(r[:,0],z[:,0])):
            dR = np.sqrt((r[:,1]-rin)**2+(z[:,1]-zin)**2)
            dr = np.min(dR)
            if dr < dr_min:
                dr_min = dr
                i_in = i
                i_out = np.argmin(dR)
        return l[i_in],l[i_out]
    
    def cross(self,L):
        r,z = np.zeros(2),np.zeros(2)
        for i,(loop,l) in enumerate(zip(self.loops,L)):
            r[i],z[i] = self.interp(loop,l)
        err = (r[0]-r[1])**2+(z[0]-z[1])**2
        return err
        
    def index(self,loop,l):
        rp,zp = self.interp(loop,l)  # point
        r,z = self.get_segment(loop)
        i = np.argmin((r-rp)**2+(z-zp)**2)
        return i
        
    def close_loop(self):
        for loop in self.loops:
            r,z = self.get_segment(loop)
            if (r[0]-r[-1])**2+(z[0]-z[-1])**2 != 0:
                r = np.append(r,r[0])
                z = np.append(z,z[0])
                self.set_segment(loop,r,z)
        
    def concentric(self,rin,zin,rout,zout):
        points = geom.inloop(rout,zout,rin,zin,side='out')
        if np.shape(points)[1] == 0:
            return True
        else:
            return False
        
    def fill(self,plot=False,color=colors[0]):  # minimization focused search
        rin,zin = self.get_segment('inner')
        rout,zout = self.get_segment('outer')        
        concentric = self.concentric(rin,zin,rout,zout)
        if concentric:
            self.close_loop()
            rin,zin = self.get_segment('inner')
            rout,zout = self.get_segment('outer') 
        self.interpolate()  # construct interpolators
        self.indx = {'inner':np.array([0,len(rin)],dtype=int),
                     'outer':np.array([0,len(rout)],dtype=int)}
        if not concentric:
            self.indx = {'inner':np.zeros(2,dtype=int),
                         'outer':np.zeros(2,dtype=int)}
            for i,dl in enumerate([[0,0.5],[0.5,1]]):  # low feild / high feild
                lo = self.sead(dl)
                L = minimize(self.cross,lo,method='L-BFGS-B',
                             bounds=([0,1],[0,1])).x 
                for loop,l in zip(self.loops,L):
                    self.indx[loop][i] = self.index(loop,l)

        r = np.append(rout[self.indx['outer'][0]:self.indx['outer'][1]],
                      rin[self.indx['inner'][0]:self.indx['inner'][1]][::-1])
        z = np.append(zout[self.indx['outer'][0]:self.indx['outer'][1]],
                      zin[self.indx['inner'][0]:self.indx['inner'][1]][::-1])
        self.patch = {'r':r,'z':z}
        r = np.append(np.append(rin[:self.indx['inner'][0]],
                      rout[self.indx['outer'][0]:self.indx['outer'][1]]),
                      rin[self.indx['inner'][1]:])
        z = np.append(np.append(zin[:self.indx['inner'][0]],
                      zout[self.indx['outer'][0]:self.indx['outer'][1]]),
                      zin[self.indx['inner'][1]:])
        self.segment = {'r':r,'z':z}
        if plot:
            self.plot(color=color)
        return self.segment,self.patch
            
    def plot(self,plot_loops=False,color=colors[0]):
        geom.polyfill(self.patch['r'],self.patch['z'],color=color)
        if plot_loops:
            pl.plot(self.segment['r'],self.segment['z'],color=0.75*np.ones(3))
            for loop in self.loops:
                r,z = self.get_segment(loop)
                pl.plot(r,z,'-')
    
    

        


