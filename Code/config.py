import pylab as pl
import numpy as np
import shelve
from scipy.interpolate import interp1d

class layout(object):
    
    def __init__(self, geom):
        self.geom = geom
        self.coil_keys = geom.coil.keys()
        self.flip = 1
        
    def set_keys(self,keys):
        self.coil_keys = keys
        
    def get_keys(self,keys):
        self.coil_keys = self.geom.coil.keys()
        
    def plasma_coils(self):
        coil = self.geom.coil
        for name in self.coil_keys:
            if 'plasma' in name:
                pl.plot(self.flip*coil[name]['r'], coil[name]['z'], 'bx')
        
    def sol(self):
        plasma = self.geom.plasma
        pl.plot(self.flip*plasma['sol']['r'], plasma['sol']['z'], 'r-', linewidth=2)

        
    def P6_support(self, coil):
        data = shelve.open(self.path+self.file+'_'+'sol')
        r_sol = data['data-1'][:,0]
        z_sol = data['data-1'][:,1]
        data.close()
        
        # upper pinch
        r = [r for r,z in zip(r_sol,z_sol) if r>9 and z>-7.5]
        z = [z for r,z in zip(r_sol,z_sol) if r>9 and z>-7.5]
        ind = z.index(min(z))
        r_min = r[ind]
        z_min = z[ind]
        
        # lower pinch
        r = [r for r,z in zip(r_sol,z_sol) if r>9 and z<-7.5]
        z = [z for r,z in zip(r_sol,z_sol) if r>9 and z<-7.5]
        ind = z.index(max(z))
        r_max = r[ind]
        z_max = z[ind]
        
        Lpinch = 0.9*(z_max-z_min)
        
        r = np.mean([r_min,r_max])
        z = np.mean([z_min,z_max])

        pl.plot(self.flip*r,z_min, 'rx')
        pl.plot(self.flip*r,z_max, 'rx')
        pl.plot(self.flip*r,z, 'ko')
        
        pl.plot(self.flip*[r,coil['P6']['r']],[z,coil['P6']['z']], 'k-')
        
    def plates(self):
        structure = self.geom.structure
        key = ['divertorL', 'divertorU']
        for k in key:
            if k in structure.keys():
                pl.plot(self.flip*structure[k]['r'], structure[k]['z'], 'go-')
    
    def radial_cord(self, r_in, z_in):
        cylin = {}        
        r = ((r_in-self.geom.rc)**2+(z_in-self.geom.zc)**2)**0.5
        theta = np.arctan2(z_in-self.geom.zc, r_in-self.geom.rc)
        theta = np.unwrap(theta)
        index = theta.argsort()
        r,theta = r[index],theta[index]                   
        theta_i = np.linspace(np.min(theta), np.max(theta), 100)
        fr_i = interp1d(theta, r, kind='linear')
        r_i = fr_i(theta_i)
        cylin['r'] = self.geom.rc+r_i*np.cos(theta_i)
        cylin['z'] = self.geom.zc+r_i*np.sin(theta_i)
        cylin['r'] = np.append(cylin['r'], cylin['r'][0])
        cylin['z'] = np.append(cylin['z'], cylin['z'][0])
        
        return cylin
            
    def TF(self):
        structure = self.geom.structure
        key = ['TFin', 'TFout']
        cylin = {}
        for k in key:  
            cylin[k] = {}
            cylin[k] = self.radial_cord(structure[k]['r'], structure[k]['z'])
            pl.plot(self.flip*cylin[k]['r'], cylin[k]['z'], 'k-')
 
        for i in range(len(cylin['TFin']['r'])-1):
            r_fill = [cylin['TFin']['r'][i], cylin['TFout']['r'][i], 
                      cylin['TFout']['r'][i+1], cylin['TFin']['r'][i+1]]
            z_fill = [cylin['TFin']['z'][i], cylin['TFout']['z'][i], 
                      cylin['TFout']['z'][i+1], cylin['TFin']['z'][i+1]]          
            pl.fill(r_fill,z_fill,facecolor='b',alpha=0.1, edgecolor='none')
    
    def FW(self):
        fw = self.geom.fw
        structure = self.geom.structure
        cylin_fw = self.radial_cord(fw['r'], fw['z'])
        cylin_tf = self.radial_cord(structure['TFin']['r'], structure['TFin']['z'])

        for i in range(len(cylin_fw['r'])-1):
            r_fill = [cylin_fw['r'][i], cylin_tf['r'][i], 
                      cylin_tf['r'][i+1], cylin_fw['r'][i+1]]
            z_fill = [cylin_fw['z'][i], cylin_tf['z'][i], 
                      cylin_tf['z'][i+1], cylin_fw['z'][i+1]]          
            pl.fill(r_fill,z_fill,facecolor='b',alpha=0.2, edgecolor='none')            
            
        pl.plot(self.flip*cylin_fw['r'], cylin_fw['z'], 'k-.')
            
    def plasma(self):
        rc, zc = self.geom.rc, self.geom.zc
        r,z = [],[]
        pl.plot(self.flip*rc, zc, 'go')
        for radius,theta in zip(self.geom.radius, self.geom.theta):
            r.append(rc+radius*np.cos(theta))
            z.append(zc+radius*np.sin(theta))
        pl.plot(self.flip*r, z, 'r--')
        
    def coil_fill(self):
        fig = pl.gcf()
        coil = self.geom.coil
        for name in self.coil_keys:
            if 'plasma' not in name:
                if coil[name]['I'] == 0:
                    color = 'k'
                elif coil[name]['I'] > 0:
                    color = 'r'
                else:
                    color = 'b'    
                circle=pl.Circle((self.flip*coil[name]['r'],
                                  coil[name]['z']),coil[name]['rc'],
                                 color=color)
                fig.gca().add_artist(circle)
            
    def coil_circ(self):
        coil = self.geom.coil
        theta = np.linspace(0,2*np.pi)
        for name in self.coil_keys:
            #pl.plot(self.flip*coil[name]['r'], coil[name]['z'], 'rx')  # coil centers
            x_circ = coil[name]['rc']*np.cos(theta)
            y_circ = coil[name]['rc']*np.sin(theta)
            if coil[name]['I'] > 0:
                line = 'r-'
            else:
                line = 'b-'
            pl.plot(self.flip*coil[name]['r']+x_circ, coil[name]['z']+y_circ, line) 
    
    def coil_sheild(self):
        coil = self.geom.coil
        theta = np.linspace(0,2*np.pi)
        for name in self.coil_keys: 
            if 'plasma' not in name:
                x_circ = (coil[name]['rc']+0.8)*np.cos(theta)
                y_circ = (coil[name]['rc']+0.8)*np.sin(theta)
                pl.plot(self.flip*coil[name]['r']+x_circ, coil[name]['z']+y_circ, 'k--') 
        
    def coil_label(self, coil_keys=[]):
        coil = self.geom.coil
        if not coil_keys:  # empty list
          coil_keys = coil.keys()
          
        for name in self.coil_keys:
            if 'plasma' not in name and name in coil_keys:
                ax = pl.gca()
                bbox_props = dict(boxstyle="round", fc="white", ec="k", lw=2)
                ax.text(coil[name]['r']+coil[name]['rc']*np.cos(np.pi/4),
                        coil[name]['z']+coil[name]['rc']*np.sin(np.pi/4), 
                        '  '+name,
                        horizontalalignment='left',verticalalignment='bottom',
                        bbox=bbox_props)
                        
    def coil_current(self,xlim,zlim):
        coil = self.geom.coil
        for name in self.coil_keys:
            if 'plasma' not in name and \
            (coil[name]['r']>xlim[0] and coil[name]['r']<xlim[1]) and \
            (coil[name]['z']>zlim[0] and coil[name]['z']<zlim[1]):
                ax = pl.gca()
                bbox_props = dict(boxstyle="round", fc="white", ec="y", lw=2)
                ax.text(coil[name]['r'],coil[name]['z']-0.5, 
                        '{:1.1f}'.format(coil[name]['I']/1e6)+'MA',
                        horizontalalignment='center',verticalalignment='top',
                        bbox=bbox_props)
                        
    def coil_force(self,xlim,zlim):
        coil = self.geom.coil
        for name in self.coil_keys:
            if 'plasma' not in name and \
            (coil[name]['r']>xlim[0] and coil[name]['r']<xlim[1]) and \
            (coil[name]['z']>zlim[0] and coil[name]['z']<zlim[1]):
                ax = pl.gca()
                bbox_props = dict(boxstyle="round", fc="white", ec="g", lw=2)
                ax.text(coil[name]['r'],coil[name]['z']+0.65, 
                        '{:1.1f}'.format(coil[name]['Fz']/1e6)+'MN',
                        horizontalalignment='center',verticalalignment='bottom',
                        bbox=bbox_props)
        
        


