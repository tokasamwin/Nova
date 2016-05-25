    def dotTF(self,xCoil):
        dsum = 0
        for Side in ['in','out']:
            dsum += self.dot_diffrence(xCoil,Side) 
        print(dsum)    
        Zfrac = 0.9  # fraction of plasma covered by straight TF
        top = xCoil[-1]+xCoil[0]  # top of straight section
        bottom = xCoil[-1]-xCoil[0]  # bottom of straight section
        Zplasma = np.max(self.Zp)-np.min(self.Zp)
        #if np.max(self.Zp)-top > (1-Zfrac)/2*Zplasma:  # straigh at top
        #    dsum -= np.max(self.Zp)-top-(1-Zfrac)/2*Zplasma
        #if bottom-np.min(self.Zp) > (1-Zfrac)/2*Zplasma:  # straight at bottom
        #    dsum -= bottom-np.min(self.Zp)-(1-Zfrac)/2*Zplasma
        '''
        Zplasma = np.max(self.Zp)-self.sf.Mpoint[1]
        Zstraight = xCoil[0]+self.xo[1]-self.sf.Mpoint[1]#-xCoil[4])  # xCoil[0]-(self.xo[1]-xCoil[4])
        if 0.65*Zplasma > Zstraight:
            dsum -= 0.65*Zplasma-Zstraight  
        '''

        return dsum
