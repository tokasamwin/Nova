    def inflate_divertor(self,Rb,Zb):
        #Rb,Zb = self.midplane_loop(Rb,Zb)  # boundary clockwise from LFS
        Rb,Zb = Rb[::-1],Zb[::-1]  # flip boundary (to clockwise)
        Rb,Zb = self.rzSLine(Rb,Zb,Np=len(Rb))  # smooth
        index = np.zeros(3)
        for i,leg in enumerate(['inner','outer']):  # targets
            index[i] = np.argmin((self.targets[leg]['Ro']-Rb)**2+
                               (self.targets[leg]['Zo']-Zb)**2)
        iplus = np.argmin((self.sf.xo[0]-Rb[index[1]:])**2+
                          (self.sf.xo[1]-Zb[index[1]:])**2)
        Rplus,Zplus = Rb[index[1]:][iplus],Zb[index[1]:][iplus]
        index[i+1] = np.argmin((Rplus-Rb)**2+(Zplus-Zb)**2)  # X-point end
        R,Z = self.offset(Rb,Zb,self.conf.inflate)  # offset diveror rejoin
        for i,start,end,trans in zip(range(2),[0.1,0.1],[0.9,0.9],[0.4,0.4]):
            w = self.window(int(index[i+1]-index[i]),start,end,trans)
            Rb[index[i]:index[i+1]] = w*R[index[i]:index[i+1]]+(1-w)\
            *Rb[index[i]:index[i+1]]
            Zb[index[i]:index[i+1]] = w*Z[index[i]:index[i+1]]+(1-w)\
            *Zb[index[i]:index[i+1]]
        return Rb,Zb
