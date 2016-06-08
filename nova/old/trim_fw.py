    def trim_fw(self,R,Z):
        r,z = R[0],Z[0]
        r1,r2 = r[r>self.Xpoint[0]],r[r<=self.Xpoint[0]]
        z1,z2 = z[r>self.Xpoint[0]],z[r<=self.Xpoint[0]]
        arg1 = np.argmin((r1-self.Xpoint[0])**2+(z1-self.Xpoint[1])**2)
        arg2 = np.argmin((r2-self.Xpoint[0])**2+(z2-self.Xpoint[1])**2)
        r = np.append(r1[arg1:],r2[:arg2+1])
        z = np.append(z1[arg1:],z2[:arg2+1])
        return r,z
