    def trim_legs(self):
        for leg in self.legs.keys():
            for layer_index in range(self.legs[leg]['i']):
                if 'core' not in leg:
                    R,Z = self.snip(leg,layer_index)
                    self.legs[leg]['R'][layer_index] = R
                    self.legs[leg]['Z'][layer_index] = Z
