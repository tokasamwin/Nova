import h5py
import numpy as np

with h5py.File("mytestfile.hdf5", "w") as f:
    dset = f.create_dataset("mydataset", (100,), dtype='float')
    dset[...] = np.arange(100)
    
    print(dset[:10])