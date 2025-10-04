### Place Dirac string geometries at various separations on various lattice sizes

import numpy as np

def make_dirac_str(nrows, ncols, r):
    erows = 2*nrows
    ecols = 4*ncols
    out = np.zeros((erows, ecols), dtype=np.uint8)
    row = 2*(nrows//2) + 1
    
