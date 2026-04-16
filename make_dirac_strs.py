### Place Dirac string geometries at various separations on various lattice sizes

import numpy as np

def make_dirac_str(nrows, ncols, r, fname):
    erows = 2*nrows
    ecols = 4*ncols
    out = np.zeros((erows, ecols), dtype=np.uint8)
    row = 2*(nrows//2) + 1
    xi = (ecols - 4*r)//2
    xf = xi + 4*(r-1)
    assert xi >= 0 and xf < ecols
    out[row, xi:xf+1] = 1
    print(f'Saving to {fname}')
    out.tofile(fname)

def make_dirac_str2(nrows, ncols, r, fname):
    """Same geometry deformed dirac str for testing"""
    erows = 2*nrows
    ecols = 4*ncols
    out = np.zeros((erows, ecols), dtype=np.uint8)
    row = 2*(nrows//2) + 1
    xi = (ecols - 4*r)//2
    xf = xi + 4*(r-1)
    assert xi >= 0 and xf < ecols
    out[row-2, xi:xf+1] = 1
    out[row-1, xi-1] = 1
    out[row-1, xf+1] = 1
    print(f'Saving to {fname}')
    out.tofile(fname)

def make_init2(nrows, ncols, r, fname):
    """Initialization for second dirac str for testing"""
    erows = 2*nrows
    ecols = 4*ncols
    out = np.zeros((erows, ecols), dtype=np.uint8)
    row = 2*(nrows//2) + 1
    xi = (ecols - 4*r)//2
    xf = xi + 4*(r-1)
    assert xi >= 0 and xf < ecols
    out[row-1, xi-1:xf+1] = 1
    out[row-2, xi-1:xf+1] = 1
    print(f'Saving to {fname}')
    out.tofile(fname)

def main():
    L = 64
    for r in range(8, 48, 4):
        make_dirac_str(L, L, r, f'geoms/{L}_{L}_r{r}_dirac_str.dat')
        make_dirac_str2(L, L, r, f'geoms/{L}_{L}_r{r}_dirac_str2.dat')
        make_init2(L, L, r, f'geoms/{L}_{L}_r{r}_init2.dat')
    L = 8
    for r in [4, 6]:
        make_dirac_str(L, L, r, f'geoms/{L}_{L}_r{r}_dirac_str.dat')
        make_dirac_str2(L, L, r, f'geoms/{L}_{L}_r{r}_dirac_str2.dat')
        make_init2(L, L, r, f'geoms/{L}_{L}_r{r}_init2.dat')

if __name__ == '__main__':
    main()
