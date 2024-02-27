import analysis as al
import matplotlib.pyplot as plt
import numpy as np
import paper_plt
paper_plt.load_basic_config()

def main():
    NT = 16
    NX = 4
    NY = 8
    ens = np.fromfile('tmp.out', dtype=np.uint8).reshape(-1, NT, NX, NY)
    MA, MB = measure_M(ens)
    fig, ax = plt.subplots(1,1)
    ax.hist2d(MA, MB, bins=30, range=[[-0.5, 1.5], [-0.5, 1.5]])
    ax.set_xlabel(r'$M_A$')
    ax.set_ylabel(r'$M_B$')
    ax.set_aspect(1)
    plt.show()

def measure_M(ens):
    coordsA = [(0,0), (0,4), (2,2)]
    coordsB = [(0,2), (2,0), (2,4)]
    MA = np.mean([ens[(...,*cA)] for cA in coordsA], axis=(0,2))
    MB = np.mean([ens[(...,*cB)] for cB in coordsB], axis=(0,2))
    assert len(MA.shape) == 1
    assert MA.shape == MB.shape
    return MA, MB

if __name__ == '__main__':
    main()
