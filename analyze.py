import analysis as al
import matplotlib.pyplot as plt
import numpy as np
import paper_plt
paper_plt.load_basic_config()

KP = 1.0
KE = 1.0

def main():
    NT = 64
    NX = 4
    NY = 8

    ens = np.fromfile('tmp.out.ens.dat', dtype=np.uint8).reshape(-1, NT, NX, NY)

    print(ens[0,0])
    print(ens[1,0])
    print(ens[2,0])
    print(ens[3,0])
    
    HT = np.fromfile('tmp.out.HT.dat', dtype=np.float64)
    HP = np.fromfile('tmp.out.HP.dat', dtype=np.float64)
    HE = np.fromfile('tmp.out.HE.dat', dtype=np.float64)

    MA, MB = measure_M(ens)

    fig, ax = plt.subplots(1,1)
    bins = np.linspace(0.0, 1.0, num=51, endpoint=True)
    ax.hist2d(MA, MB, bins=bins) #, range=[[-0.5, 1.5], [-0.5, 1.5]])
    ax.set_xlabel(r'$M_A$')
    ax.set_ylabel(r'$M_B$')
    ax.set_aspect(1)

    fig, axes = plt.subplots(2,1, sharex=True)
    axes[0].plot(MA, label='$M_A$')
    axes[0].legend()
    axes[1].plot(MB, label='$M_B$')
    axes[1].legend()

    fig, axes = plt.subplots(3, 2, figsize=(6,6), sharey='row')
    for (axl, axr), Hi, label in zip(axes, [HT, HP, HE], ['$H_T$', '$H_P$', '$H_E$']):
        xs = np.arange(len(Hi))
        axl.plot(xs[::100], Hi[::100], color='0.8')
        axl.plot(*al.bin_data(Hi, binsize=100), label=label)
        axl.legend()
        axr.hist(Hi, bins=30, orientation='horizontal')

    N_TRI = NX * NY / 4
    N_PET = NX * NY*3 / 8
    H_est = al.bootstrap(-HT*N_TRI - (KP*HP + KE*HE)*N_PET, Nboot=1000, f=al.rmean)
    print(f'Ground state energy: {H_est}')
        
    plt.show()

def measure_M(ens):
    coordsA = [(0,0), (0,4), (2,2)]
    coordsB = [(0,2), (2,0), (2,4)]
    MA = np.mean([ens[(...,*cA)] for cA in coordsA], axis=(0,2))
    MB = np.mean([ens[(...,*cB)] for cB in coordsB], axis=(0,2))
    # MA += 0.05*np.random.normal(size=MA.shape)
    # MB += 0.05*np.random.normal(size=MA.shape)
    assert len(MA.shape) == 1
    assert MA.shape == MB.shape
    return MA, MB

if __name__ == '__main__':
    main()
