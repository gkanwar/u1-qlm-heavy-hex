import analysis as al
import matplotlib.pyplot as plt
import numpy as np
import paper_plt
paper_plt.load_basic_config()

Ts = np.array([16, 32, 64, 128])
dts = 0.01 * np.arange(1, 5+1)
KP = 1.0
def main():
    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(6,3))
    axes[0].set_title(r'$\left< H_T \right>$')
    axes[1].set_title(r'$\left< H_P \right>$')
    axes[2].set_title(r'$\left< H \right>$')
    style = dict(marker='o', fillstyle='none', markersize=4)
    for dt in dts:
        LTs = []
        HTs = []
        HPs = []
        Hs = []
        for T in Ts:
            prefix = f'KP{KP:.2f}_T{T}_dt{dt:.2f}'
            HT = np.fromfile(f'{prefix}.HT.dat', dtype=np.int64).reshape(-1, 2)[:,1]
            HP = np.fromfile(f'{prefix}.HP.dat', dtype=np.int64).reshape(-1, 2)[:,1]
            HT = al.bin_data(-HT[int(0.1*len(HT)):], binsize=100)[1]
            HP = al.bin_data(-KP*HP[int(0.1*len(HP)):], binsize=100)[1]
            HTs.append(al.bootstrap(HT, Nboot=1000, f=al.rmean))
            HPs.append(al.bootstrap(HP, Nboot=1000, f=al.rmean))
            Hs.append(al.bootstrap(HT+HP, Nboot=1000, f=al.rmean))
            LTs.append(dt*T)
        LTs = np.array(LTs)
        HTs = np.stack(HTs, axis=1)
        HPs = np.stack(HPs, axis=1)
        Hs = np.stack(Hs, axis=1)

        al.add_errorbar(HTs, xs=1/LTs, ax=axes[0], label=f'$dt={dt:.2f}$', **style)
        al.add_errorbar(HPs, xs=1/LTs, ax=axes[1], label=f'$dt={dt:.2f}$', **style)
        al.add_errorbar(Hs, xs=1/LTs, ax=axes[2], label=f'$dt={dt:.2f}$', **style)

    for ax in axes:
        ax.legend()
        ax.set_xlabel(r'Temp')

    fig.suptitle('6T6P6O')
    fig.set_tight_layout(True)
    fig.savefig('figs/analyze_sweep.pdf')

if __name__ == '__main__':
    main()
