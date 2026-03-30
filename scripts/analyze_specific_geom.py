import analysis as al
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import paper_plt
import scipy as sp
import struct
paper_plt.load_basic_config()

def x_over_tanh(x, c):
    return np.where(
        np.isclose(c * x, 0.0),
        np.ones_like(x) / c, # small arg limit
        x / np.tanh(c * x)
    )

def bin_search(xs, ys, y0):
    """Binary search and lin interp to x0 st f(x0) = y0"""
    i = 0
    j = len(ys)-1
    while j > i+1:
        i_mid = (j+i)//2
        if ys[i_mid] > y0:
            j = i_mid
        else:
            i = i_mid
    alpha = (y0-ys[i])/(ys[j]-ys[i])
    return (1-alpha)*xs[i] + alpha*xs[j]

def analyze(KP_str, geom, E0):
    Ts = np.array([16, 32, 64, 128])
    dts = 0.01 * np.arange(1, 5+1)
    data_dir = f'data_new/{geom}'
    KT = 1.0 # fixed in hamiltonian
    KP = float(KP_str)
    ed_fname = f'{data_dir}/ED_thermal_expt_KP1.0.txt'
    if os.path.exists(ed_fname):
        ed_data = np.loadtxt(ed_fname)
        print(f'{ed_data.shape=}')
    else:
        ed_data = None
    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(6,3))
    axes[0].set_title(r'$\left< H_T \right>$')
    axes[1].set_title(r'$\left< H_P \right>$')
    axes[2].set_title(r'$\left< H \right>$')
    style = dict(marker='o', fillstyle='none', markersize=4, linestyle='')
    # style = dict(linestyle='', linewidth=0.5, capsize=1.5)
    cmap = plt.get_cmap('viridis')
    for i,dt in enumerate(dts):
        LTs = []
        HTs = []
        HPs = []
        Hs = []
        for T in Ts:
            prefix = f'{data_dir}/KP{KP:.2f}_T{T}_dt{dt:.2f}'
            # load energies and analyze
            HT = np.fromfile(f'{prefix}.HT.dat', dtype=np.int64).reshape(-1, 2)
            HP = np.fromfile(f'{prefix}.HP.dat', dtype=np.int64).reshape(-1, 2)
            # instead of this transfer matrix approach we estimate assuming a step
            # is proportional to (1 - dt H), so we extract the components of the transition
            # part proportional to dt KP and dt KT respectively.
            HT = -(KT*HT[...,0]*np.tanh(dt * KT) + HT[...,1]*x_over_tanh(KT, dt)) / T
            HP = -(KP*HP[...,0]*np.tanh(dt * KP) + HP[...,1]*x_over_tanh(KP, dt)) / T
            # HT = -HT[...,1]/(dt*T)
            # HP = -HP[...,1]/(dt*T)
            HT = al.bin_data(HT[int(0.1*len(HT)):], binsize=100)[1]
            HP = al.bin_data(HP[int(0.1*len(HP)):], binsize=100)[1]
            HTs.append(al.bootstrap(HT, Nboot=1000, f=al.rmean))
            HPs.append(al.bootstrap(HP, Nboot=1000, f=al.rmean))
            Hs.append(al.bootstrap(HT+HP, Nboot=1000, f=al.rmean))
            LTs.append(dt*T)
        LTs = np.array(LTs)
        HTs = np.stack(HTs, axis=1)
        HPs = np.stack(HPs, axis=1)
        Hs = np.stack(Hs, axis=1)

        al.add_errorbar(HTs, xs=1/LTs, ax=axes[0], label=f'$dt={dt:.2f}$', color=cmap(i/len(dts)), **style)
        al.add_errorbar(HPs, xs=1/LTs, ax=axes[1], label=f'$dt={dt:.2f}$', color=cmap(i/len(dts)), **style)
        al.add_errorbar(Hs, xs=1/LTs, ax=axes[2], label=f'$dt={dt:.2f}$', color=cmap(i/len(dts)), **style)

        # spline interpolate to target T at second-smallest dt
        if i == 1 and E0 is not None:
            H0 = -KP*E0
            T0 = bin_search(np.flip(1/LTs), np.flip(Hs[0]), H0)
            axes[2].axhline(H0, color='k', linestyle='--', label=f'$H_0={H0:.2f}$')
            axes[2].axvline(T0, color='k', linestyle='--', label=f'$T_0={T0:.2f}$')

    if ed_data is not None:
        axes[0].plot(ed_data[:,0], ed_data[:,1], color='k', label='ED', zorder=3)
        axes[1].plot(ed_data[:,0], ed_data[:,2], color='k', label='ED', zorder=3)
        axes[2].plot(ed_data[:,0], ed_data[:,3], color='k', label='ED', zorder=3)
    for ax in axes:
        ax.legend(loc='lower right', fontsize=7)
        ax.set_xlabel(r'Temp')

    fig.suptitle(rf'{geom} ($K_P = {KP_str}$)')
    fig.set_tight_layout(True)
    os.makedirs(f'{data_dir}/figs', exist_ok=True)
    fig.savefig(f'{data_dir}/figs/analyze_sweep_KP{KP_str}.pdf')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--KP', type=str, required=True)
    parser.add_argument('--geom', type=str, required=True)
    parser.add_argument('--E0', type=int)
    args = parser.parse_args()
    # for KP in ['0.40', '0.70', '2.00']:
    analyze(args.KP, args.geom, args.E0)

if __name__ == '__main__':
    main()
