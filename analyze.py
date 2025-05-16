import argparse
import analysis as al
import matplotlib.pyplot as plt
import numpy as np
import os
import paper_plt
paper_plt.load_basic_config()
import struct

def x_over_tanh(x, c):
    return np.where(
        np.isclose(c * x, 0.0),
        np.ones_like(x) / c, # small arg limit
        x / np.tanh(c * x)
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, required=True)
    args = parser.parse_args()
    prefix = args.prefix
    prefix_dir, prefix_tail = os.path.split(prefix)
    figs_dir = os.path.join(prefix_dir, 'figs')
    os.makedirs(figs_dir, exist_ok=True)
    figs_prefix = os.path.join(figs_dir, prefix_tail)
    
    with open(f'{prefix}.meta.dat', 'rb') as f:
        meta_bytes = f.read()
    header = 'dddiiiI'
    h_size = struct.calcsize(header)
    (dt, KP, KE, NT, NX, NY, seed) = struct.unpack(header, meta_bytes[:h_size])
    print(f'Couplings: {dt=} {KP=} {KE=}')
    print(f'Shape: {NT=} {NX=} {NY=}')
    print(f'Seed: {seed=}')
    geom = np.frombuffer(meta_bytes[h_size:], dtype=np.uint8).reshape(NX, NY)
    # ens = np.fromfile(f'{prefix}.ens.dat', dtype=np.uint8).reshape(-1, NT, NX, NY)
    mx = np.fromfile(f'{prefix}.Mx.dat', dtype=np.float64).reshape(NX, NY)
    Ex = np.fromfile(f'{prefix}.Ex.dat', dtype=np.float64).reshape(NX, NY)
    HTx, HPx, HEx = np.fromfile(f'{prefix}.Hx.dat', dtype=np.float64).reshape(3, 2, NX, NY)
    free_tri_mask = np.zeros_like(geom)
    free_tri_mask[::2,::2][geom[::2,::2] == 0xff] = 1
    free_pet_even_mask = np.zeros_like(geom)
    free_pet_odd_mask = np.zeros_like(geom)
    free_pet_even_mask[:,1::2][geom[:,1::2] == 0xff] = 1
    free_pet_odd_mask[1::2,:][geom[1::2,:] == 0xff] = 1
    free_pet_mask = free_pet_even_mask | free_pet_odd_mask
    N_FREE_TRI = np.sum(free_tri_mask)
    N_FREE_PET = np.sum(free_pet_mask)
    print(f'Geom:\n{geom}')
    print(f'{N_FREE_TRI=} {N_FREE_PET=}')

    # print(ens[0,0])
    # print(ens[1,0])
    # print(ens[2,0])
    # print(ens[3,0])
    
    HT = np.fromfile(f'{prefix}.HT.dat', dtype=np.int64).reshape(-1, 2)
    HP = np.fromfile(f'{prefix}.HP.dat', dtype=np.int64).reshape(-1, 2)
    HE = np.fromfile(f'{prefix}.HE.dat', dtype=np.int64).reshape(-1, 2)
    assert HT.shape[0] == HP.shape[0] == HE.shape[0]
    print(f'{HT.shape=}')
    # DB conventions
    KT = 4
    KP *= 2
    # HT = (KT*HT[...,0]*np.tanh(dt * KT) + HT[...,1]*x_over_tanh(KT, dt)) / NT
    # HP = (KP*HP[...,0]*np.tanh(dt * KP) + HP[...,1]*x_over_tanh(KP, dt)) / NT
    HT = -KT*HT[...,1] / NT
    HP = -KP*HP[...,1] / NT
    HE = (KE * (HE[...,0] - HE[...,1])) / NT

    # HTx = (KT*HTx[0]*np.tanh(dt * KT) + HTx[1]*x_over_tanh(KT, dt)) / NT
    # HPx = (KP*HPx[0]*np.tanh(dt * KP) + HPx[1]*x_over_tanh(KP, dt)) / NT
    HTx = -KT*HTx[1] / NT
    HPx = -KP*HPx[1]
    HHEx = (HEx[0] - HEx[1]) / NT
    HEx = (KE * (HEx[0] - HEx[1])) / NT
    Hx = HTx + HPx + HEx
    HTx[free_tri_mask != 1] = float('nan')
    HPx[free_pet_mask != 1] = float('nan')
    HEx[free_pet_mask != 1] = float('nan')
    HHEx[free_pet_mask != 1] = float('nan')
    Hx[(free_tri_mask != 1) & (free_pet_mask != 1)] = float('nan')

    cmap = plt.get_cmap('PiYG').copy()
    cmap.set_bad(color='w')
    cmap2 = plt.get_cmap('PuOr').copy()
    cmap2.set_bad(color='w', alpha=0.0)

    MTx = mx - 0.5
    MTx[free_tri_mask != 1] = float('nan')
    MPx = mx - 0.5
    MPx[free_pet_mask != 1] = float('nan')
    Ex[free_pet_mask != 1] = float('nan')

    fig, axes = plt.subplots(
        2, 2, layout='compressed', figsize=(10,6), squeeze=False)
    ax = axes[0,0]
    # cax = axes[0,1]
    ax.set_title(r'$M_T(x)$')
    vmax = np.nanmax(np.abs(MTx))
    cs = ax.imshow(MTx, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    # fig.colorbar(cs, cax=cax)
    fig.colorbar(cs, ax=ax)
    ax = axes[0,1]
    # cax = axes[0,3]
    ax.set_title(r'$M_P(x)$')
    vmax = np.nanmax(np.abs(MPx))
    cs = ax.imshow(MPx, cmap=cmap2, vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    # fig.colorbar(cs, cax=cax)
    fig.colorbar(cs, ax=ax)

    # stacked images
    ax = axes[1,0]
    ax.set_title(r'$M_T(x)$, $M_P(x)$')
    vmax = np.nanmax(np.abs(MTx))
    ax.imshow(MTx, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation='nearest')
    vmax = np.nanmax(np.abs(MPx))
    ax.imshow(MPx, cmap=cmap2, vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_aspect(1)
    # Ex
    ax = axes[1,1]
    ax.set_title(r'$E_P(x)$')
    vmax = np.nanmax(np.abs(Ex))
    cs = ax.imshow(Ex, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation='nearest')
    fig.colorbar(cs, ax=ax)

    fig.suptitle(prefix)
    fig.savefig(f'{figs_prefix}.Mx.pdf', dpi=600)

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad(color='k')

    fig, axes = plt.subplots(
        2, 4, layout='compressed', figsize=(10,6),
        gridspec_kw=dict(width_ratios=[0.45, 0.02, 0.45, 0.02]))
    ax = axes[0,0]
    cax = axes[0,1]
    ax.set_title(r'$H_T(x)$')
    cs = ax.imshow(HTx, cmap=cmap, interpolation='nearest') # vmin=-0.5, vmax=0.5, 
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    fig.colorbar(cs, cax=cax)
    ax = axes[0,2]
    cax = axes[0,3]
    ax.set_title(r'$H_P(x)$')
    cs = ax.imshow(HPx, cmap=cmap, interpolation='nearest') # vmin=-0.5, vmax=0.5, 
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    fig.colorbar(cs, cax=cax)
    ax = axes[1,0]
    cax = axes[1,1]
    ax.set_title(r'$H_E(x)/K_E$')
    cs = ax.imshow(HHEx, cmap=cmap, interpolation='nearest') # vmin=-0.5, vmax=0.5, 
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    fig.colorbar(cs, cax=cax)
    ax = axes[1,2]
    cax = axes[1,3]
    ax.set_title(r'$H(x)$')
    cs = ax.imshow(Hx, cmap=cmap, interpolation='nearest') # vmin=-0.5, vmax=0.5, 
    ax.set_aspect(1)
    bounds = np.nonzero((geom != 0xff) & (geom != 0xaa))
    values = geom[bounds]
    ax.scatter(bounds[1], bounds[0], c=['r' if v == 0 else 'b' for v in values], marker='o')
    fig.colorbar(cs, cax=cax)
    fig.suptitle(prefix)
    fig.savefig(f'{figs_prefix}.Hx.pdf', dpi=600)

    # fig, ax = plt.subplots(1,1)
    # ens_f = ens.astype(np.float64)
    # ens_f[:,:,geom != 0xff] = float('nan')
    # cs = ax.imshow(ens_f.mean(axis=(0, 1)) - 0.5, vmin=-0.5, vmax=0.5, cmap=cmap, interpolation='nearest')
    # ax.set_aspect(1)
    # fig.colorbar(cs)
    # fig.savefig(f'{figs_prefix}.Mx_ens.pdf', dpi=600)

    # TODO: extract MT, MP from Mx data instead of ens
    # MA, MB, MP = measure_M(ens, geom)
    # MA, MB, MP = split_M(m, geom)
    # MT = (MA + MB)/2

    MT = np.fromfile(f'{prefix}.MT.dat', dtype=np.int64)
    MP = np.fromfile(f'{prefix}.MP.dat', dtype=np.int64)
    MT = MT / float(N_FREE_TRI * NT) - 0.5
    MP = MP / float(N_FREE_PET * NT) - 0.5

    fig, ax = plt.subplots(1,1)
    bins = np.linspace(-0.5, 0.5, num=51, endpoint=True)
    ax.hist2d(MT, MP, bins=bins) #, range=[[-0.5, 1.5], [-0.5, 1.5]])
    ax.set_xlabel(r'$M_T$')
    ax.set_ylabel(r'$M_P$')
    ax.set_aspect(1)
    fig.savefig(f'{figs_prefix}.M_hist.pdf')

    # fig, axes = plt.subplots(4,1, sharex=True)
    # axes[0].plot(MA, label='$M_A$')
    # axes[0].legend()
    # axes[1].plot(MB, label='$M_B$')
    # axes[1].legend()
    # axes[2].plot(MT, label='$M_T$')
    # axes[2].legend()
    # axes[3].plot(MP, label='$M_P$')
    # axes[3].legend()
    # fig.savefig(f'{figs_prefix}.M_trace.pdf')

    fig, axes = plt.subplots(3, 2, figsize=(6,6), sharey='row')
    for (axl, axr), Hi, label in zip(axes, [HT, HP, HE], ['$H_T$', '$H_P$', '$H_E$']):
        xs = np.arange(len(Hi))
        axl.plot(xs[::10], Hi[::10], color='0.8')
        axl.plot(*al.bin_data(Hi, binsize=10), label=label)
        axl.legend()
        axr.hist(Hi, bins=30, orientation='horizontal')
    fig.savefig(f'{figs_prefix}.H_trace.pdf')

    # N_TRI = NX * NY / 4
    # N_PET = NX * NY*3 / 8
    # H_est = al.bootstrap(-HT - HP - HE, Nboot=1000, f=al.rmean)
    # DB convention
    H_est = al.bootstrap(-HT - HP - HE, Nboot=1000, f=al.rmean)
    print(f'Ground state energy: {H_est}')
        
    # plt.show()

def get_sublattice(x, y):
    if (x + y) % 2 == 1:
        return 'P'
    if (x + y) % 4 == 2:
        return 'A'
    else:
        assert (x+y) % 4 == 0
        return 'B'

def split_M(mag, geom):
    MA = 0
    MB = 0
    MP = 0
    nA = 0
    nB = 0
    nP = 0
    for x in range(geom.shape[0]):
        for y in range(geom.shape[1]):
            if geom[x,y] != 0xff:
                continue
            sub = get_sublattice(x,y)
            m = mag[x,y]
            if sub == 'A':
                MA += m
                nA += 1
            elif sub == 'B':
                MB += m
                nB += 1
            else:
                assert sub == 'P'
                MP += m
                nP += 1
    MA /= nA
    MB /= nB
    MP /= nP
    MA -= 0.5
    MB -= 0.5
    MP -= 0.5
    return MA, MB, MP

def measure_M(ens, geom):
    MA = np.zeros(ens.shape[0], dtype=np.float64)
    MB = np.zeros(ens.shape[0], dtype=np.float64)
    MP = np.zeros(ens.shape[0], dtype=np.float64)
    nA = 0
    nB = 0
    nP = 0
    for x in range(geom.shape[0]):
        for y in range(geom.shape[1]):
            if geom[x,y] != 0xff:
                continue
            sub = get_sublattice(x,y)
            m = np.mean(ens[:,:,x,y], axis=1)
            if sub == 'A':
                MA += m
                nA += 1
            elif sub == 'B':
                MB += m
                nB += 1
            else:
                assert sub == 'P'
                MP += m
                nP += 1
    MA /= nA
    MB /= nB
    MP /= nP
    MA -= 0.5
    MB -= 0.5
    MP -= 0.5
    return MA, MB, MP

if __name__ == '__main__':
    main()
