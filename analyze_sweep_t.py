import analysis as al
import argparse 
import matplotlib.pyplot as plt
import numpy as np
import paper_plt
paper_plt.load_basic_config()
import scipy as sp
import scipy.optimize
import struct

def x_over_tanh(x, c):
    return np.where(
        np.isclose(c * x, 0.0),
        np.ones_like(x) / c, # small arg limit
        x / np.tanh(c * x)
    )

def load(prefix):
    with open(f'{prefix}.meta.dat', 'rb') as f:
        meta_bytes = f.read()
    header = 'dddiiiI'
    h_size = struct.calcsize(header)
    (dt, KP, KE, NT, NX, NY, seed) = struct.unpack(header, meta_bytes[:h_size])
    print(f'Couplings: {dt=} {KP=} {KE=}')
    print(f'Shape: {NT=} {NX=} {NY=}')
    print(f'Seed: {seed=}')
    geom = np.frombuffer(meta_bytes[h_size:], dtype=np.uint8).reshape(NX, NY)
    ens = np.fromfile(f'{prefix}.ens.dat', dtype=np.uint8).reshape(-1, NT, NX, NY)
    N_FREE_TRI = np.sum(geom[::2,::2] == 0xff)
    N_FREE_PET = np.sum(geom[1::2,:] == 0xff) + np.sum(geom[:,1::2] == 0xff)
    print(f'Geom:\n{geom}')
    print(f'{N_FREE_TRI=} {N_FREE_PET=}')

    print(ens[0,0])
    print(ens[1,0])
    print(ens[2,0])
    print(ens[3,0])
    
    HT = np.fromfile(f'{prefix}.HT.dat', dtype=np.int32).reshape(-1, 2)
    HP = np.fromfile(f'{prefix}.HP.dat', dtype=np.int32).reshape(-1, 2)
    HE = np.fromfile(f'{prefix}.HE.dat', dtype=np.int32).reshape(-1, 2)
    assert HT.shape[0] == HP.shape[0] == HE.shape[0]
    print(f'{HT.shape=}')
    # DB conventions
    KT = 4
    KP *= 2
    HT = (KT*HT[...,0] * np.tanh(dt * KT) + HT[...,1]*x_over_tanh(KT, dt)) / NT
    HP = (KP*HP[...,0]*np.tanh(dt * KP) + HP[...,1]*x_over_tanh(KP, dt)) / NT
    HE = (KE * (HE[...,0] - HE[...,1])) / NT

    return dict(ens=ens, HT=HT, HP=HP, HE=HE)

def fit_dt2(xs, ys):
    def f(dt, A, B):
        return A + B*dt**2
    A0 = ys[0,0]
    B0 = 0
    p0 = [A0, B0]
    popt, pcov = sp.optimize.curve_fit(f, xs, ys[0], sigma=ys[1], p0=p0)
    fopt = lambda xs: f(xs, *popt)
    yfit = fopt(xs)
    resid = ys[0] - yfit
    chisq = np.sum(resid**2 / ys[1]**2)
    chisq_per_dof = chisq / (ys.shape[-1] - len(popt))
    min_x = np.min(xs)
    max_x = np.max(xs)
    range_x = max_x - min_x
    fit_xs = np.linspace(min_x - 0.02*range_x, max_x + 0.02*range_x, num=100, endpoint=True)
    fit_ys = fopt(fit_xs)
    return dict(
        popt=popt, pcov=pcov, fopt=fopt,
        chisq=chisq, chisq_per_dof=chisq_per_dof,
        trace=(fit_xs, fit_ys)
    )

def do_fit(traces, dt_max=None):
    xs = []
    ys = []
    for (xs_i, ys_i) in traces.values():
        xs.append(xs_i)
        ys.append(ys_i)
    xs = np.concatenate(xs, axis=-1)
    ys = np.concatenate(ys, axis=-1)
    if dt_max is not None:
        inds = np.nonzero(xs <= dt_max)[0]
        xs = xs[inds]
        ys = ys[:,inds]
    print(f'Fitting {xs=} {ys=}')
    return fit_dt2(xs, ys)

ALL_NT = [64, 80, 96, 128]
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, required=True)
    parser.add_argument('--dts', type=float, nargs='+', required=True)
    parser.add_argument('--ED', type=float)
    parser.add_argument('--fit_dt_max', type=float)
    args = parser.parse_args()

    traces = {}
    # all_data = {}
    for base_dt in args.dts:
        xs = []
        # data = []
        trace = []
        for NT in ALL_NT:
            dt_str = f'{base_dt*NT/ALL_NT[0]:.4f}'
            if dt_str[0] == "0": # strip leading 0
                dt_str = dt_str[1:]
            dt = float(dt_str)
            prefix = args.prefix + f'.T{NT}_dt{dt_str}'
            res = load(prefix)
            # H_est = al.bootstrap(-HT - HP - HE, Nboot=1000, f=al.rmean)
            # DB convention
            H = -res['HT'] - res['HP'] - res['HE']
            H_est = al.bootstrap(al.bin_data(H, binsize=10)[1], Nboot=1000, f=al.rmean)
            print(f'Ground state energy: {H_est}')
            # data.append(H)
            trace.append(H_est)
            xs.append(dt)
        xs = np.array(xs)
        # data = np.transpose(data)
        trace = np.transpose(trace)
        T = base_dt * ALL_NT[0]
        traces[T] = (xs, trace)
        # all_data[T] = (xs, data)

    # def fit_E_boot(inds):
    #     xs = []
    #     ys = []
    #     for xs_i,data_i in all_data.items():
    #         H_est_i = al.bootstrap(data_i[inds], Nboot=100, f=al.rmean)
    #         xs.append(xs_i)
    #         ys.append(H_est_i)
    #     xs = np.concatenate(xs, axis=-1)
    #     ys = np.concatenate(ys, axis=-1)
    #     inds = np.nonzero(xs <= args.fit_dt_max)[0]
    #     xs = xs[inds]
    #     ys = ys[inds]
    #     fit_t2(xs, ys)
    res = do_fit(traces, dt_max=args.fit_dt_max)
    print(f'Fit params = {res["popt"]}')
    print(f'Fit chisq/dof = {res["chisq_per_dof"]}')
    print(f'Fit err = {np.sqrt(res["pcov"][0,0])}')

    fig, ax = plt.subplots(1,1, figsize=(6,4))
    style = dict(marker='o', markersize=5, fillstyle='none', capsize=2, linestyle='')
    max_x = 0.0
    for T,(xs,trace) in traces.items():
        max_x = max(np.max(xs), max_x)
        al.add_errorbar(trace, xs=xs, ax=ax, label=rf'$E$ [$T={T}$]', **style)
    ax.plot(*res['trace'], color='k')
    E = res['popt'][0]
    E_err = np.sqrt(res['pcov'][0,0])
    ax.axhline(E, color='k', linestyle='--', label=rf'$E={E:.2f} \pm {E_err:.2f}$')
    ax.fill_between([0, max_x], 2*[E - E_err], 2*[E + E_err], color='0.8', linestyle='')

    if args.ED is not None:
        style_ED = dict(color='r', linestyle='--', linewidth=1.5, zorder=3)
        ax.axhline(args.ED, **style_ED, label=rf'$E={args.ED:.2f}$ [ED]')
        
    ax.set_xlim(0.0, max_x)
    ax.set_xlabel(r'$dt$')
    ax.legend()
    fig.set_tight_layout(True)
    plt.show()

if __name__ == '__main__':
    main()
