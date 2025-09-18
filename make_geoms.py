### Creates the geometries relevant for our runs on IBMQ.

import matplotlib.pyplot as plt
import numpy as np

def plot_geom_and_init(geom, init, *, title=None):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    geom_latt = np.where(geom != DUMMY, geom, float('nan'))
    ax.imshow(geom_latt, interpolation='none', cmap='Greys')
    cy, cx = np.meshgrid(np.arange(geom.shape[0]), np.arange(geom.shape[1]), indexing='ij')
    coords = np.stack((cx[init == 0], cy[init == 0]))
    ax.plot(*coords, marker='o', color='b', linestyle='', label='init 0')
    coords = np.stack((cx[init == 1], cy[init == 1]))
    ax.plot(*coords, marker='o', color='r', linestyle='', label='init 1')
    coords_fix = np.stack((cx[geom == 0], cy[geom == 0]))
    ax.plot(*coords_fix, marker='o', color='b', linestyle='', fillstyle='none', markersize=12, label='fixed 0')
    coords_fix = np.stack((cx[geom == 1], cy[geom == 1]))
    ax.plot(*coords_fix, marker='o', color='r', linestyle='', fillstyle='none', markersize=12, label='fixed 1')
    fig.legend(loc='center right')
    if title is not None:
        ax.set_title(title)
    return fig, ax

FREE = 0xff
DUMMY = 0xaa

### String representation to geom values or init values
### {-, |} -> fixed 0 (outside boundaries)
### {0, 1} -> fixed to that value (boundaries)
### {x, y} -> initialized to {0, 1} (free bulk)
### {' '} -> dummy (not on lattice)
def geom_char(c):
    if c == '-' or c == '|':
        return 0
    elif c == ' ':
        return DUMMY
    elif c == '0':
        return 0
    elif c == '1':
        return 1
    elif c == 'x' or c == 'y':
        return FREE
    else:
        raise ValueError(c)
def init_char(c):
    if c == '-' or c == '|' or c == ' ':
        return FREE
    elif c == 'x' or c == '0':
        return 0
    elif c == 'y' or c == '1':
        return 1
    else:
        raise ValueError(c)

def make_16T29P10O():
    # 16T29P10O is embedded in a 8x4 geometry
    NROWS = 8
    NCOLS = 4
    EROWS = 2*NROWS
    ECOLS = 4*NCOLS
    signature = [
        '----------------',
        '|   |   |   |   ',
        '------0---------',
        '  |   x   |   | ',
        '--0xxxxxxx0-----',
        '|   x   x   |   ',
        '0xxxxxxxxxxx0---',
        '  x   x   x   | ',
        '1yyyyyyyyyyy1---',
        '|   y   y   |   ',
        '--1yyyyyyy1-----',
        '  |   y   |   | ',
        '------1---------',
        '|   |   |   |   ',
        '----------------',
        '  |   |   |   | ',
    ]
    geom = np.stack(list(map(
        lambda line: np.array(list(map(geom_char, line)), dtype=np.uint8), signature
    )))
    init = np.stack(list(map(
        lambda line: np.array(list(map(init_char, line)), dtype=np.uint8), signature
    )))
    assert geom.shape == (EROWS, ECOLS)
    assert init.shape == geom.shape
    return geom, init

def make_42T72P18O():
    # 42T72P18O is also embedded in a 8x4 geometry
    NROWS = 8
    NCOLS = 5
    EROWS = 2*NROWS
    ECOLS = 4*NCOLS
    signature = [
        '----0---0---0-------',
        '|   x   x   x   |   ',
        '0xxxxxxxxxxxxxxx0---',
        '  x   x   x   x   | ',
        '0xxxxxxxxxxxxxxx0---',
        '    x   x   x       ',
        '0xxxxxxxxxxxxxxx0---',
        '  x   x   x   x   | ',
        '1yyyyyyyyyyyyyyy1---',
        '    y   y   y       ',
        '1yyyyyyyyyyyyyyy1---',
        '  y   y   y   y   | ',
        '1yyyyyyyyyyyyyyy1---',
        '|   y   y   y   |   ',
        '----1---1---1-------',
        '  |   |   |   |   | ',
    ]
    geom = np.stack(list(map(
        lambda line: np.array(list(map(geom_char, line)), dtype=np.uint8), signature
    )))
    init = np.stack(list(map(
        lambda line: np.array(list(map(init_char, line)), dtype=np.uint8), signature
    )))
    assert geom.shape == (EROWS, ECOLS)
    assert init.shape == geom.shape
    return geom, init
    
if __name__ == '__main__':
    geom, init = make_16T29P10O()
    geom.tofile('geoms/16T29P10O_geom.dat')
    init.tofile('geoms/16T29P10O_init.dat')
    fig, ax = plot_geom_and_init(geom, init, title='16T29P10O')
    fig.savefig('geoms/16T29P10O.pdf')
    geom, init = make_42T72P18O()
    geom.tofile('geoms/42T72P18O_geom.dat')
    init.tofile('geoms/42T72P18O_init.dat')
    fig, ax = plot_geom_and_init(geom, init, title='42T72P18O')
    fig.savefig('geoms/42T72P18O.pdf')
    plt.show()
