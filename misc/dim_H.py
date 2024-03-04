# We can count the dimension of the gauge-invariant subspace by considering the
# valid height-variable assignments within one hexagon per possible frozen
# boundary configuration, then "gluing" together these templates.
# One hexagon, with frozen boundaries, looks like:
#       y
# y-o-x-o-x-o-y
#   x       x
# y-o-x-o-x-o-y
#       y
# Where o = triangle, x = unfrozen petal, y = frozen petal.
# Valid height configurations are determined by:
#  * Applying unbounded integers to all intermediate variables
#  * Triangles = even, petals = odd
#  * Fixing the top external petal = +/- 1
#  * Fixing the other external petals to odd numbers [-9, 9]

import itertools
import numpy as np
import tqdm

'''
# variables in matrix order in the geometry drawn above
triangle_labels = [
    'o1', 'o2', 'o3',
    'o4', 'o5', 'o6'
]
petal_inner_labels = [
    'x1', 'x2',
    'x3', 'x4',
    'x5', 'x6'
]
petal_outer_labels = [
    'y1',
    'y2', 'y3',
    'y4', 'y5',
    'y6'
]


nodes = triangle_labels + petal_inner_labels + petal_outer_labels
triangles = [nodes.index(a) for a in triangle_labels]
petals_inner = [nodes.index(a) for a in petal_inner_labels]
petals_outer = [nodes.index(a) for a in petal_outer_labels]

edge_labels = [
    ('o1', 'y2'), ('o1', 'x1'), ('o1', 'x3'),
    ('o2', 'x1'), ('o2', 'y1'), ('o2', 'x2'),
    ('o3', 'x2'), ('o3', 'y3'), ('o3', 'x4'),
    ('o4', 'x3'), ('o4', 'y4'), ('o4', 'x5'),
    ('o5', 'x5'), ('o5', 'x6'), ('o5', 'y6'),
    ('o6', 'x6'), ('o6', 'x4'), ('o6', 'y5'),
]
edges = [
    (nodes.index(a), nodes.index(b)) for a,b in edge_labels
]
'''

def count_1hex():
    counts: dict[tuple[int,...], int] = {}
    o1 = 0
    for o2,o4 in itertools.product([o1-2,o1,o1+2], repeat=2):
        for o3 in [o2-2, o2, o2+2]:
            for o5 in [o4-2, o4, o4+2]:
                for o6 in [o5-2, o5, o5+2]:
                    if abs(o3-o6) > 2: continue
                    count = 2**(
                        (o1 == o2) + (o2 == o3) + (o3 == o6) +
                        (o6 == o5) + (o5 == o4) + (o4 == o1) )
                    for dy1,dy2,dy3,dy4,dy5,dy6 in itertools.product([-1,+1], repeat=6):
                        key = (
                            o1 + dy1, o2 + dy2, o3 + dy3,
                            o4 + dy4, o5 + dy5, o6 + dy6 )
                        if key not in counts:
                            counts[key] = 0
                        counts[key] += count
    return counts

def count_2hex(counts_1hex):
    # assume hexagons A and B are glued together such that y1A == y6B
    counts: dict[tuple[int,...], int] = {}
    for keyA, countA in tqdm.tqdm(counts_1hex.items()):
        petal_AB = keyA[5]
        outer_keyA = (keyA[0], keyA[1], keyA[2], keyA[3], keyA[4])
        # FORNOW: specialize for no external flux
        ext = outer_keyA[0]
        if not all_equal_to(outer_keyA, ext):
            continue
        for keyB, countB in counts_1hex.items():
            petal_BA = keyB[0]
            dAB = petal_AB - petal_BA
            outer_keyB = tuple(kB + dAB for kB in (keyB[1], keyB[2], keyB[3], keyB[4], keyB[5]))
            # FORNOW: specialize for no external flux
            if not all_equal_to(outer_keyB, ext):
                continue
            key = outer_keyA + outer_keyB
            if key not in counts:
                counts[key] = 0
            counts[key] += countA * countB
    return counts

def count_3hex_v1(counts_1hex):
    # assume hexagons A and B are glued together such that y1A == y6B, hexagons
    # B and C are glued together such that y5B == y2C, hexagons C and A are
    # glued together such that y3C == y4A
    counts: dict[tuple[int,...], int] = {}
    for keyA, countA in tqdm.tqdm(counts_1hex.items()):
        petal_AB = keyA[0]
        petal_AC = keyA[3]
        outer_keyA = (keyA[1], keyA[2], keyA[4], keyA[5])
        # FORNOW: specialize for no external flux
        ext = outer_keyA[0]
        if not all_equal_to(outer_keyA, ext):
            continue
        for keyB, countB in counts_1hex.items():
            petal_BA = keyB[5]
            petal_BC = keyB[4]
            dAB = petal_AB - petal_BA
            outer_keyB = tuple(kB + dAB for kB in (keyB[0], keyB[1], keyB[2], keyB[3]))
            # FORNOW: specialize for no external flux
            if not all_equal_to(outer_keyB, ext):
                continue
            for keyC, countC in counts_1hex.items():
                petal_CB = keyC[1]
                petal_CA = keyC[2]
                dBC = petal_BC - petal_CB
                dCA = petal_CA - petal_AC
                # Interior Gauss Law
                if dAB + dBC + dCA != 0: continue
                outer_keyC = tuple(kC - dCA for kC in (keyC[0], keyC[3], keyC[4], keyC[5]))
                # FORNOW: specialize for no external flux
                if not all_equal_to(outer_keyC, ext):
                    continue
                key = (
                    outer_keyA + 
                    outer_keyB +
                    outer_keyC
                )
                if key not in counts:
                    counts[key] = 0
                counts[key] += countA * countB * countC
    return counts

def count_3hex_v2(counts_1hex):
    # assume hexagons A and B are glued together such that y1A == y6B, hexagons
    # B and C are glued together such that y4B == y3C
    counts: dict[tuple[int,...], int] = {}
    for keyA, countA in tqdm.tqdm(counts_1hex.items()):
        petal_AB = keyA[0]
        outer_keyA = (keyA[1], keyA[2], keyA[3], keyA[4], keyA[5])
        # FORNOW: specialize for no external flux
        ext = outer_keyA[0]
        if not all_equal_to(outer_keyA, ext):
            continue
        for keyB, countB in counts_1hex.items():
            petal_BA = keyB[5]
            petal_BC = keyB[3]
            dAB = petal_AB - petal_BA
            outer_keyB = tuple(kB + dAB for kB in (keyB[0], keyB[1], keyB[2], keyB[4]))
            # FORNOW: specialize for no external flux
            if not all_equal_to(outer_keyB, ext):
                continue
            for keyC, countC in counts_1hex.items():
                petal_CB = keyC[2]
                dBC = petal_BC - petal_CB
                outer_keyC = tuple(kC + dAB + dBC for kC in (keyC[0], keyC[1], keyC[3], keyC[4], keyC[5]))
                # FORNOW: specialize for no external flux
                if not all_equal_to(outer_keyC, ext):
                    continue
                key = (
                    outer_keyA + 
                    outer_keyB +
                    outer_keyC
                )
                if key not in counts:
                    counts[key] = 0
                counts[key] += countA * countB * countC
    return counts

def all_equal_to(tup, value):
    return tup == (value,) * len(tup)
def all_equal(tup):
    return all_equal_to(tup, tup[0])

def count_open(counts):
    tot = 0
    for key, value in counts.items():
        if not all_equal(key): continue
        tot += value
    return tot

def to_numpy(counts):
    tab = []
    for k,v in counts.items():
        tab.append(k + (v,))
    return np.stack(tab)

if __name__ == '__main__':
    counts_1hex = count_1hex()
    np.save('counts_1hex.npy', to_numpy(counts_1hex))
    print(f'== 1 Hexagon ==')
    print(f'num keys = {len(counts_1hex)}')
    print(f'total = {sum(counts_1hex.values())}')
    print(f'open = {count_open(counts_1hex)}')

    counts_2hex = count_2hex(counts_1hex)
    print(f'== 2 Hexagons ==')
    print(f'num keys = {len(counts_2hex)}')
    print(f'total = {sum(counts_2hex.values())}')
    print(f'open = {count_open(counts_2hex)}')
    
    counts_3hex = count_3hex_v1(counts_1hex)
    print(f'== 3 Hexagons (triangle) ==')
    print(f'num keys = {len(counts_3hex)}')
    print(f'total = {sum(counts_3hex.values())}')
    print(f'open = {count_open(counts_3hex)}')

    counts_3hex = count_3hex_v2(counts_1hex)
    print(f'== 3 Hexagons (line) ==')
    print(f'num keys = {len(counts_3hex)}')
    print(f'total = {sum(counts_3hex.values())}')
    print(f'open = {count_open(counts_3hex)}')
