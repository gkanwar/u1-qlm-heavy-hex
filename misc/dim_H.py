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
from typing import Optional

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
    keyA: tuple[int,...]
    for keyA, countA in tqdm.tqdm(counts_1hex.items()):
        petal_AB = keyA[5]
        outer_keyA = (keyA[0], keyA[1], keyA[2], keyA[3], keyA[4])
        # FORNOW: specialize for no external flux
        ext = outer_keyA[0]
        if not all_equal_to(outer_keyA, ext):
            continue
        keyB: tuple[int,...]
        for keyB, countB in counts_1hex.items():
            petal_BA = keyB[0]
            dAB = petal_AB - petal_BA
            outer_keyB = tuple(kB + dAB for kB in (keyB[1], keyB[2], keyB[3], keyB[4], keyB[5]))
            # FORNOW: specialize for no external flux
            if not all_equal_to(outer_keyB, ext):
                continue
            key: tuple[int,...]
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

def valid_adjacency(adj):
    for i,neighbors in adj.items():
        for j in neighbors:
            if i not in adj[j]:
                return False
    return True

# TODO: hex_adjacency and fixed_values can be packed into one structure, because
# they take the value None in complementary cases.
def count_generic(
        counts_1hex, hexes: list[str], hex_adjacency: dict[str, list[Optional[str]]],
        fixed_values: dict[str, list[Optional[int]]]):
    assert valid_adjacency(hex_adjacency)
    tot = 0
    first: str = hexes[0]
    adj_first: list[Optional[str]]
    fixed_first: list[Optional[int]]
    adj_first = hex_adjacency[first]
    fixed_first = fixed_values[first]

    # base case
    if len(hexes) == 1:
        assert not any(x is None for x in fixed_first)
        key = tuple(fixed_first)
        return counts_1hex.get(key, 0)

    # recursive case
    it = tqdm.tqdm(counts_1hex.items(), position=len(hexes)-1, leave=False, desc=first)
    for key, count in it:
        offset = None
        assert len(fixed_first) == len(key)
        skip = False
        for value, fixed in zip(key, fixed_first):
            if fixed is not None:
                if offset is None:
                    offset = fixed - value
                if value + offset != fixed:
                    skip = True
                    break
        if skip: continue
        assert offset is not None, 'must pass hexes ordered from boundary to interior'
        # NOTE: this hex always fixes the same values, so no need to undo this
        for i,neighbor in enumerate(adj_first):
            if neighbor not in hexes or neighbor is None:
                assert fixed_first[i] is not None
                continue
            ind = hex_adjacency[neighbor].index(first)
            fixed_values[neighbor][ind] = key[i] + offset
        tot += count * count_generic(counts_1hex, hexes[1:], hex_adjacency, fixed_values)
    return tot

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

    generic_2hex = 0
    for env in [-1, 1]:
        hexes = ['A', 'B']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, None, None, None],
            'B': [None, None, None, None, None, 'A'],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, env, env, env],
            'B': [env, env, env, env, env, None],
        }
        generic_2hex += count_generic(counts_1hex, hexes, hex_adj, fixed_values)
    print(f'open (generic) = {generic_2hex}')
    
    counts_3hex = count_3hex_v1(counts_1hex)
    print(f'== 3 Hexagons (triangle) ==')
    print(f'num keys = {len(counts_3hex)}')
    print(f'total = {sum(counts_3hex.values())}')
    print(f'open = {count_open(counts_3hex)}')

    generic_3hex = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, 'C', None, None],
            'B': [None, None, None, None, 'C', 'A'],
            'C': [None, 'B', 'A', None, None, None],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, None, env, env],
            'B': [env, env, env, env, None, None],
            'C': [env, None, None, env, env, env]
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic_3hex += count_env
    print(f'open (generic) = {generic_3hex}')

    counts_3hex = count_3hex_v2(counts_1hex)
    print(f'== 3 Hexagons (line) ==')
    print(f'num keys = {len(counts_3hex)}')
    print(f'total = {sum(counts_3hex.values())}')
    print(f'open = {count_open(counts_3hex)}')

    generic_3hex = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, None, None, None],
            'B': [None, None, None, 'C', None, 'A'],
            'C': [None, None, 'B', None, None, None],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, env, env, env],
            'B': [env, env, env, None, env, None],
            'C': [env, env, None, env, env, env]
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic_3hex += count_env
    print(f'open (generic) = {generic_3hex}')

    print(f'== 4 Hexagons (rhombus) ==')
    generic = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C', 'D']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, 'C', None, None],
            'B': [None, None, None, 'D', 'C', 'A'],
            'C': ['D', 'B', 'A', None, None, None],
            'D': [None, None, 'B', None, None, 'C'],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, None, env, env],
            'B': [env, env, env, None, None, None],
            'C': [None, None, None, env, env, env],
            'D': [env, env, None, env, env, None],
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic += count_env
    print(f'open (generic) = {generic}')

    print(f'== 5 Hexagons (trapezoid) ==')
    generic = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C', 'E', 'D']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, 'C', None, None],
            'B': [None, None, None, 'D', 'C', 'A'],
            'C': ['D', 'B', 'A', 'E', None, None],
            'D': [None, None, 'B', None, 'E', 'C'],
            'E': [None, 'D', 'C', None, None, None],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, None, env, env],
            'B': [env, env, env, None, None, None],
            'C': [None, None, None, None, env, env],
            'D': [env, env, None, env, None, None],
            'E': [env, None, None, env, env, env],
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic += count_env
    print(f'open (generic) = {generic}')

    print(f'== 6 Hexagons (crescent) ==')
    generic = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C', 'E', 'F', 'D']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, 'C', None, None],
            'B': [None, None, None, 'D', 'C', 'A'],
            'C': ['D', 'B', 'A', 'E', 'F', None],
            'D': [None, None, 'B', None, 'E', 'C'],
            'E': [None, 'D', 'C', None, None, 'F'],
            'F': ['E', 'C', None, None, None, None],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, None, env, env],
            'B': [env, env, env, None, None, None],
            'C': [None, None, None, None, None, env],
            'D': [env, env, None, env, None, None],
            'E': [env, None, None, env, env, None],
            'F': [None, None, env, env, env, env],
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic += count_env
    print(f'open (generic) = {generic}')

    print(f'== 7 Hexagons (hex) ==')
    generic = 0
    for env in [-1, 1]:
        hexes = ['A', 'B', 'C', 'E', 'F', 'G', 'D']
        hex_adj: dict[str, list[Optional[str]]] = {
            'A': ['B', None, None, 'C', 'G', None],
            'B': [None, None, None, 'D', 'C', 'A'],
            'C': ['D', 'B', 'A', 'E', 'F', 'G'],
            'D': [None, None, 'B', None, 'E', 'C'],
            'E': [None, 'D', 'C', None, None, 'F'],
            'F': ['E', 'C', 'G', None, None, None],
            'G': ['C', 'A', None, 'F', None, None],
        }
        fixed_values: dict[str, list[Optional[int]]] = {
            'A': [None, env, env, None, None, env],
            'B': [env, env, env, None, None, None],
            'C': [None, None, None, None, None, None],
            'D': [env, env, None, env, None, None],
            'E': [env, None, None, env, env, None],
            'F': [None, None, None, env, env, env],
            'G': [None, None, env, None, env, env],
        }
        count_env = count_generic(counts_1hex, hexes, hex_adj, fixed_values)
        print(f'... {env} -> {count_env}')
        generic += count_env
    print(f'open (generic) = {generic}')
