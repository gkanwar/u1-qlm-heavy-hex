# U(1) QLM on the Heavy Hex lattice
Quantum Monte Carlo code for the U(1) QLM on the Heavy Hex lattice. This repo
implements a cluster algorithm on the dual height variables.

The simulation is implemented under `src/main.c` and can be built from the
top level by running
```
make NCOLS=<NCOLS> NROWS=<NROWS> NT=<NT> (release|debug)
```
This will produce a binary at
```
bin/cluster.<NT>_<NROWS>_<NCOLS>.(release|debug)
```

# Geometry
Everything is set up in terms of a mapping of the hexagonal lattice into a
larger rectangular grid. This is a little wasteful in that many elements of
the embedding array are not dynamical sites, but makes calculations a bit
simpler.

For example, for `NROWS = 2` and `NCOLS = 2`, we have the following periodic
geometry:
```
o-x-o-x-o-x-o-x-
x       x
o-x-o-x-o-x-o-x-
    x       x
```
The embedding array has geometry `EROWS x ECOLS`, where `EROWS = 2*NROWS` and
`ECOLS = 4*NCOLS` and only 5/8 of the array is populated. The triangular
plaquettes are given by the 'o's at the (even,even) sites of the embedding
array. The petal plaquettes are given by the 'x's at the (even,odd),
(1 mod 4, 0 mod 4), and (3 mod 4, 2 mod 4) sites of the embedding array.

Note that `NCOLS` can be odd, but `NROWS` should be even to avoid weird gluing
at the bottom/top boundary.

In the temporal direction, we simply stack this geometry `NT` times. In layer
t, the petals are considered to be at slice t-1/2 while the triangles are
considered to be at slice t.

# Dirac strings
Dirac strings are necessary to introduce odd static charges. We restrict to only
allow Dirac strings on petals only:

  * For odd petals they modify the interaction with the triangle above.
  * For even petals they modify the interaction with the left triangle.
 
In either case, the petal and triangle see each other as the opposite
of their value.

# Command-line arguments
The cluster simulation binary takes several command line arguments:

  * `-t <dt>`: temporal lattice spacing
  * `-p <KP>`: petal interaction coupling
  * `-e <KE>`: electric term coupling
  * `-f <prefix>`: prefix for output files
  * `-s <save_freq>`: steps between saved configurations (set this large, since configs take space)
  * `-m <meas_freq>`: steps between measurements (usually smaller than `save_freq`)
  * `-i <n_iter>`: number of total MCMC iterations
  * `-r <seed>`: seed for RNG
  * `-b (pbc|obc|...|file)`: set the geometry / boundary conditions
  * `-c (cold|cold_str|hot|file)`: set the initialization strategy
  * `-y <fname>`: optional file to set geometry (if `-b file`)
  * `-z <fname>`: optional file to initialize spins (if `-c file`)
  * `-d <fname>`: optional file to set dirac string

See `make_geoms.py` and `make_dirac_strs.py` for some hints on the structure of
the data files to be passed above. They are essentially just packed binary
arrays with appropriate values assigned to the positions of the `EROWS x ECOLS`
embedding array.
