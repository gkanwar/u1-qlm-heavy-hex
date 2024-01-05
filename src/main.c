#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef NROWS
#define NROWS 4
#endif
#ifndef NCOLS
#define NCOLS 4
#endif
#ifndef NT
#define NT 8
#endif

/// We can wastefully map the heavy hex into a rectangular grid. For example,
/// for NROWS = 2 and NCOLS = 2, we have the following periodic geometry:
/// o-x-o-x-o-x-o-x-
/// x       x       
/// o-x-o-x-o-x-o-x-
///     x       x
///
/// The embedding array has geometry (2*NROWS) x (4*NCOLS) and only 5/8 of the
/// array is populated. The triangular plaquettes are given by the 'o's at the
/// (even,even) sites of the embedding array. The petal plaquettes are given by
/// the 'x's at the (even,odd), (1 mod 4, 0 mod 4), and (3 mod 4, 2 mod 4) sites
/// of the embedding array.
///
/// Note that NCOLS can be odd, but NROWS should be even to avoid weird gluing
/// at the bottom/top boundary.
///
/// In the temporal direction, we simply stack this geometry NT times. In layer
/// t, the petals are considered to be at slice t-1/2 while the triangles are
/// considered to be at slice t.

#define EROWS (2*NROWS)
#define ECOLS (4*NCOLS)
static bool lattice[NT * EROWS * ECOLS];

static inline int ind_tri(int t, int i, int j) {
  assert(t < NT && i < NROWS && j < 2*NCOLS);
  return (t*EROWS + 2*i)*ECOLS + 2*j;
}
static inline int ind_pet_even(int t, int i, int j) {
  assert(t < NT && i < NROWS && j < 2*NCOLS);
  return (t*EROWS + 2*i)*ECOLS + 2*j + 1;
}
static inline int ind_pet_odd(int t, int i, int j) {
  assert(t < NT && i < NROWS && j < NCOLS);
  return (t*EROWS + 2*i+1)*ECOLS + 4*j + 2*(i % 2);
}

static inline int wrap_i(int i) {
  return (NROWS+i) % NROWS;
}
static inline int wrap_sj(int j) {
  return (NCOLS+j) % NCOLS;
}
static inline int wrap_dj(int j) {
  return (2*NCOLS+j) % (2*NCOLS);
}
static inline int wrap_t(int t) {
  return (NT+t) % NT;
}

// bit mask 00 | tfwd | tbwd | left | right | up | down
//  - up/down bonds are irrelevant for petals
//  - left/right bonds are three-way for petals
typedef enum {
  DOWN,
  UP,
  RIGHT,
  LEFT,
  T_BWD,
  T_FWD
} bond_bit;
typedef uint8_t bond_t;
static inline bool bond_d(bond_t x) {
  return (x >> DOWN) & 1;
}
static inline bool bond_u(bond_t x) {
  return (x >> UP) & 1;
}
static inline bool bond_r(bond_t x) {
  return (x >> RIGHT) & 1;
}
static inline bool bond_l(bond_t x) {
  return (x >> LEFT) & 1;
}
static inline bool bond_tbwd(bond_t x) {
  return (x >> T_BWD) & 1;
}
static inline bool bond_tfwd(bond_t x) {
  return (x >> T_FWD) & 1;
}

static struct {
  bond_t tri[NT][NROWS][2*NCOLS];
  bond_t pet_even[NT][NROWS][2*NCOLS];
  bond_t pet_odd[NT][NROWS][NCOLS];
} bonds;

static struct {
  bool tri[NT][NROWS][2*NCOLS];
  bool pet_even[NT][NROWS][2*NCOLS];
  bool pet_odd[NT][NROWS][NCOLS];
} seen;


// statically alloc'd global queue
// capacity = maximum number of sites in a sublattice
#define QUEUE_CAP NROWS * 3*NCOLS
typedef enum {
  ODD, EVEN
} parity;
typedef struct {
  int t, i, j;
  parity p;
} coord;
static struct {
  coord values[QUEUE_CAP];
  int size;
} queue;
static inline int empty() {
  return queue.size == 0;
}
static inline coord pop() {
  assert(!empty());
  return queue.values[--queue.size];
}
static inline void push(coord idx) {
  assert(queue.size < QUEUE_CAP);
  queue.values[queue.size++] = idx;
}


static inline bool rand_bool() {
  assert(RAND_MAX % 2 == 1);
  return rand() % 2;
}
static inline double rand_double() {
  // TODO: may have discretization artifacts if RAND_MAX too small
  return rand() / (double)RAND_MAX;
}


// couplings betaX = t J_x
static double betaT = 1.0;
static double betaP = 1.0;
static double betaPP = 0.0;

void init_cold() {
  memset(lattice, 0, sizeof(lattice));
}
bool tri_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  double w1 = exp(-betaPP)*sinh(betaP) + delta_b*(exp(-betaPP-betaP) - exp(betaPP));
  double w0 = delta_b*exp(betaPP);
  assert(w1 >= 0 && w0 >= 0);
  double p = w1 / (w1 + w0);
  return rand_double() < p;
}

bool tri_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  if (!delta_b) {
    return true;
  }
  double w1 = exp(-betaT);
  double w0 = sinh(betaT);
  assert(w1 >= 0 && w0 >= 0);
  double p = w1 / (w1 + w0);
  return rand_double() < p;
}

bool pet_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  double w1 = sinh(betaT) + delta_a*(exp(-betaT) - 1.0);
  double w0 = (double)delta_a;
  assert(w1 >= 0 && w0 >= 0);
  double p = w1 / (w1 + w0);
  return rand_double() < p;
}

bool pet_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  if (!delta_a) {
    return true;
  }
  double w1 = exp(-betaPP-betaP);
  double w0 = exp(-betaPP)*sinh(betaP);
  assert(w1 >= 0 && w0 >= 0);
  double p = w1 / (w1 + w0);
  return rand_double() < p;
}

void sample_tri_bonds() {
  memset(bonds.tri, 0, sizeof(bonds.tri));
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        // relevant triangles
        bool tri = lattice[ind_tri(t, i, j)];
        int t_fwd = wrap_t(t+1);
        bool tri_fwd = lattice[ind_tri(t_fwd, i, j)];
        int j_l = wrap_dj(j-1);
        int j_r = wrap_dj(j+1);
        bool tri_l = lattice[ind_tri(t, i, j_l)];
        bool tri_r = lattice[ind_tri(t, i, j_r)];
        int tri_i_ud = wrap_i(((i+j) % 2 == 0) ? i+1 : i-1);
        int pet_i_ud = wrap_i(((i+j) % 2 == 0) ? i : i-1);
        bond_bit UD = ((i+j) % 2 == 0) ? DOWN : UP;
        bond_bit FLIP_UD = ((i+j) % 2 == 0) ? UP : DOWN;
        bool tri_ud = lattice[ind_tri(t, tri_i_ud, j)];
        
        // relevant petals
        bool pet_bwd_l = lattice[ind_pet_even(t, i, j_l)];
        bool pet_bwd_r = lattice[ind_pet_even(t, i, j)];
        bool pet_bwd_ud = lattice[ind_pet_odd(t, pet_i_ud, j)];
        bool pet_fwd_l = lattice[ind_pet_even(t_fwd, i, j_l)];
        bool pet_fwd_r = lattice[ind_pet_even(t_fwd, i, j)];
        bool pet_fwd_ud = lattice[ind_pet_odd(t_fwd, pet_i_ud, j)];

        // bonds
        // T_FWD / T_BWD
        {
          bool delta_a = (tri == tri_fwd);
          bool delta_b = (pet_fwd_l == pet_fwd_r && pet_fwd_r == pet_fwd_ud);
          if (tri_bond_temporal(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << T_FWD;
            bonds.tri[t_fwd][i][j] |= 1 << T_BWD;
          }
        }
        if (i+j % 2 != 0) {
          continue;
        }
        // UP / DOWN
        {
          bool delta_a = (tri == tri_ud);
          bool delta_b = (pet_fwd_ud == pet_bwd_ud);
          if (tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << UD;
            bonds.tri[t][tri_i_ud][j] |= 1 << FLIP_UD;
          }
        }
        // RIGHT
        {
          bool delta_a = (tri == tri_r);
          bool delta_b = (pet_fwd_r == pet_bwd_r);
          if (tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << RIGHT;
            bonds.tri[t][i][j_r] |= 1 << LEFT;
          }
        }
        // LEFT
        {
          bool delta_a = (tri == tri_l);
          bool delta_b = (pet_fwd_l == pet_bwd_l);
          if (tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << LEFT;
            bonds.tri[t][i][j_l] |= 1 << RIGHT;
          }
        }
      }
    }
  }
}

void flood_fill_tri(bool spin) {
  while (!empty()) {
    coord x = pop();
    int t = x.t, i = x.i, j = x.j;
    assert(!seen.tri[t][i][j]);
    seen.tri[t][i][j] = true;
    int idx = ind_tri(t, i, j);
    lattice[idx] = spin;
    bond_t bond = bonds.tri[t][i][j];
    // TFWD
    if (bond_tfwd(bond)) {
      int tfwd = wrap_t(t+1);
      if (!seen.tri[tfwd][i][j]) {
        push((coord){tfwd, i, j});
      }
    }
    // TBWD
    if (bond_tbwd(bond)) {
      int tbwd = wrap_t(t-1);
      if (!seen.tri[tbwd][i][j]) {
        push((coord){tbwd, i, j});
      }
    }
    // LEFT
    if (bond_l(bond)) {
      int j_l = wrap_dj(j-1);
      if (!seen.tri[t][i][j_l]) {
        push((coord){t, i, j_l});
      }
    }
    // RIGHT
    if (bond_r(bond)) {
      int j_r = wrap_dj(j-1);
      if (!seen.tri[t][i][j_r]) {
        push((coord){t, i, j_r});
      }
    }
    // UP/DOWN
    if (bond_u(bond) || bond_d(bond)) {
      int i_ud = wrap_i(((i+j) % 2 == 0) ? i+1 : i-1);
      if (!seen.tri[t][i_ud][j]) {
        push((coord){t, i_ud, j});
      }
    }
  }
}

void update_tri_spins() {
  memset(seen.tri, 0, sizeof(seen.tri));
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.tri[t][i][j]) {
          continue;
        }
        push((coord){t, i, j});
        flood_fill_tri(rand_bool());
      }
    }
  }
}

void sample_pet_bonds() {
  memset(bonds.pet_even, 0, sizeof(bonds.pet_even));
  memset(bonds.pet_odd, 0, sizeof(bonds.pet_odd));
  // even petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        // relevant petals
        bool pet = lattice[ind_pet_even(t, i, j)];
        int t_fwd = wrap_t(t+1);
        bool pet_fwd = lattice[ind_pet_even(t_fwd, i, j)];
        int j_l = j;
        int j_r = wrap_dj(j+1);

        // relevant triangles
        bool tri_l = lattice[ind_tri(t, i, j_l)];
        bool tri_r = lattice[ind_tri(t, i, j_r)];

        // bonds
        // T_FWD / T_BWD
        {
          bool delta_b = (pet == pet_fwd);
          bool delta_a = (tri_l == tri_r);
          if (pet_bond_temporal(delta_a, delta_b)) {
            bonds.pet_even[t][i][j] |= 1 << T_FWD;
            bonds.pet_even[t_fwd][i][j] |= 1 << T_BWD;
          }
        }
      }
    }
  }
  // odd petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        // relevant petals
        bool pet = lattice[ind_pet_odd(t, i, j)];
        int t_fwd = wrap_t(t+1);
        int t_bwd = wrap_t(t-1);
        bool pet_fwd = lattice[ind_pet_odd(t_fwd, i, j)];
        int i_u = i;
        int i_d = wrap_i(i+1);
        int j_even = 2*j + 2*(i % 2);
        int j_even_r = j_even;
        int j_even_l = wrap_dj(j_even - 1);
        bool pet_ur = lattice[ind_pet_even(t, i_u, j_even_r)];
        bool pet_ul = lattice[ind_pet_even(t, i_u, j_even_l)];
        bool pet_dr = lattice[ind_pet_even(t, i_d, j_even_r)];
        bool pet_dl = lattice[ind_pet_even(t, i_d, j_even_l)];

        // relevant triangles
        bool tri_u = lattice[ind_tri(t, i_u, j_even)];
        bool tri_d = lattice[ind_tri(t, i_d, j_even)];
        bool tri_bwd_u = lattice[ind_tri(t_bwd, i_u, j_even)];
        bool tri_bwd_d = lattice[ind_tri(t_bwd, i_d, j_even)];

        // bonds
        // T_FWD / T_BWD
        {
          bool delta_b = (pet == pet_fwd);
          bool delta_a = (tri_u == tri_d);
          if (pet_bond_temporal(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << T_FWD;
            bonds.pet_odd[t_fwd][i][j] |= 1 << T_BWD;
          }
        }
        // UP
        {
          bool delta_b = (pet == pet_ur && pet_ur == pet_ul);
          bool delta_a = (tri_u == tri_bwd_u);
          if (pet_bond_spatial(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << UP;
            bonds.pet_even[t][i_u][j_even_r] |= 1 << LEFT;
            bonds.pet_even[t][i_u][j_even_l] |= 1 << RIGHT;
          }
        }
        // DOWN
        {
          bool delta_b = (pet == pet_dr && pet_dr == pet_dl);
          bool delta_a = (tri_d == tri_bwd_d);
          if (pet_bond_spatial(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << DOWN;
            bonds.pet_even[t][i_d][j_even_r] |= 1 << LEFT;
            bonds.pet_even[t][i_d][j_even_l] |= 1 << RIGHT;
          }
        }
      }
    }
  }
}

void flood_fill_pet(bool spin) {
  while (!empty()) {
    coord x = pop();
    int t = x.t, i = x.i, j = x.j;
    // even petals
    if (x.p == EVEN) {
      assert(!seen.pet_even[t][i][j]);
      seen.pet_even[t][i][j] = true;
      int idx = ind_pet_even(t, i, j);
      lattice[idx] = spin;
      int bond = bonds.pet_even[t][i][j];
      // TFWD
      if (bond_tfwd(bond)) {
        int tfwd = wrap_t(t+1);
        if (!seen.pet_even[tfwd][i][j]) {
          push((coord){tfwd, i, j, EVEN});
        }
      }
      // TBWD
      if (bond_tbwd(bond)) {
        int tbwd = wrap_t(t-1);
        if (!seen.pet_even[tbwd][i][j]) {
          push((coord){tbwd, i, j, EVEN});
        }
      }
      // LEFT
      if (bond_l(bond)) {
        int j_le = wrap_dj(j-1);
        int j_lo = j/2;
        int i_l = ((i+j) % 2 == 0) ? i : wrap_i(i-1);
        if (!seen.pet_even[t][i][j_le]) {
          push((coord){t, i, j_le, EVEN});
        }
        if (!seen.pet_odd[t][i_l][j_lo]) {
          push((coord){t, i_l, j_lo, ODD});
        }
      }
      // RIGHT
      if (bond_r(bond)) {
        int j_re = wrap_dj(j+1);
        int j_ro = wrap_dj(j+1)/2;
        int i_r = ((i+j) % 2 == 0) ? wrap_i(i-1) : i;
        if (!seen.pet_even[t][i][j_re]) {
          push((coord){t, i, j_re, EVEN});
        }
        if (!seen.pet_odd[t][i_r][j_ro]) {
          push((coord){t, i_r, j_ro, ODD});
        }
      }
      assert(!bond_u(bond));
      assert(!bond_d(bond));
    }
    // odd petals
    else {
      assert(!seen.pet_odd[t][i][j]);
      seen.pet_odd[t][i][j] = true;
      int idx = ind_pet_odd(t, i, j);
      lattice[idx] = spin;
      int bond = bonds.pet_odd[t][i][j];
      // TFWD
      if (bond_tfwd(bond)) {
        int tfwd = wrap_t(t+1);
        if (!seen.pet_odd[tfwd][i][j]) {
          push((coord){tfwd, i, j, ODD});
        }
      }
      // TBWD
      if (bond_tbwd(bond)) {
        int tbwd = wrap_t(t-1);
        if (!seen.pet_odd[tbwd][i][j]) {
          push((coord){tbwd, i, j, ODD});
        }
      }
      // UP
      if (bond_u(bond)) {
        int j_l = 2*j;
        int j_r = wrap_dj(2*j+1);
        if (!seen.pet_even[t][i][j_l]) {
          push((coord){t, i, j_l, EVEN});
        }
        if (!seen.pet_even[t][i][j_r]) {
          push((coord){t, i, j_r, EVEN});
        }
      }
      // DOWN
      if (bond_d(bond)) {
        int j_l = 2*j;
        int j_r = wrap_dj(2*j+1);
        int i_d = wrap_i(i+1);
        if (!seen.pet_even[t][i_d][j_l]) {
          push((coord){t, i_d, j_l, EVEN});
        }
        if (!seen.pet_even[t][i_d][j_r]) {
          push((coord){t, i_d, j_r, EVEN});
        }
      }
      assert(!bond_l(bond));
      assert(!bond_r(bond));
    }
  }
}

void update_pet_spins() {
  memset(seen.pet_even, 0, sizeof(seen.pet_even));
  memset(seen.pet_odd, 0, sizeof(seen.pet_odd));
  // outer loop even
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.pet_even[t][i][j]) {
          continue;
        }
        push((coord){t, i, j, EVEN});
        flood_fill_pet(rand_bool());
      }
    }
  }
  // outer loop odd
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (seen.pet_odd[t][i][j]) {
          continue;
        }
        push((coord){t, i, j, ODD});
        flood_fill_pet(rand_bool());
      }
    }
  }
}


int main() {
  int n_iter = 100;
  srand(1234);
  init_cold();
  for (int i = 0; i < n_iter; ++i) {
    // triangle sublattice
    sample_tri_bonds();
    update_tri_spins();
    // petal sublattice
    sample_pet_bonds();
    update_pet_spins();
  }
}
