#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef NROWS
#define NROWS 2
#endif
#ifndef NCOLS
#define NCOLS 2
#endif
#ifndef NT
#define NT 64
#endif

static const double dt = 0.1;
static const double KP = 1.0;
static const double KE = 1.0;
// derived couplings betaX = dt K_x
static const double betaT = dt;
static const double betaP = dt * KP;
static const double betaE = dt * KE;


// settings for debugging
#define ENABLE_PET_SPATIAL 1
#define ENABLE_TRI_SPATIAL 1


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
typedef uint8_t spin_t;
static spin_t lattice[NT * EROWS * ECOLS];

#define FREE 0xff
static spin_t fixed[EROWS * ECOLS];

static void init_fixed_obc() {
  for (int x = 0; x < EROWS; ++x) {
    for (int y = 0; y < ECOLS; ++y) {
      int i = x*ECOLS + y;
      if (x == EROWS-1 || y >= ECOLS-3) {
        fixed[i] = false;
      }
      else {
        fixed[i] = FREE;
      }
    }
  }
}

static void init_fixed_pbc() {
  for (int x = 0; x < EROWS; ++x) {
    for (int y = 0; y < ECOLS; ++y) {
      int i = x*ECOLS + y;
      fixed[i] = FREE;
    }
  }
}

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
// static inline int wrap_sj(int j) {
//   return (NCOLS+j) % NCOLS;
// }
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
typedef uint8_t bond_dirs_t;
static inline bool bond_d(bond_dirs_t x) {
  return (x >> DOWN) & 1;
}
static inline bool bond_u(bond_dirs_t x) {
  return (x >> UP) & 1;
}
static inline bool bond_r(bond_dirs_t x) {
  return (x >> RIGHT) & 1;
}
static inline bool bond_l(bond_dirs_t x) {
  return (x >> LEFT) & 1;
}
static inline bool bond_tbwd(bond_dirs_t x) {
  return (x >> T_BWD) & 1;
}
static inline bool bond_tfwd(bond_dirs_t x) {
  return (x >> T_FWD) & 1;
}

static struct {
  bond_dirs_t tri[NT][NROWS][2*NCOLS];
  bond_dirs_t pet_even[NT][NROWS][2*NCOLS];
  bond_dirs_t pet_odd[NT][NROWS][NCOLS];
} bonds;

static struct {
  bool tri[NT][NROWS][2*NCOLS];
  bool pet_even[NT][NROWS][2*NCOLS];
  bool pet_odd[NT][NROWS][NCOLS];
} seen;


// statically alloc'd global queue
// capacity = maximum number of sites in a sublattice
#define QUEUE_CAP (NT * NROWS * 3*NCOLS)
typedef enum {
  ODD, EVEN
} pet_parity;
typedef struct {
  int t, i, j;
  pet_parity p;
} coord_t;
static struct {
  coord_t values[QUEUE_CAP];
  int size;
} queue;
static inline void q_init() {
  queue.size = 0;
}
static inline int q_empty() {
  return queue.size == 0;
}
static inline coord_t q_pop() {
  assert(!q_empty());
  return queue.values[--queue.size];
}
static inline void q_push(coord_t idx) {
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


void init_lat_cold() {
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        lattice[ind_tri(t, i, j)] = false;
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        lattice[ind_pet_even(t, i, j)] = false;
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        lattice[ind_pet_odd(t, i, j)] = false;
      }
    }
  }
}

void init_lat_hot_tri() {
  init_lat_cold();
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        lattice[ind_tri(t, i, j)] = rand_bool();
      }
    }
  }
}

void init_lat_hot_petal() {
  init_lat_cold();
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        lattice[ind_pet_even(t, i, j)] = rand_bool();
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        lattice[ind_pet_odd(t, i, j)] = rand_bool();
      }
    }
  }
}

void init_lat_apply_fixed() {
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        spin_t value = fixed[ind_tri(0, i, j)];
        if (value != FREE) {
          lattice[ind_tri(t, i, j)] = value;
        }
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        spin_t value = fixed[ind_pet_even(0, i, j)];
        if (value != FREE) {
          lattice[ind_pet_even(t, i, j)] = value;
        }
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        spin_t value = fixed[ind_pet_odd(0, i, j)];
        if (value != FREE) {
          lattice[ind_pet_odd(t, i, j)] = value;
        }
      }
    }
  }
}

bool tri_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  if (!delta_b) {
    return true;
  }
  double p = 1.0 - exp(-2*betaE) / cosh(betaP);
  assert(0.0 <= p && p <= 1.0);
  return rand_double() < p;
}

bool tri_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  if (!delta_b) {
    return true;
  }
  double p = 1.0 - tanh(betaT);
  assert(0.0 <= p && p <= 1.0);
  return rand_double() < p;
}

bool pet_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  if (!delta_a) {
    return true;
  }
  double p = 1.0 - 1.0 / cosh(betaT);
  assert(0.0 <= p && p <= 1.0);
  return rand_double() < p;
}

bool pet_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  if (!delta_a) {
    return true;
  }
  double p = 1.0 - tanh(betaP);
  assert(0.0 <= p && p <= 1.0);
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
        // even triangles have a downwards neighbor, while
        // odd triangles have an upwards neighbor
        bool even_site = (i+j) % 2 == 0;
        int tri_i_ud = wrap_i(even_site ? i+1 : i-1);
        int pet_i_ud = wrap_i(even_site ? i : i-1);
        bond_bit UD = even_site ? DOWN : UP;
        bond_bit FLIP_UD = even_site ? UP : DOWN;
        bool tri_ud = lattice[ind_tri(t, tri_i_ud, j)];

        // relevant petals
        int j_odd = j/2;
        bool pet_bwd_l = lattice[ind_pet_even(t, i, j_l)];
        bool pet_bwd_r = lattice[ind_pet_even(t, i, j)];
        bool pet_bwd_ud = lattice[ind_pet_odd(t, pet_i_ud, j_odd)];
        bool pet_fwd_l = lattice[ind_pet_even(t_fwd, i, j_l)];
        bool pet_fwd_r = lattice[ind_pet_even(t_fwd, i, j)];
        bool pet_fwd_ud = lattice[ind_pet_odd(t_fwd, pet_i_ud, j_odd)];
        bool pet_l_fixed = fixed[ind_pet_even(0, i, j_l)] != FREE;
        bool pet_r_fixed = fixed[ind_pet_even(0, i, j)] != FREE;
        bool pet_ud_fixed = fixed[ind_pet_odd(0, pet_i_ud, j_odd)] != FREE;

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
        if (!even_site) {
          continue;
        }
        // UP / DOWN
        if (ENABLE_TRI_SPATIAL) {
          bool delta_a = (tri == tri_ud);
          bool delta_b = (pet_fwd_ud == pet_bwd_ud);
          if (!pet_ud_fixed && tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << UD;
            bonds.tri[t][tri_i_ud][j] |= 1 << FLIP_UD;
          }
        }
        // RIGHT
        if (ENABLE_TRI_SPATIAL) {
          bool delta_a = (tri == tri_r);
          bool delta_b = (pet_fwd_r == pet_bwd_r);
          if (!pet_r_fixed && tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << RIGHT;
            bonds.tri[t][i][j_r] |= 1 << LEFT;
          }
        }
        // LEFT
        if (ENABLE_TRI_SPATIAL) {
          bool delta_a = (tri == tri_l);
          bool delta_b = (pet_fwd_l == pet_bwd_l);
          if (!pet_l_fixed && tri_bond_spatial(delta_a, delta_b)) {
            bonds.tri[t][i][j] |= 1 << LEFT;
            bonds.tri[t][i][j_l] |= 1 << RIGHT;
          }
        }
      }
    }
  }
}

// flood fill expanding from sites already pushed on the queue
void flood_fill_tri(spin_t spin) {
  while (!q_empty()) {
    coord_t x = q_pop();
    int t = x.t, i = x.i, j = x.j;
    assert(seen.tri[t][i][j]);
    lattice[ind_tri(t, i, j)] = spin;
    bond_dirs_t bond_dirs = bonds.tri[t][i][j];
    // TFWD
    if (bond_tfwd(bond_dirs)) {
      int tfwd = wrap_t(t+1);
      if (!seen.tri[tfwd][i][j]) {
        seen.tri[tfwd][i][j] = true;
        q_push((coord_t){tfwd, i, j});
      }
    }
    // TBWD
    if (bond_tbwd(bond_dirs)) {
      int tbwd = wrap_t(t-1);
      if (!seen.tri[tbwd][i][j]) {
        seen.tri[tbwd][i][j] = true;
        q_push((coord_t){tbwd, i, j});
      }
    }
    // LEFT
    if (bond_l(bond_dirs)) {
      int j_l = wrap_dj(j-1);
      if (!seen.tri[t][i][j_l]) {
        seen.tri[t][i][j_l] = true;
        q_push((coord_t){t, i, j_l});
      }
    }
    // RIGHT
    if (bond_r(bond_dirs)) {
      int j_r = wrap_dj(j+1);
      if (!seen.tri[t][i][j_r]) {
        seen.tri[t][i][j_r] = true;
        q_push((coord_t){t, i, j_r});
      }
    }
    // UP/DOWN
    if (bond_u(bond_dirs) || bond_d(bond_dirs)) {
      bool even_site = (i+j) % 2 == 0;
      int i_ud = wrap_i(even_site ? i+1 : i-1);
      if (!seen.tri[t][i_ud][j]) {
        seen.tri[t][i_ud][j] = true;
        q_push((coord_t){t, i_ud, j});
      }
    }
  }
}

void update_tri_spins() {
  memset(seen.tri, 0, sizeof(seen.tri));
  // outer loop, fixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.tri[t][i][j]) {
          continue;
        }
        if (fixed[ind_tri(0, i, j)] == FREE) {
          continue;
        }
        seen.tri[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j});
        flood_fill_tri(fixed[ind_tri(0, i, j)]);
        assert(q_empty());
      }
    }
  }
  // outer loop, unfixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.tri[t][i][j]) {
          continue;
        }
        seen.tri[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j});
        flood_fill_tri(rand_bool());
        assert(q_empty());
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

        // relevant triangles
        int j_l = j;
        int j_r = wrap_dj(j+1);
        bool tri_fwd_l = lattice[ind_tri(t, i, j_l)];
        bool tri_fwd_r = lattice[ind_tri(t, i, j_r)];

        // bonds
        // T_FWD / T_BWD
        {
          bool delta_b = (pet == pet_fwd);
          bool delta_a = (tri_fwd_l == tri_fwd_r);
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
        int j_even = 2*j + (i % 2);
        int j_even_r = j_even;
        int j_even_l = wrap_dj(j_even - 1);
        bool pet_ur = lattice[ind_pet_even(t, i_u, j_even_r)];
        bool pet_ul = lattice[ind_pet_even(t, i_u, j_even_l)];
        bool pet_dr = lattice[ind_pet_even(t, i_d, j_even_r)];
        bool pet_dl = lattice[ind_pet_even(t, i_d, j_even_l)];

        // relevant triangles
        bool tri_fwd_u = lattice[ind_tri(t, i_u, j_even)];
        bool tri_fwd_d = lattice[ind_tri(t, i_d, j_even)];
        bool tri_bwd_u = lattice[ind_tri(t_bwd, i_u, j_even)];
        bool tri_bwd_d = lattice[ind_tri(t_bwd, i_d, j_even)];
        bool tri_u_fixed = fixed[ind_tri(0, i_u, j_even)] != FREE;
        bool tri_d_fixed = fixed[ind_tri(0, i_d, j_even)] != FREE;

        // bonds
        // T_FWD / T_BWD
        {
          bool delta_b = (pet == pet_fwd);
          bool delta_a = (tri_fwd_u == tri_fwd_d);
          if (pet_bond_temporal(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << T_FWD;
            bonds.pet_odd[t_fwd][i][j] |= 1 << T_BWD;
          }
        }
        // UP (3-way)
        if (ENABLE_PET_SPATIAL) {
          bool delta_b = (pet == pet_ur && pet_ur == pet_ul);
          bool delta_a = (tri_fwd_u == tri_bwd_u);
          if (!tri_u_fixed && pet_bond_spatial(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << UP;
            bonds.pet_even[t][i_u][j_even_r] |= 1 << LEFT;
            bonds.pet_even[t][i_u][j_even_l] |= 1 << RIGHT;
          }
        }
        // DOWN (3-way)
        if (ENABLE_PET_SPATIAL) {
          bool delta_b = (pet == pet_dr && pet_dr == pet_dl);
          bool delta_a = (tri_fwd_d == tri_bwd_d);
          if (!tri_d_fixed && pet_bond_spatial(delta_a, delta_b)) {
            bonds.pet_odd[t][i][j] |= 1 << DOWN;
            bonds.pet_even[t][i_d][j_even_r] |= 1 << LEFT;
            bonds.pet_even[t][i_d][j_even_l] |= 1 << RIGHT;
          }
        }
      }
    }
  }
}

// flood fill expanding from sites already pushed on the queue
void flood_fill_pet(spin_t spin) {
  while (!q_empty()) {
    coord_t x = q_pop();
    int t = x.t, i = x.i, j = x.j;
    // even petals
    if (x.p == EVEN) {
      assert(seen.pet_even[t][i][j]);
      lattice[ind_pet_even(t, i, j)] = spin;
      bond_dirs_t bond_dirs = bonds.pet_even[t][i][j];
      // TFWD
      if (bond_tfwd(bond_dirs)) {
        int tfwd = wrap_t(t+1);
        if (!seen.pet_even[tfwd][i][j]) {
          seen.pet_even[tfwd][i][j] = true;
          q_push((coord_t){tfwd, i, j, EVEN});
        }
      }
      // TBWD
      if (bond_tbwd(bond_dirs)) {
        int tbwd = wrap_t(t-1);
        if (!seen.pet_even[tbwd][i][j]) {
          seen.pet_even[tbwd][i][j] = true;
          q_push((coord_t){tbwd, i, j, EVEN});
        }
      }
      // LEFT (3-way)
      if (bond_l(bond_dirs)) {
        int j_le = wrap_dj(j-1);
        int j_lo = j/2;
        int i_l = ((i+j) % 2 == 0) ? i : wrap_i(i-1);
        if (!seen.pet_even[t][i][j_le]) {
          seen.pet_even[t][i][j_le] = true;
          q_push((coord_t){t, i, j_le, EVEN});
        }
        if (!seen.pet_odd[t][i_l][j_lo]) {
          seen.pet_odd[t][i_l][j_lo] = true;
          q_push((coord_t){t, i_l, j_lo, ODD});
        }
      }
      // RIGHT (3-way)
      if (bond_r(bond_dirs)) {
        int j_re = wrap_dj(j+1);
        int j_ro = wrap_dj(j+1)/2;
        int i_r = ((i+j) % 2 == 0) ? wrap_i(i-1) : i;
        if (!seen.pet_even[t][i][j_re]) {
          seen.pet_even[t][i][j_re] = true;
          q_push((coord_t){t, i, j_re, EVEN});
        }
        if (!seen.pet_odd[t][i_r][j_ro]) {
          seen.pet_odd[t][i_r][j_ro] = true;
          q_push((coord_t){t, i_r, j_ro, ODD});
        }
      }
      assert(!bond_u(bond_dirs));
      assert(!bond_d(bond_dirs));
    }
    // odd petals
    else {
      assert(seen.pet_odd[t][i][j]);
      lattice[ind_pet_odd(t, i, j)] = spin;
      bond_dirs_t bond_dirs = bonds.pet_odd[t][i][j];
      // TFWD
      if (bond_tfwd(bond_dirs)) {
        int tfwd = wrap_t(t+1);
        if (!seen.pet_odd[tfwd][i][j]) {
          seen.pet_odd[tfwd][i][j] = true;
          q_push((coord_t){tfwd, i, j, ODD});
        }
      }
      // TBWD
      if (bond_tbwd(bond_dirs)) {
        int tbwd = wrap_t(t-1);
        if (!seen.pet_odd[tbwd][i][j]) {
          seen.pet_odd[tbwd][i][j] = true;
          q_push((coord_t){tbwd, i, j, ODD});
        }
      }
      // UP (3-way)
      if (bond_u(bond_dirs)) {
        bool even_row = i % 2 == 0;
        int j_r = wrap_dj(2*j + (even_row ? 0 : 1));
        int j_l = wrap_dj(j_r-1);
        if (!seen.pet_even[t][i][j_l]) {
          seen.pet_even[t][i][j_l] = true;
          q_push((coord_t){t, i, j_l, EVEN});
        }
        if (!seen.pet_even[t][i][j_r]) {
          seen.pet_even[t][i][j_r] = true;
          q_push((coord_t){t, i, j_r, EVEN});
        }
      }
      // DOWN (3-way)
      if (bond_d(bond_dirs)) {
        bool even_row = i % 2 == 0;
        int j_r = wrap_dj(2*j + (even_row ? 0 : 1));
        int j_l = wrap_dj(j_r-1);
        int i_d = wrap_i(i+1);
        if (!seen.pet_even[t][i_d][j_l]) {
          seen.pet_even[t][i_d][j_l] = true;
          q_push((coord_t){t, i_d, j_l, EVEN});
        }
        if (!seen.pet_even[t][i_d][j_r]) {
          seen.pet_even[t][i_d][j_r] = true;
          q_push((coord_t){t, i_d, j_r, EVEN});
        }
      }
      assert(!bond_l(bond_dirs));
      assert(!bond_r(bond_dirs));
    }
  }
}

void update_pet_spins() {
  memset(seen.pet_even, 0, sizeof(seen.pet_even));
  memset(seen.pet_odd, 0, sizeof(seen.pet_odd));
  // outer loop even, fixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.pet_even[t][i][j]) {
          continue;
        }
        if (fixed[ind_pet_even(0, i, j)] == FREE) {
          continue;
        }
        seen.pet_even[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, EVEN});
        flood_fill_pet(fixed[ind_pet_even(0, i, j)]);
        assert(q_empty());
      }
    }
  }
  // outer loop even, unfixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.pet_even[t][i][j]) {
          continue;
        }
        seen.pet_even[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, EVEN});
        flood_fill_pet(rand_bool());
        assert(q_empty());
      }
    }
  }
  // outer loop odd, fixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (seen.pet_odd[t][i][j]) {
          continue;
        }
        if (fixed[ind_pet_odd(0, i, j)] == FREE) {
          continue;
        }
        seen.pet_odd[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, ODD});
        flood_fill_pet(fixed[ind_pet_odd(0, i, j)]);
        assert(q_empty());
      }
    }
  }
  // outer loop odd, unfixed
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (seen.pet_odd[t][i][j]) {
          continue;
        }
        seen.pet_odd[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, ODD});
        flood_fill_pet(rand_bool());
        assert(q_empty());
      }
    }
  }
}

void measure_HP_HE(double* HP, double* HE) {
  int tot_HP = 0;
  int tot_HE = 0;
  // even petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (fixed[ind_pet_even(0, i, j)] != FREE) {
          continue;
        }
        // relevant petals
        bool pet = lattice[ind_pet_even(t, i, j)];
        int t_fwd = wrap_t(t+1);
        bool pet_fwd = lattice[ind_pet_even(t_fwd, i, j)];
        int j_l = j;
        int j_r = wrap_dj(j+1);

        // relevant triangles
        bool tri_l = lattice[ind_tri(t, i, j_l)];
        bool tri_r = lattice[ind_tri(t, i, j_r)];

        bool delta_b = (pet == pet_fwd);
        bool delta_a = (tri_l == tri_r);

        tot_HP += (!delta_b) && delta_a;
        tot_HE += 2*delta_a - 1;
      }
    }
  }
  // odd petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (fixed[ind_pet_odd(0, i, j)] != FREE) {
          continue;
        }
        // relevant petals
        bool pet = lattice[ind_pet_odd(t, i, j)];
        int t_fwd = wrap_t(t+1);
        bool pet_fwd = lattice[ind_pet_odd(t_fwd, i, j)];
        int i_u = i;
        int i_d = wrap_i(i+1);
        int j_even = 2*j + (i % 2);

        // relevant triangles
        bool tri_u = lattice[ind_tri(t, i_u, j_even)];
        bool tri_d = lattice[ind_tri(t, i_d, j_even)];

        bool delta_b = (pet == pet_fwd);
        bool delta_a = (tri_u == tri_d);

        tot_HP += (!delta_b) && delta_a;
        tot_HE += 2*delta_a - 1;
      }
    }
  }
  const int N_PETALS = NT * NROWS * 3*NCOLS;
  *HP = tot_HP / (double)N_PETALS;
  *HE = tot_HE / (double)N_PETALS;
}

void measure_HT(double* HT) {
  int tot_HT = 0;
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (fixed[ind_tri(0, i, j)] != FREE) {
          continue;
        }
        // relevant triangles
        bool tri = lattice[ind_tri(t, i, j)];
        int t_fwd = wrap_t(t+1);
        bool tri_fwd = lattice[ind_tri(t_fwd, i, j)];
        int j_l = wrap_dj(j-1);
        int pet_i_ud = wrap_i(((i+j) % 2 == 0) ? i : i-1);

        // relevant petals
        int j_odd = j/2;
        bool pet_fwd_l = lattice[ind_pet_even(t_fwd, i, j_l)];
        bool pet_fwd_r = lattice[ind_pet_even(t_fwd, i, j)];
        bool pet_fwd_ud = lattice[ind_pet_odd(t_fwd, pet_i_ud, j_odd)];

        bool delta_a = (tri == tri_fwd);
        bool delta_b = (pet_fwd_l == pet_fwd_r && pet_fwd_r == pet_fwd_ud);
        tot_HT += (!delta_a) && delta_b;
      }
    }
  }
  const int N_TRIANGLES = NT * NROWS * 2*NCOLS;
  *HT = tot_HT / (double)N_TRIANGLES;
}

void write_lattice(FILE *f) {
  assert(f != NULL);
  fwrite(lattice, 1, sizeof(lattice), f);
}

typedef enum {
  E_OK = 0,
  E_ARGS,
  E_OUT_FILE,
} error_t;

const size_t STRLEN = 256;

typedef struct {
  int n_iter;
  int save_freq;
  int meas_freq;
  unsigned seed;
  const char* prefix;
  // derived
  int n_meas;
  char fname_ens[STRLEN], fname_HT[STRLEN], fname_HP[STRLEN], fname_HE[STRLEN];
} config_t;

int parse_args(int argc, char** argv, config_t* cfg) {
  if (argc < 2) {
    printf("Usage: %s <prefix>\n", argv[0]);
    return E_ARGS;
  }
  const char* prefix = argv[1];
  size_t len = strlen(prefix);
  if (len >= STRLEN - 50) {
    printf("Invalid prefix, too long\n");
    return E_ARGS;
  }
  strncpy(cfg->fname_ens, prefix, STRLEN);
  strncpy(cfg->fname_HT, prefix, STRLEN);
  strncpy(cfg->fname_HP, prefix, STRLEN);
  strncpy(cfg->fname_HE, prefix, STRLEN);
  strcpy(cfg->fname_ens + len, ".ens.dat");
  strcpy(cfg->fname_HT + len, ".HT.dat");
  strcpy(cfg->fname_HP + len, ".HP.dat");
  strcpy(cfg->fname_HE + len, ".HE.dat");

  // TODO
  cfg->n_iter = 100000;
  cfg->save_freq = 100;
  cfg->meas_freq = 1;
  cfg->n_meas = (cfg->n_iter + 1) / cfg->meas_freq;
  return 0;
}

int main(int argc, char** argv) {
  config_t cfg;
  int ret = parse_args(argc, argv, &cfg);
  if (ret != E_OK) {
    return ret;
  }

  FILE *f = fopen(cfg.fname_ens, "wb");
  FILE *f_HT = fopen(cfg.fname_HT, "wb");
  FILE *f_HP = fopen(cfg.fname_HP, "wb");
  FILE *f_HE = fopen(cfg.fname_HE, "wb");
  if (f == NULL || f_HT == NULL || f_HP == NULL || f_HE == NULL) {
    printf("Failed to open output file\n");
    return E_OUT_FILE;
  }

  memset(lattice, 0xff, sizeof(lattice));
  memset(fixed, 0xff, sizeof(fixed));
  q_init();

  double* HT = malloc(cfg.n_meas * sizeof(double));
  double* HP = malloc(cfg.n_meas * sizeof(double));
  double* HE = malloc(cfg.n_meas * sizeof(double));
  srand(59263491);

  // init_fixed_pbc();
  init_fixed_obc();

  // init_lat_cold();
  init_lat_hot_tri();
  init_lat_apply_fixed();

  for (int i = 0; i < cfg.n_iter; ++i) {
    if ((i+1) % 1000 == 0) {
      printf("Iter %d / %d\n", i+1, cfg.n_iter);
    }
    // triangle sublattice
    sample_tri_bonds();
    update_tri_spins();
    // petal sublattice
    sample_pet_bonds();
    update_pet_spins();

    if ((i+1) % cfg.meas_freq == 0) {
      const int ind = ((i+1) / cfg.meas_freq) - 1;
      measure_HT(&HT[ind]);
      measure_HP_HE(&HP[ind], &HE[ind]);
    }

    if ((i+1) % cfg.save_freq == 0) {
      write_lattice(f);
    }
  }

  fwrite(HT, sizeof(double), cfg.n_meas, f_HT);
  fwrite(HP, sizeof(double), cfg.n_meas, f_HP);
  fwrite(HE, sizeof(double), cfg.n_meas, f_HE);

  free(HT);
  free(HP);
  free(HE);
  fclose(f);
  fclose(f_HT);
  fclose(f_HP);
  fclose(f_HE);
}
