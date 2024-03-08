#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


#ifndef NROWS
#define NROWS 2
#endif
#ifndef NCOLS
#define NCOLS 2
#endif
#ifndef NT
#define NT 64
#endif

// static double betaT;
// static double betaP;
// static double betaE;
// static double tanh_betaT;
// static double cosh_betaT;
// static double tanh_betaP;
// static double cosh_betaP;
// static double exp_m2_betaE;
static double p_tri_spatial;
static double p_tri_temporal;
static double p_pet_spatial;
static double p_pet_temporal;

void init_couplings(double dt, double KP, double KE) {
  // derived couplings betaX = dt K_x
  // double betaT = dt; // GK conventions
  // double betaP = dt * KP; // GK conventions
  double betaE = dt * KE;
  double betaT = dt * 4; // DB conventions
  double betaP = dt * KP * 2; // DB conventions
  // derived bond probs
  {
    double p = 1.0 - exp(-2*betaE) / cosh(betaP);
    assert(0.0 <= p && p <= 1.0);
    p_tri_spatial = p;
  }
  {
    double p = 1.0 - tanh(betaT);
    assert(0.0 <= p && p <= 1.0);
    p_tri_temporal = p;
  }
  {
    double p = 1.0 - 1.0 / cosh(betaT);
    assert(0.0 <= p && p <= 1.0);
    p_pet_spatial = p;
  }
  {
    double p = 1.0 - tanh(betaP);
    assert(0.0 <= p && p <= 1.0);
    p_pet_temporal = p;
  }
}

// settings for debugging
#define ENABLE_PET_SPATIAL 1
#define ENABLE_TRI_SPATIAL 1

// error codes
typedef enum {
  E_OK = 0,
  E_ARGS,
  E_OUT_FILE,
  E_VALUE,
} error_t;


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
#define DUMMY 0xaa
static spin_t geom[EROWS * ECOLS];


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


typedef enum {OBC, PBC} bc_t;
typedef enum {COLD, HOT} init_t;

static void init_geom_obc() {
  // fix left and bottom zones
  for (int x = 0; x < EROWS; ++x) {
    for (int y = 0; y < ECOLS; ++y) {
      int i = x*ECOLS + y;
      if (x == EROWS-1 || y >= ECOLS-3) {
        geom[i] = false;
      }
    }
  }
}

static void init_geom_pbc() {
  // no need to change anything
}

static int init_geom(bc_t bc_kind) {
  // mark hex lattice as FREE and other sites as DUMMY
  memset(geom, DUMMY, sizeof(geom));
  for (int i = 0; i < NROWS; ++i) {
    for (int j = 0; j < 2*NCOLS; ++j) {
      geom[ind_tri(0, i, j)] = FREE;
      geom[ind_pet_even(0, i, j)] = FREE;
    }
    for (int j = 0; j < NCOLS; ++j) {
      geom[ind_pet_odd(0, i, j)] = FREE;
    }
  }
  // fill in fixed sites based on BCs
  if (bc_kind == OBC) {
    init_geom_obc();
  }
  else if (bc_kind == PBC) {
    init_geom_pbc();
  }
  else {
    return E_VALUE;
  }
  return E_OK;
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

int init_lat(init_t init_kind) {
  memset(lattice, 0xff, sizeof(lattice));
  if (init_kind == COLD) {
    init_lat_cold();
  }
  else if (init_kind == HOT) {
    init_lat_hot_tri();
  }
  else {
    return E_VALUE;
  }
  return E_OK;
}

void init_lat_apply_geom() {
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        spin_t value = geom[ind_tri(0, i, j)];
        if (value != FREE) {
          lattice[ind_tri(t, i, j)] = value;
        }
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        spin_t value = geom[ind_pet_even(0, i, j)];
        if (value != FREE) {
          lattice[ind_pet_even(t, i, j)] = value;
        }
      }
    }
  }
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        spin_t value = geom[ind_pet_odd(0, i, j)];
        if (value != FREE) {
          lattice[ind_pet_odd(t, i, j)] = value;
        }
      }
    }
  }
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


bool tri_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  if (!delta_b) {
    return true;
  }
  return rand_double() < p_tri_spatial;
}

bool tri_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_a) {
    return false;
  }
  if (!delta_b) {
    return true;
  }
  return rand_double() < p_tri_temporal;
}

bool pet_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  if (!delta_a) {
    return true;
  }
  return rand_double() < p_pet_spatial;
}

bool pet_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_b) {
    return false;
  }
  if (!delta_a) {
    return true;
  }
  return rand_double() < p_pet_temporal;
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
        bool pet_l_fixed = geom[ind_pet_even(0, i, j_l)] != FREE;
        bool pet_r_fixed = geom[ind_pet_even(0, i, j)] != FREE;
        bool pet_ud_fixed = geom[ind_pet_odd(0, pet_i_ud, j_odd)] != FREE;

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
        if (geom[ind_tri(0, i, j)] == FREE) {
          continue;
        }
        seen.tri[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j});
        flood_fill_tri(geom[ind_tri(0, i, j)]);
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
        bool tri_u_fixed = geom[ind_tri(0, i_u, j_even)] != FREE;
        bool tri_d_fixed = geom[ind_tri(0, i_d, j_even)] != FREE;

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
        if (geom[ind_pet_even(0, i, j)] == FREE) {
          continue;
        }
        seen.pet_even[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, EVEN});
        flood_fill_pet(geom[ind_pet_even(0, i, j)]);
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
        if (geom[ind_pet_odd(0, i, j)] == FREE) {
          continue;
        }
        seen.pet_odd[t][i][j] = true;
        assert(q_empty());
        q_push((coord_t){t, i, j, ODD});
        flood_fill_pet(geom[ind_pet_odd(0, i, j)]);
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

void measure_HP_HE(int* HP_plus, int* HP_minus, int* HE_plus, int* HE_minus) {
  int tot_HP_plus = 0;
  int tot_HP_minus = 0;
  int tot_HE_plus = 0;
  int tot_HE_minus = 0;
  // even petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (geom[ind_pet_even(0, i, j)] != FREE) {
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

        tot_HP_plus += delta_b && delta_a;
        tot_HP_minus += (!delta_b) && delta_a;
        tot_HE_plus += delta_a;
        tot_HE_minus += !delta_a;
      }
    }
  }
  // odd petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (geom[ind_pet_odd(0, i, j)] != FREE) {
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

        tot_HP_plus += delta_b && delta_a;
        tot_HP_minus += (!delta_b) && delta_a;
        tot_HE_plus += delta_a;
        tot_HE_minus += !delta_a;
      }
    }
  }
  *HP_plus = tot_HP_plus;
  *HP_minus = tot_HP_minus;
  *HE_plus = tot_HE_plus;
  *HE_minus = tot_HE_minus;
}

void measure_HT(int* HT_plus, int* HT_minus) {
  int tot_HT_plus = 0;
  int tot_HT_minus = 0;
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (geom[ind_tri(0, i, j)] != FREE) {
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
        tot_HT_plus += delta_a && delta_b;
        tot_HT_minus += (!delta_a) && delta_b;
      }
    }
  }
  *HT_plus = tot_HT_plus;
  *HT_minus = tot_HT_minus;
}

void write_lattice(FILE *f) {
  assert(f != NULL);
  fwrite(lattice, 1, sizeof(lattice), f);
}

// const size_t STRLEN = 256;
#define STRLEN 256

typedef struct {
  double dt, KP, KE;
  int n_iter;
  int save_freq;
  int meas_freq;
  unsigned seed;
  const char* prefix;
  // derived
  int n_meas;
  char fname_ens[STRLEN], fname_meta[STRLEN],
    fname_HT[STRLEN], fname_HP[STRLEN], fname_HE[STRLEN];
} config_t;

void usage(const char* prog) {
  printf("Usage: %s -t <dt> -p <KP> -e <KE> -f <out_prefix> "
         "-i <n_iter> -s <save_freq> -m <meas_freq> "
         "[-b (obc|pbc)] [-c (cold|hot)] [-r <seed>]\n", prog);
}

int parse_args(int argc, char** argv, config_t* cfg) {
  if (argc < 2) {
    usage(argv[0]);
    return E_ARGS;
  }

  bool set_dt = false, set_KP = false, set_KE = false,
      set_prefix = false, set_n_iter = false, set_save_freq = false,
      set_meas_freq = false;
  double dt, KP, KE;
  bc_t bc_kind = PBC;
  init_t init_kind = HOT;
  const char* prefix;
  char c;
  while ((c = getopt(argc, argv, "t:p:e:f:i:s:m:b:c:r:")) != -1) {
    if (c == 't') {
      dt = atof(optarg);
      set_dt = true;
    }
    else if (c == 'p') {
      KP = atof(optarg);
      set_KP = true;
    }
    else if (c == 'e') {
      KE = atof(optarg);
      set_KE = true;
    }
    else if (c == 'f') {
      prefix = optarg;
      set_prefix = true;
    }
    else if (c == 'i') {
      cfg->n_iter = atoi(optarg);
      set_n_iter = true;
    }
    else if (c == 's') {
      cfg->save_freq = atoi(optarg);
      set_save_freq = true;
    }
    else if (c == 'm') {
      cfg->meas_freq = atoi(optarg);
      set_meas_freq = true;
    }
    else if (c == 'r') {
      cfg->seed = atoi(optarg);
    }
    else if (c == 'b') {
      if (strcmp(optarg, "pbc") == 0) {
        bc_kind = PBC;
      }
      else if (strcmp(optarg, "obc") == 0) {
        bc_kind = OBC;
      }
      else {
        usage(argv[0]);
        printf("Invalid -b argument: %s\n", optarg);
        return E_ARGS;
      }
    }
    else if (c == 'c') {
      if (strcmp(optarg, "cold") == 0) {
        init_kind = COLD;
      }
      else if (strcmp(optarg, "hot") == 0) {
        init_kind = HOT;
      }
      else {
        usage(argv[0]);
        printf("Invalid -c argument: %s\n", optarg);
        return E_ARGS;
      }
    }
    else if (c == '?') {
      usage(argv[0]);
      printf("Missing arg or invalid -%c\n", optopt);
      return E_ARGS;
    }
    else {
      usage(argv[0]);
      return E_ARGS;
    }
  }

  // validate args
  if (!set_dt || !set_KP || !set_KE || !set_prefix ||
      !set_n_iter || !set_save_freq || !set_meas_freq) {
    usage(argv[0]);
    printf("Missing a flag\n");
    return E_ARGS;
  }

  if (cfg->meas_freq <= 0 || cfg->save_freq <= 0) {
    printf("Invalid meas_freq or save_freq\n");
    return E_ARGS;
  }
  if (dt <= 0.0) {
    printf("Invalid dt = %f\n", dt);
    return E_ARGS;
  }

  size_t len = strlen(prefix);
  if (len >= STRLEN - 50) {
    printf("Invalid prefix, too long\n");
    return E_ARGS;
  }

  // derived
  strncpy(cfg->fname_ens, prefix, STRLEN);
  strncpy(cfg->fname_meta, prefix, STRLEN);
  strncpy(cfg->fname_HT, prefix, STRLEN);
  strncpy(cfg->fname_HP, prefix, STRLEN);
  strncpy(cfg->fname_HE, prefix, STRLEN);
  strncpy(cfg->fname_ens + len, ".ens.dat", STRLEN-len);
  strncpy(cfg->fname_meta + len, ".meta.dat", STRLEN-len);
  strncpy(cfg->fname_HT + len, ".HT.dat", STRLEN-len);
  strncpy(cfg->fname_HP + len, ".HP.dat", STRLEN-len);
  strncpy(cfg->fname_HE + len, ".HE.dat", STRLEN-len);

  // derived
  if (cfg->meas_freq > 1) {
    cfg->n_meas = (cfg->n_iter + 1) / cfg->meas_freq;
  }
  else {
    assert(cfg->meas_freq == 1);
    cfg->n_meas = cfg->n_iter;
  }

  // init global state
  int ret;
  srand(cfg->seed);
  cfg->dt = dt;
  cfg->KP = KP;
  cfg->KE = KE;
  init_couplings(dt, KP, KE);
  ret = init_geom(bc_kind);
  assert(ret == E_OK);
  ret = init_lat(init_kind);
  assert(ret == E_OK);
  init_lat_apply_geom();

  return E_OK;
}

typedef struct {
  double dt, KP, KE;
  int nt, erows, ecols;
  unsigned seed;
  spin_t geom[EROWS * ECOLS];
} meta_t;

int main(int argc, char** argv) {
  config_t cfg;
  cfg.seed = time(NULL);
  q_init();

  int ret = parse_args(argc, argv, &cfg);
  if (ret != E_OK) {
    return ret;
  }

  FILE *f = fopen(cfg.fname_ens, "wb");
  FILE *f_meta = fopen(cfg.fname_meta, "wb");
  FILE *f_HT = fopen(cfg.fname_HT, "wb");
  FILE *f_HP = fopen(cfg.fname_HP, "wb");
  FILE *f_HE = fopen(cfg.fname_HE, "wb");
  if (f == NULL || f_meta == NULL || f_HT == NULL || f_HP == NULL || f_HE == NULL) {
    printf("Failed to open output file\n");
    return E_OUT_FILE;
  }
  meta_t meta = {
    .dt = cfg.dt, .KP = cfg.KP, .KE = cfg.KE,
    .nt = NT, .erows = EROWS, .ecols = ECOLS,
    .seed = cfg.seed
  };
  assert(sizeof(meta.geom) == sizeof(geom));
  memcpy(meta.geom, geom, sizeof(geom));
  fwrite(&meta, 1, sizeof(meta), f_meta);
  fclose(f_meta);


  // Hamiltonian terms stored as pairs [(H+, H-)[0], (H+, H-)[1], ...]
  size_t len_Hi = 2 * cfg.n_meas;
  int* HT = malloc(len_Hi * sizeof(int));
  int* HP = malloc(len_Hi * sizeof(int));
  int* HE = malloc(len_Hi * sizeof(int));

  clock_t start = clock();
  for (int i = 0; i < cfg.n_iter; ++i) {
    // progress
    if ((i+1) % 1000 == 0) {
      clock_t ticks = clock() - start;
      double elapsed = ticks / (double)CLOCKS_PER_SEC;
      double expected = (elapsed * cfg.n_iter) / (i+1);
      double rate = (i+1) / elapsed;
      printf(
          "Iter %d / %d (%.2f / %.2fs | %.2f it/s)\n",
          i+1, cfg.n_iter, elapsed, expected, rate);
    }
    // triangle sublattice
    sample_tri_bonds();
    update_tri_spins();
    // petal sublattice
    sample_pet_bonds();
    update_pet_spins();

    if ((i+1) % cfg.meas_freq == 0) {
      const int ind = ((i+1) / cfg.meas_freq) - 1;
      measure_HT(&HT[2*ind], &HT[2*ind+1]);
      measure_HP_HE(&HP[2*ind], &HP[2*ind+1], &HE[2*ind], &HE[2*ind+1]);
    }

    if ((i+1) % cfg.save_freq == 0) {
      write_lattice(f);
    }
  }

  fwrite(HT, sizeof(int), len_Hi, f_HT);
  fwrite(HP, sizeof(int), len_Hi, f_HP);
  fwrite(HE, sizeof(int), len_Hi, f_HE);

  free(HT);
  free(HP);
  free(HE);
  fclose(f);
  fclose(f_HT);
  fclose(f_HP);
  fclose(f_HE);
}
