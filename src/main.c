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
#error "NROWS must be defined"
#endif
#ifndef NCOLS
#error "NCOLS must be defined"
#endif
#ifndef NT
#error "NT must be defined"
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
  double betaT = dt; // GK conventions
  double betaP = dt * KP; // GK conventions
  double betaE = dt * KE;
  // double betaT = dt * 4; // DB conventions
  // double betaP = dt * KP * 2; // DB conventions
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

/// Dirac strings can be introduced on petals only.
///  - for odd petals they modify the interaction with the triangle above
///  - for even petals they modify the interaction with the left triangle
/// In either case, the petal and triangle see each other as the opposite
/// of their value.

// Which odd petals carry a Dirac string (other sites ignored)
static uint8_t dirac_str[EROWS * ECOLS];

// Gauss Law defects, just for tracking and validation
static int defects[NROWS * NCOLS];

static uint64_t mag_accum[EROWS * ECOLS];
static uint64_t mag_n = 0;
static int64_t Ex_accum[EROWS * ECOLS];
static uint64_t Ex_n = 0;
static uint64_t HT_plus_accum[EROWS * ECOLS];
static uint64_t HT_minus_accum[EROWS * ECOLS];
static uint64_t HP_plus_accum[EROWS * ECOLS];
static uint64_t HP_minus_accum[EROWS * ECOLS];
static uint64_t HE_plus_accum[EROWS * ECOLS];
static uint64_t HE_minus_accum[EROWS * ECOLS];
static uint64_t H_n = 0;


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

static inline int wrap_mod4(int x) {
  return ((x%4)+6)%4 - 2;
}

// DEBUG: measure gauss law violations across the lattice
static void measure_Gx(bool verify) {
  bool fail = false;
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        int djA = 2*j, djB = wrap_dj(djA+1), djC = wrap_dj(djA+2);
        int sjA = j, sjB = wrap_sj(j+1);
        int i_fwd = wrap_i(i+1);
        if (i % 2 == 1) {
          djA = wrap_dj(djA+1);
          djB = wrap_dj(djB+1);
          djC = wrap_dj(djC+1);
        }
        // triangles are {0,2}, petals are {1,3}
        int s[12] = {
          2*lattice[ind_tri(t, i, djA)],
          -1+2*lattice[ind_pet_even(t, i, djA)],
          2*lattice[ind_tri(t, i, djB)],
          -1+2*lattice[ind_pet_even(t, i, djB)],
          2*lattice[ind_tri(t, i, djC)],
          -1+2*lattice[ind_pet_odd(t, i, sjB)],
          2*lattice[ind_tri(t, i_fwd, djC)],
          -1+2*lattice[ind_pet_even(t, i_fwd, djB)],
          2*lattice[ind_tri(t, i_fwd, djB)],
          -1+2*lattice[ind_pet_even(t, i_fwd, djA)],
          2*lattice[ind_tri(t, i_fwd, djA)],
          -1+2*lattice[ind_pet_odd(t, i, sjA)]
        };
        int dirac[12] = {
          0, dirac_str[ind_pet_even(t, i, djA)],
          0, dirac_str[ind_pet_even(t, i, djB)],
          0, dirac_str[ind_pet_odd(t, i, sjB)],
          0, dirac_str[ind_pet_even(t, i_fwd, djB)],
          0, dirac_str[ind_pet_even(t, i_fwd, djA)],
          0, dirac_str[ind_pet_odd(t, i, sjA)]
        };
        // printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
        //        s[0], s[1], s[2], s[3],
        //        s[4], s[5], s[6], s[7],
        //        s[8], s[9], s[10], s[11]);
        int G = 0;
        for (int k = 0; k < 12; ++k) {
          // wrap difference into {-1, 1}
          G += wrap_mod4(s[k] - s[(k+1)%12] + dirac[k]);
        }
        if (t == 0) {
          if (!verify) {
            printf("G[%d,%d] = %d\n", i, j, G);
            defects[i*NCOLS + j] = G;
          }
        }
        if (defects[i*NCOLS + j] != G) {
          fail = true;
        }
      }
    }
  }
  assert(!fail);
}



typedef enum {OBC, PBC, RHOMB, RHOMB_STRING, RHOMB2_STRING, GEOM_FILE} bc_t;
typedef enum {COLD, HOT, COLD_STRING, INIT_FILE} init_t;

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

static void init_geom_rhomb(bool with_string) {
  assert(NROWS % 2 == 0);
  for (int x = 0; x < EROWS; ++x) {
    for (int y = 0; y < ECOLS; ++y) {
      int i = x*ECOLS + y;
      // trim top right tail
      if (x == 0 && y >= ECOLS-7) {
        geom[i] = with_string;
      }
      // trim bottom left tail
      else if (x == EROWS-4 && y < 2) {
        geom[i] = false;
      }
      // last 3 rows clear
      else if (x >= EROWS-3) {
        if (x >= EROWS-1) {
          geom[i] = with_string;
        }
        else {
          geom[i] = false;
        }
      }
      // last 5 cols clear
      else if (y >= ECOLS-5) {
        if (y >= ECOLS-3) {
          geom[i] = false;
        }
        else {
          geom[i] = with_string;
        }
      }
    }
  }
}

static void init_geom_rhomb2(bool with_string) {
  assert(NROWS % 2 == 0);
  for (int x = 0; x < EROWS; ++x) {
    for (int y = 0; y < ECOLS; ++y) {
      int i = x*ECOLS + y;
      if (x == 0 || y >= ECOLS-2) {
        geom[i] = with_string;
      }
      if (x >= EROWS-2 || y == 0) {
        geom[i] = false;
      }
    }
  }
}

static void init_geom_pbc() {
  // no need to change anything
}

static void init_geom_file(const char* fname) {
  FILE *f = fopen(fname, "rb");
  if (f == NULL) {
    printf("Failed to open geom file %s\n", fname);
    abort();
  }
  size_t geom_size = sizeof(geom)/sizeof(spin_t);
  size_t nelts = fread(geom, sizeof(spin_t), geom_size, f);
  fclose(f);
  if (nelts != geom_size) {
    printf("Unable to load geom file %s\n", fname);
    abort();
  }
  // check that all values are valid
  bool valid = true;
  for (size_t i = 0; i < geom_size; ++i) {
    if (geom[i] != FREE && geom[i] != DUMMY && geom[i] >= 2) {
      printf("Invalid geom %d (site %lu)\n", geom[i], i);
      valid = false;
    }
  }
  if (!valid) {
    abort();
  }
}

static int init_geom(bc_t bc_kind, const char* fname) {
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
  else if (bc_kind == RHOMB) {
    init_geom_rhomb(false);
  }
  else if (bc_kind == RHOMB_STRING) {
    init_geom_rhomb(true);
  }
  else if (bc_kind == RHOMB2_STRING) {
    init_geom_rhomb2(true);
  }
  else if (bc_kind == GEOM_FILE) {
    init_geom_file(fname);
  }
  else {
    return E_VALUE;
  }
  return E_OK;
}

static void print_geom() {
  for (int i = 0; i < EROWS; ++i) {
    for (int j = 0; j < ECOLS; ++j) {
      spin_t g = geom[i*ECOLS + j];
      char c = (g == DUMMY) ? ' ' : '.';
      printf("%c", c);
    }
    printf("\n");
  }
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

// add a string to a cold-init lattice by placing a defect-antidefect pair
// at positions
//  - (NROWS/2, (NCOLS+string_sep)/2)
//  - (NROWS/2, (NCOLS-string_sep)/2)
void init_lat_add_string(int string_sep) {
  printf("== BEFORE STRING ==\n");
  measure_Gx(false);
  int i = NROWS/2;
  if (i % 2 == 1) {
    i += 1;
  }
  // dual lattice column indices are 2x primal
  int jm = NCOLS - string_sep;
  int jp = NCOLS + string_sep;
  for (int t = 0; t < NT; ++t) {
    for (int j = jm; j < jp-1; ++j) {
      lattice[ind_tri(t, i, j)] = true;
      if (j % 2 == 0) {
        lattice[ind_pet_odd(t, i, j/2)] = true;
      }
    }
  }
  printf("== AFTER STRING ==\n");
  measure_Gx(false);
}

// TODO: needs to be updated to satisfy Gauss' Law
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

// file init uses a time-independent array
static spin_t init_file_buffer[EROWS * ECOLS];
void init_lat_file(const char* fname) {
  FILE *f = fopen(fname, "rb");
  if (f == NULL) {
    printf("Failed to open lattice init file %s\n", fname);
    abort();
  }
  size_t init_size = sizeof(init_file_buffer)/sizeof(spin_t);
  size_t nelts = fread(init_file_buffer, sizeof(spin_t), init_size, f);
  fclose(f);
  if (nelts != init_size) {
    printf("Unable to load lattice init file %s\n", fname);
    abort();
  }
  // check that init matches geom
  bool valid = true;
  assert(sizeof(geom) == sizeof(init_file_buffer));
  for (size_t i = 0; i < init_size; ++i) {
    if (geom[i] == DUMMY) {
      if (init_file_buffer[i] != 0xff) {
        printf("Init file set dummy site %lu to %d\n", i, init_file_buffer[i]);
        valid = false;
      }
    }
    else if (geom[i] == FREE) {
      if (init_file_buffer[i] >= 2) {
        printf("Init file set site %lu to bad value %d\n", i, init_file_buffer[i]);
        valid = false;
      }
    }
    else {
      assert(geom[i] <= 2);
      if (init_file_buffer[i] == 0xff) {
        init_file_buffer[i] = geom[i];
      }
      if (init_file_buffer[i] != geom[i]) {
        printf(
            "Init file overwrote fixed geom site %lu from %d to %d\n",
            i, geom[i], init_file_buffer[i]);
        valid = false;
      }
    }
  }
  if (!valid) {
    abort();
  }
  // write into lattice across all timeslices
  assert(sizeof(lattice) == NT * sizeof(init_file_buffer));
  for (int t = 0; t < NT; ++t) {
    memcpy(&lattice[t * init_size], init_file_buffer, sizeof(init_file_buffer));
  }
}

int init_lat(init_t init_kind, int string_sep, const char* fname) {
  memset(lattice, 0xff, sizeof(lattice));
  memset(mag_accum, 0.0, sizeof(mag_accum));
  memset(Ex_accum, 0.0, sizeof(Ex_accum));
  memset(HT_plus_accum, 0.0, sizeof(HT_plus_accum));
  memset(HT_minus_accum, 0.0, sizeof(HT_minus_accum));
  memset(HP_plus_accum, 0.0, sizeof(HP_plus_accum));
  memset(HP_minus_accum, 0.0, sizeof(HP_minus_accum));
  memset(HE_plus_accum, 0.0, sizeof(HE_plus_accum));
  memset(HE_minus_accum, 0.0, sizeof(HE_minus_accum));
  if (init_kind == COLD) {
    init_lat_cold();
  }
  else if (init_kind == HOT) {
    // TODO
    printf("HOT init not supported yet\n");
    abort();
    init_lat_hot_tri();
  }
  else if (init_kind == COLD_STRING) {
    init_lat_cold();
    init_lat_add_string(string_sep);
  }
  else if (init_kind == INIT_FILE) {
    init_lat_file(fname);
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

static int init_dirac_str(const char* fname) {
  if (fname == NULL) {
    memset(dirac_str, 0, sizeof(dirac_str));
    return E_OK;
  }
  FILE *f = fopen(fname, "rb");
  if (f == NULL) {
    printf("Failed to open dirac str file %s\n", fname);
    abort();
  }
  size_t size = sizeof(dirac_str)/sizeof(uint8_t);
  size_t nelts = fread(dirac_str, sizeof(uint8_t), size, f);
  fclose(f);
  if (nelts != size) {
    printf("Unabled to load dirac str file %s\n", fname);
    abort();
  }
  return E_OK;
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

// which directions each point is bonded in for cluster formation
static struct {
  bond_dirs_t tri[NT][NROWS][2*NCOLS];
  bond_dirs_t pet_even[NT][NROWS][2*NCOLS];
  bond_dirs_t pet_odd[NT][NROWS][NCOLS];
} bonds;

// how many times we visited each node already when flipping or unflipping
static struct {
  int tri[NT][NROWS][2*NCOLS];
  int pet_even[NT][NROWS][2*NCOLS];
  int pet_odd[NT][NROWS][NCOLS];
} seen;

// which temporal sheet we are on (to avoid wrapping for Gauss' Law)
static struct {
  int tri[NT][NROWS][2*NCOLS];
  int pet_even[NT][NROWS][2*NCOLS];
  int pet_odd[NT][NROWS][NCOLS];
} sheet;


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
    assert(delta_b);
    return false;
  }
  if (!delta_b) {
    return true;
  }
  return rand_double() < p_tri_spatial;
}

bool tri_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_a) {
    assert(delta_b);
    return false;
  }
  if (!delta_b) {
    return true;
  }
  return rand_double() < p_tri_temporal;
}

bool pet_bond_spatial(bool delta_a, bool delta_b) {
  if (!delta_b) {
    assert(delta_a);
    return false;
  }
  if (!delta_a) {
    return true;
  }
  return rand_double() < p_pet_spatial;
}

bool pet_bond_temporal(bool delta_a, bool delta_b) {
  if (!delta_b) {
    assert(delta_a);
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
        // dirac string between tri and petals to the right or below
        bool has_dirac_str_r = dirac_str[ind_pet_even(0, i, j)];
        bool has_dirac_str_ud = dirac_str[ind_pet_odd(0, pet_i_ud, j_odd)] && (UD == DOWN);
        bool pet_bwd_l = lattice[ind_pet_even(t, i, j_l)];
        bool pet_bwd_r = (lattice[ind_pet_even(t, i, j)] + has_dirac_str_r) % 2;
        bool pet_bwd_ud = (lattice[ind_pet_odd(t, pet_i_ud, j_odd)] + has_dirac_str_ud) % 2;
        bool pet_fwd_l = lattice[ind_pet_even(t_fwd, i, j_l)];
        bool pet_fwd_r = (lattice[ind_pet_even(t_fwd, i, j)] + has_dirac_str_r) % 2;
        bool pet_fwd_ud = (lattice[ind_pet_odd(t_fwd, pet_i_ud, j_odd)] + has_dirac_str_ud) % 2;
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
void flood_fill_tri(spin_t spin, int gen, bool *temporal_wrap, int *count) {
  *temporal_wrap = false;
  *count = 0;
  while (!q_empty()) {
    *count += 1;
    coord_t x = q_pop();
    int t = x.t, i = x.i, j = x.j;
    assert(seen.tri[t][i][j] == gen);
    int prev_sheet = sheet.tri[t][i][j];
    lattice[ind_tri(t, i, j)] = spin;
    bond_dirs_t bond_dirs = bonds.tri[t][i][j];
    // TFWD
    if (bond_tfwd(bond_dirs)) {
      int tfwd = wrap_t(t+1);
      int next_sheet = (t+1 >= NT) ? prev_sheet+1 : prev_sheet;
      assert(bond_tbwd(bonds.tri[tfwd][i][j]));
      if (seen.tri[tfwd][i][j] < gen) {
        seen.tri[tfwd][i][j] = gen;
        sheet.tri[tfwd][i][j] = next_sheet;
        q_push((coord_t){tfwd, i, j, 0});
      }
      else if (sheet.tri[tfwd][i][j] != next_sheet) {
        *temporal_wrap = true;
      }
    }
    // TBWD
    if (bond_tbwd(bond_dirs)) {
      int tbwd = wrap_t(t-1);
      int next_sheet = (t-1 < 0) ? prev_sheet-1 : prev_sheet;
      assert(bond_tfwd(bonds.tri[tbwd][i][j]));
      if (seen.tri[tbwd][i][j] < gen) {
        seen.tri[tbwd][i][j] = gen;
        sheet.tri[tbwd][i][j] = next_sheet;
        q_push((coord_t){tbwd, i, j, 0});
      }
      else if (sheet.tri[tbwd][i][j] != next_sheet) {
        *temporal_wrap = true;
      }
    }
    // LEFT
    if (bond_l(bond_dirs)) {
      int j_l = wrap_dj(j-1);
      assert(bond_r(bonds.tri[t][i][j_l]));
      if (seen.tri[t][i][j_l] < gen) {
        seen.tri[t][i][j_l] = gen;
        sheet.tri[t][i][j_l] = prev_sheet;
        q_push((coord_t){t, i, j_l, 0});
      }
      else if (sheet.tri[t][i][j_l] != prev_sheet) {
        *temporal_wrap = true;
      }
    }
    // RIGHT
    if (bond_r(bond_dirs)) {
      int j_r = wrap_dj(j+1);
      assert(bond_l(bonds.tri[t][i][j_r]));
      if (seen.tri[t][i][j_r] < gen) {
        seen.tri[t][i][j_r] = gen;
        sheet.tri[t][i][j_r] = prev_sheet;
        q_push((coord_t){t, i, j_r, 0});
      }
      else if (sheet.tri[t][i][j_r] != prev_sheet) {
        *temporal_wrap = true;
      }
    }
    // UP/DOWN
    if (bond_u(bond_dirs) || bond_d(bond_dirs)) {
      bool even_site = (i+j) % 2 == 0;
      int i_ud = wrap_i(even_site ? i+1 : i-1);
      if (even_site) {
        assert(bond_d(bond_dirs));
        assert(bond_u(bonds.tri[t][i_ud][j]));
      }
      else {
        assert(bond_u(bond_dirs));
        assert(bond_d(bonds.tri[t][i_ud][j]));
      }
      if (seen.tri[t][i_ud][j] < gen) {
        seen.tri[t][i_ud][j] = gen;
        sheet.tri[t][i_ud][j] = prev_sheet;
        q_push((coord_t){t, i_ud, j, 0});
      }
      else if (sheet.tri[t][i_ud][j] != prev_sheet) {
        *temporal_wrap = true;
      }
    }
  }
}

double update_tri_spins() {
  memset(seen.tri, 0, sizeof(seen.tri));
  memset(sheet.tri, 0, sizeof(sheet.tri));
  bool temporal_wrap;
  int count;
  int acc = 0, rej = 0;
  // outer loop, fixed
  // TODO: do we need to run this loop?
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.tri[t][i][j]) {
          continue;
        }
        if (geom[ind_tri(0, i, j)] == FREE) {
          continue;
        }
        int gen = 1;
        seen.tri[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, 0});
        flood_fill_tri(geom[ind_tri(0, i, j)], gen, &temporal_wrap, &count);
        assert(q_empty());
        rej += count;
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
        int gen = 1;
        seen.tri[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, 0});
        spin_t prev_val = lattice[ind_tri(t, i, j)];
        flood_fill_tri(rand_bool(), gen, &temporal_wrap, &count);
        // printf("New cluster (%d)\n", count);
        assert(q_empty());
        if (temporal_wrap) { // unflip if wrapped
          // printf("temporal wrap, unflipping!\n");
          int prev_count = count;
          gen += 1;
          seen.tri[t][i][j] = gen;
          q_push((coord_t){t, i, j, 0});
          flood_fill_tri(prev_val, gen, &temporal_wrap, &count);
          assert(q_empty());
          assert(count == prev_count);
          rej += count;
        }
        else {
          acc += count;
        }
      }
    }
  }
  return acc / (double)(acc + rej);
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
        bool has_dirac_str = dirac_str[ind_pet_even(0, i, j)];
        int t_fwd = wrap_t(t+1);
        bool pet_fwd = lattice[ind_pet_even(t_fwd, i, j)];

        // relevant triangles (LEFT triangles shifted by dirac str)
        int j_l = j;
        int j_r = wrap_dj(j+1);
        bool tri_fwd_l = (lattice[ind_tri(t, i, j_l)] + has_dirac_str) % 2;
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
        bool has_dirac_str = dirac_str[ind_pet_odd(0, i, j)];
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

        // relevant triangles (UP triangles shifted by dirac str)
        bool tri_fwd_u = (
            lattice[ind_tri(t, i_u, j_even)] + has_dirac_str) % 2;
        bool tri_fwd_d = lattice[ind_tri(t, i_d, j_even)];
        bool tri_bwd_u = (
            lattice[ind_tri(t_bwd, i_u, j_even)] + has_dirac_str) % 2;
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
void flood_fill_pet(spin_t spin, int gen, bool *temporal_wrap, int *count) {
  *temporal_wrap = false;
  *count = 0;
  while (!q_empty()) {
    *count += 1;
    coord_t x = q_pop();
    int t = x.t, i = x.i, j = x.j;
    // even petals
    if (x.p == EVEN) {
      assert(seen.pet_even[t][i][j] == gen);
      int prev_sheet = sheet.pet_even[t][i][j];
      lattice[ind_pet_even(t, i, j)] = spin;
      // printf("    (%d,%d,%d)\n",
      //        ind_pet_even(t, i, j)/(EROWS*ECOLS),
      //        (ind_pet_even(t, i, j)%(EROWS*ECOLS))/ECOLS,
      //        ind_pet_even(t, i, j)%(ECOLS));
      bond_dirs_t bond_dirs = bonds.pet_even[t][i][j];
      // TFWD
      if (bond_tfwd(bond_dirs)) {
        int tfwd = wrap_t(t+1);
        int next_sheet = (t+1 >= NT) ? prev_sheet+1 : prev_sheet;
        assert(bond_tbwd(bonds.pet_even[tfwd][i][j]));
        if (seen.pet_even[tfwd][i][j] < gen) {
          seen.pet_even[tfwd][i][j] = gen;
          sheet.pet_even[tfwd][i][j] = next_sheet;
          q_push((coord_t){tfwd, i, j, EVEN});
        }
        else if (sheet.pet_even[tfwd][i][j] != next_sheet) {
          *temporal_wrap = true;
        }
      }
      // TBWD
      if (bond_tbwd(bond_dirs)) {
        int tbwd = wrap_t(t-1);
        int next_sheet = (t-1 < 0) ? prev_sheet-1 : prev_sheet;
        assert(bond_tfwd(bonds.pet_even[tbwd][i][j]));
        if (seen.pet_even[tbwd][i][j] < gen) {
          seen.pet_even[tbwd][i][j] = gen;
          sheet.pet_even[tbwd][i][j] = next_sheet;
          q_push((coord_t){tbwd, i, j, EVEN});
        }
        else if (sheet.pet_even[tbwd][i][j] != next_sheet) {
          *temporal_wrap = true;
        }
      }
      // LEFT (3-way)
      if (bond_l(bond_dirs)) {
        int j_le = wrap_dj(j-1);
        int j_lo = j/2;
        int i_l = ((i+j) % 2 == 0) ? i : wrap_i(i-1);
        assert(bond_r(bonds.pet_even[t][i][j_le]));
        if (i_l == i) {
          assert(bond_u(bonds.pet_odd[t][i_l][j_lo]));
        }
        else {
          assert(bond_d(bonds.pet_odd[t][i_l][j_lo]));
        }
        if (seen.pet_even[t][i][j_le] < gen) {
          seen.pet_even[t][i][j_le] = gen;
          sheet.pet_even[t][i][j_le] = prev_sheet;
          q_push((coord_t){t, i, j_le, EVEN});
        }
        else if (sheet.pet_even[t][i][j_le] != prev_sheet) {
          *temporal_wrap = true;
        }
        if (seen.pet_odd[t][i_l][j_lo] < gen) {
          seen.pet_odd[t][i_l][j_lo] = gen;
          sheet.pet_odd[t][i_l][j_lo] = prev_sheet;
          q_push((coord_t){t, i_l, j_lo, ODD});
        }
        else if (sheet.pet_odd[t][i_l][j_lo] != prev_sheet) {
          *temporal_wrap = true;
        }
      }
      // RIGHT (3-way)
      if (bond_r(bond_dirs)) {
        int j_re = wrap_dj(j+1);
        int j_ro = wrap_dj(j+1)/2;
        int i_r = ((i+j) % 2 == 0) ? wrap_i(i-1) : i;
        assert(bond_l(bonds.pet_even[t][i][j_re]));
        if (i_r == i) {
          assert(bond_u(bonds.pet_odd[t][i_r][j_ro]));
        }
        else {
          assert(bond_d(bonds.pet_odd[t][i_r][j_ro]));
        }
        if (seen.pet_even[t][i][j_re] < gen) {
          seen.pet_even[t][i][j_re] = gen;
          sheet.pet_even[t][i][j_re] = prev_sheet;
          q_push((coord_t){t, i, j_re, EVEN});
        }
        else if (sheet.pet_even[t][i][j_re] != prev_sheet) {
          *temporal_wrap = true;
        }
        if (seen.pet_odd[t][i_r][j_ro] < gen) {
          seen.pet_odd[t][i_r][j_ro] = gen;
          sheet.pet_odd[t][i_r][j_ro] = prev_sheet;
          q_push((coord_t){t, i_r, j_ro, ODD});
        }
        else if (sheet.pet_odd[t][i_r][j_ro] != prev_sheet) {
          *temporal_wrap = true;
        }
      }
      assert(!bond_u(bond_dirs));
      assert(!bond_d(bond_dirs));
    }
    // odd petals
    else {
      assert(seen.pet_odd[t][i][j] == gen);
      int prev_sheet = sheet.pet_odd[t][i][j];
      lattice[ind_pet_odd(t, i, j)] = spin;
      // printf("    (%d,%d,%d) sheet %d\n",
      //        ind_pet_odd(t, i, j)/(EROWS*ECOLS),
      //        (ind_pet_odd(t, i, j)%(EROWS*ECOLS))/ECOLS,
      //        ind_pet_odd(t, i, j)%(ECOLS),
      //        prev_sheet);
      bond_dirs_t bond_dirs = bonds.pet_odd[t][i][j];
      // TFWD
      if (bond_tfwd(bond_dirs)) {
        int tfwd = wrap_t(t+1);
        int next_sheet = (t+1 >= NT) ? prev_sheet+1 : prev_sheet;
        assert(bond_tbwd(bonds.pet_odd[tfwd][i][j]));
        if (seen.pet_odd[tfwd][i][j] < gen) {
          // printf("Not seen tfwd=%d\n", tfwd);
          seen.pet_odd[tfwd][i][j] = gen;
          sheet.pet_odd[tfwd][i][j] = next_sheet;
          q_push((coord_t){tfwd, i, j, ODD});
        }
        else if (sheet.pet_odd[tfwd][i][j] != next_sheet) {
          // printf("Seen, wrap tfwd=%d\n", tfwd);
          *temporal_wrap = true;
        }
        else {
          // printf("Seen, no wrap tfwd=%d\n", tfwd);
        }
      }
      // TBWD
      if (bond_tbwd(bond_dirs)) {
        int tbwd = wrap_t(t-1);
        int next_sheet = (t-1 < 0) ? prev_sheet-1 : prev_sheet;
        assert(bond_tfwd(bonds.pet_odd[tbwd][i][j]));
        if (seen.pet_odd[tbwd][i][j] < gen) {
          // printf("Not seen tbwd=%d\n", tbwd);
          seen.pet_odd[tbwd][i][j] = gen;
          sheet.pet_odd[tbwd][i][j] = next_sheet;
          q_push((coord_t){tbwd, i, j, ODD});
        }
        else if (sheet.pet_odd[tbwd][i][j] != next_sheet) {
          // printf("Seen, wrap tbwd=%d\n", tbwd);
          *temporal_wrap = true;
        }
        else {
          // printf("Seen, no wrap tbwd=%d\n", tbwd);
        }
      }
      // UP (3-way)
      if (bond_u(bond_dirs)) {
        bool even_row = i % 2 == 0;
        int j_r = wrap_dj(2*j + (even_row ? 0 : 1));
        int j_l = wrap_dj(j_r-1);
        assert(bond_r(bonds.pet_even[t][i][j_l]));
        assert(bond_l(bonds.pet_even[t][i][j_r]));
        if (seen.pet_even[t][i][j_l] < gen) {
          seen.pet_even[t][i][j_l] = gen;
          sheet.pet_even[t][i][j_l] = prev_sheet;
          q_push((coord_t){t, i, j_l, EVEN});
        }
        else if (sheet.pet_even[t][i][j_l] != prev_sheet) {
          *temporal_wrap = true;
        }
        if (seen.pet_even[t][i][j_r] < gen) {
          seen.pet_even[t][i][j_r] = gen;
          sheet.pet_even[t][i][j_r] = prev_sheet;
          q_push((coord_t){t, i, j_r, EVEN});
        }
        else if (sheet.pet_even[t][i][j_r] != prev_sheet) {
          *temporal_wrap = true;
        }
      }
      // DOWN (3-way)
      if (bond_d(bond_dirs)) {
        bool even_row = i % 2 == 0;
        int j_r = wrap_dj(2*j + (even_row ? 0 : 1));
        int j_l = wrap_dj(j_r-1);
        int i_d = wrap_i(i+1);
        assert(bond_r(bonds.pet_even[t][i_d][j_l]));
        assert(bond_l(bonds.pet_even[t][i_d][j_r]));
        if (seen.pet_even[t][i_d][j_l] < gen) {
          seen.pet_even[t][i_d][j_l] = gen;
          sheet.pet_even[t][i_d][j_l] = prev_sheet;
          q_push((coord_t){t, i_d, j_l, EVEN});
        }
        else if (sheet.pet_even[t][i_d][j_l] != prev_sheet) {
          *temporal_wrap = true;
        }
        if (seen.pet_even[t][i_d][j_r] < gen) {
          seen.pet_even[t][i_d][j_r] = gen;
          sheet.pet_even[t][i_d][j_r] = prev_sheet;
          q_push((coord_t){t, i_d, j_r, EVEN});
        }
        else if (sheet.pet_even[t][i_d][j_r] != prev_sheet) {
          *temporal_wrap = true;
        }
      }
      assert(!bond_l(bond_dirs));
      assert(!bond_r(bond_dirs));
    }
  }
}

double update_pet_spins() {
  memset(seen.pet_even, 0, sizeof(seen.pet_even));
  memset(seen.pet_odd, 0, sizeof(seen.pet_odd));
  memset(sheet.pet_even, 0, sizeof(sheet.pet_even));
  memset(sheet.pet_odd, 0, sizeof(sheet.pet_odd));
  bool temporal_wrap;
  int count;
  int acc = 0, rej = 0;
  // outer loop even, fixed
  // TODO: do we need to run this loop?
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (seen.pet_even[t][i][j]) {
          continue;
        }
        if (geom[ind_pet_even(0, i, j)] == FREE) {
          continue;
        }
        int gen = 1;
        seen.pet_even[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, EVEN});
        flood_fill_pet(geom[ind_pet_even(0, i, j)], gen, &temporal_wrap, &count);
        assert(q_empty());
        rej += count;
      }
    }
  }
  // outer loop odd, fixed
  // TODO: do we need to run this loop?
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (seen.pet_odd[t][i][j]) {
          continue;
        }
        if (geom[ind_pet_odd(0, i, j)] == FREE) {
          continue;
        }
        int gen = 1;
        seen.pet_odd[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, ODD});
        flood_fill_pet(geom[ind_pet_odd(0, i, j)], gen, &temporal_wrap, &count);
        assert(q_empty());
        rej += count;
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
        int gen = 1;
        seen.pet_even[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, EVEN});
        spin_t prev_val = lattice[ind_pet_even(t, i, j)];
        flood_fill_pet(rand_bool(), gen, &temporal_wrap, &count);
        // printf("New even cluster (%d)\n", count);
        assert(q_empty());
        if (temporal_wrap) { // unflip if wrapped
          // printf("temporal wrap, unflipping!\n");
          int prev_count = count;
          gen += 1;
          seen.pet_even[t][i][j] = gen;
          q_push((coord_t){t, i, j, EVEN});
          flood_fill_pet(prev_val, gen, &temporal_wrap, &count);
          assert(q_empty());
          assert(count == prev_count);
          rej += count;
        }
        else {
          acc += count;
        }
        // FORNOW:
        // measure_Gx(true);
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
        int gen = 1;
        seen.pet_odd[t][i][j] = gen;
        assert(q_empty());
        q_push((coord_t){t, i, j, ODD});
        spin_t prev_val = lattice[ind_pet_odd(t, i, j)];
        flood_fill_pet(rand_bool(), gen, &temporal_wrap, &count);
        // printf("New odd cluster (%d)\n", count);
        assert(q_empty());
        if (temporal_wrap) { // unflip if wrapped
          // printf("temporal wrap, unflipping!\n");
          int prev_count = count;
          gen += 1;
          seen.pet_odd[t][i][j] = gen;
          q_push((coord_t){t, i, j, ODD});
          flood_fill_pet(prev_val, gen, &temporal_wrap, &count);
          assert(q_empty());
          assert(count == prev_count);
          rej += count;
        }
        else {
          acc += count;
        }
        // FORNOW:
        // measure_Gx(true);
      }
    }
  }
  return acc / (double)(acc + rej);
}

void measure_HP_HE_pet_even(
    int t, int i, int j,
    uint64_t* HP_plus, uint64_t* HP_minus, uint64_t* HE_plus, uint64_t* HE_minus) {
  if (geom[ind_pet_even(0, i, j)] != FREE) {
    return;
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

  *HP_plus += delta_b && delta_a;
  *HP_minus += (!delta_b) && delta_a;
  *HE_plus += delta_a;
  *HE_minus += !delta_a;
}

void measure_HP_HE_pet_odd(
    int t, int i, int j,
    uint64_t* HP_plus, uint64_t* HP_minus, uint64_t* HE_plus, uint64_t* HE_minus) {
  if (geom[ind_pet_odd(0, i, j)] != FREE) {
    return;
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

  *HP_plus += delta_b && delta_a;
  *HP_minus += (!delta_b) && delta_a;
  *HE_plus += delta_a;
  *HE_minus += !delta_a;
}

void measure_HP_HE(uint64_t* HP_plus, uint64_t* HP_minus, uint64_t* HE_plus, uint64_t* HE_minus) {
  uint64_t tot_HP_plus = 0;
  uint64_t tot_HP_minus = 0;
  uint64_t tot_HE_plus = 0;
  uint64_t tot_HE_minus = 0;
  // even petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        measure_HP_HE_pet_even(t, i, j, &tot_HP_plus, &tot_HP_minus, &tot_HE_plus, &tot_HE_minus);
      }
    }
  }
  // odd petals
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        measure_HP_HE_pet_odd(t, i, j, &tot_HP_plus, &tot_HP_minus, &tot_HE_plus, &tot_HE_minus);
      }
    }
  }
  *HP_plus = tot_HP_plus;
  *HP_minus = tot_HP_minus;
  *HE_plus = tot_HE_plus;
  *HE_minus = tot_HE_minus;
}

void measure_HT_tri(int t, int i, int j, uint64_t* HT_plus, uint64_t* HT_minus) {
  if (geom[ind_tri(0, i, j)] != FREE) {
    return;
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
  *HT_plus += delta_a && delta_b;
  *HT_minus += (!delta_a) && delta_b;
}

void measure_HT(uint64_t* HT_plus, uint64_t* HT_minus) {
  uint64_t tot_HT_plus = 0;
  uint64_t tot_HT_minus = 0;
  for (int t = 0; t < NT; ++t) {
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        measure_HT_tri(t, i, j, &tot_HT_plus, &tot_HT_minus);
      }
    }
  }
  *HT_plus = tot_HT_plus;
  *HT_minus = tot_HT_minus;
}

void measure_MT_MP(uint64_t* MT, uint64_t* MP) {
  uint64_t tot_MT = 0;
  uint64_t tot_MP = 0;
  for (int t = 0; t < NT; ++t) {
    // triangles
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (geom[ind_tri(0, i, j)] != FREE) {
          continue;
        }
        tot_MT += lattice[ind_tri(t, i, j)];
      }
    }
    // even petals
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < 2*NCOLS; ++j) {
        if (geom[ind_pet_even(0, i, j)] != FREE) {
          continue;
        }
        tot_MP += lattice[ind_pet_even(t, i, j)];
      }
    }
    // odd petals
    for (int i = 0; i < NROWS; ++i) {
      for (int j = 0; j < NCOLS; ++j) {
        if (geom[ind_pet_odd(0, i, j)] != FREE) {
          continue;
        }
        tot_MP += lattice[ind_pet_odd(t, i, j)];
      }
    }
  }
  *MT = tot_MT;
  *MP = tot_MP;
}

void measure_Ex_pet_even(int t, int i, int j, int64_t* Ex) {
  if (geom[ind_pet_even(0, i, j)] != FREE) {
    return;
  }
  // relevant petals
  bool pet = lattice[ind_pet_even(t, i, j)];
  int j_l = j;
  int j_r = wrap_dj(j+1);

  // relevant triangles
  bool tri_l = lattice[ind_tri(t, i, j_l)];
  bool tri_r = lattice[ind_tri(t, i, j_r)];

  int diff_l = wrap_mod4(2*tri_l - 1 - 2*pet);
  int diff_r = wrap_mod4(2*pet - 2*tri_r + 1);
  assert((diff_l + diff_r) % 2 == 0);
  *Ex += (diff_l + diff_r)/2;
}

void measure_Ex_pet_odd(int t, int i, int j, int64_t* Ex) {
  // relevant petals
  bool pet = lattice[ind_pet_odd(t, i, j)];
  int i_u = i;
  int i_d = wrap_i(i+1);
  int j_even = 2*j + (i % 2);

  // relevant triangles
  bool tri_u = lattice[ind_tri(t, i_u, j_even)];
  bool tri_d = lattice[ind_tri(t, i_d, j_even)];

  int diff_u = wrap_mod4(2*tri_u - 1 - 2*pet);
  int diff_d = wrap_mod4(2*pet - 2*tri_d + 1);
  assert((diff_u + diff_d) % 2 == 0);
  *Ex += (diff_u + diff_d)/2;
}

void write_lattice(FILE *f) {
  assert(f != NULL);
  fwrite(lattice, 1, sizeof(lattice), f);
}

// const size_t STRLEN = 256;
#define STRLEN 255

typedef struct {
  double dt, KP, KE;
  int n_iter;
  int save_freq;
  int meas_freq;
  int string_sep;
  unsigned seed;
  const char* prefix;
  // derived
  int n_meas;
  char fname_ens[STRLEN+1];
  char fname_meta[STRLEN+1];
  char fname_HT[STRLEN+1];
  char fname_HP[STRLEN+1];
  char fname_HE[STRLEN+1];
  char fname_MT[STRLEN+1];
  char fname_MP[STRLEN+1];
  char fname_Mx[STRLEN+1];
  char fname_Hx[STRLEN+1];
  char fname_Ex[STRLEN+1];
} config_t;

void usage(const char* prog) {
  printf("Usage: %s -t <dt> -p <KP> -e <KE> -f <out_prefix> "
         "-i <n_iter> -s <save_freq> -m <meas_freq> "
         "[-b (obc|pbc|rhomb|rhomb_str|file)] [-c (cold|hot|cold_str|file)] "
         "[-y geom_file] [-z init_file] [-d dirac_str_file] "
         "[-x string_sep] [-r <seed>]\n", prog);
}

int parse_args(int argc, char** argv, config_t* cfg) {
  if (argc < 2) {
    usage(argv[0]);
    return E_ARGS;
  }

  bool set_dt = false, set_KP = false, set_KE = false,
      set_prefix = false, set_n_iter = false, set_save_freq = false,
      set_meas_freq = false, set_string_sep = false;
  double dt, KP, KE;
  bc_t bc_kind = PBC;
  init_t init_kind = HOT;
  const char* geom_file = NULL;
  const char* init_file = NULL;
  const char* dirac_str_file = NULL;
  const char* prefix;
  char c;
  while ((c = getopt(argc, argv, "t:p:e:f:i:s:m:b:c:r:x:y:z:d:")) != -1) {
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
      else if (strcmp(optarg, "rhomb") == 0) {
        bc_kind = RHOMB;
      }
      else if (strcmp(optarg, "rhomb_str") == 0) {
        bc_kind = RHOMB_STRING;
      }
      else if (strcmp(optarg, "rhomb2_str") == 0) {
        bc_kind = RHOMB2_STRING;
      }
      else if (strcmp(optarg, "file") == 0) {
        bc_kind = GEOM_FILE;
      }
      else {
        usage(argv[0]);
        printf("Invalid -b argument: %s\n", optarg);
        return E_ARGS;
      }
    }
    else if (c == 'x') {
      cfg->string_sep = atoi(optarg);
      set_string_sep = true;
    }
    else if (c == 'c') {
      if (strcmp(optarg, "cold") == 0) {
        init_kind = COLD;
      }
      else if (strcmp(optarg, "cold_str") == 0) {
        init_kind = COLD_STRING;
      }
      else if (strcmp(optarg, "hot") == 0) {
        init_kind = HOT;
      }
      else if (strcmp(optarg, "file") == 0) {
        init_kind = INIT_FILE;
      }
      else {
        usage(argv[0]);
        printf("Invalid -c argument: %s\n", optarg);
        return E_ARGS;
      }
    }
    else if (c == 'y') {
      geom_file = optarg;
    }
    else if (c == 'z') {
      init_file = optarg;
    }
    else if (c == 'd') {
      dirac_str_file = optarg;
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
  if (init_kind == COLD_STRING && !set_string_sep) {
    printf("Must provide string_sep for cold string init\n");
    return E_ARGS;
  }
  if (init_kind == INIT_FILE && !init_file) {
    printf("Must provide init_file for file init\n");
  }
  if (bc_kind == GEOM_FILE && !geom_file) {
    printf("Must provide geom_file for file geom\n");
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
  strncpy(cfg->fname_MT, prefix, STRLEN);
  strncpy(cfg->fname_MP, prefix, STRLEN);
  strncpy(cfg->fname_Mx, prefix, STRLEN);
  strncpy(cfg->fname_Hx, prefix, STRLEN);
  strncpy(cfg->fname_Ex, prefix, STRLEN);
  strncpy(cfg->fname_ens + len, ".ens.dat", STRLEN-len);
  strncpy(cfg->fname_meta + len, ".meta.dat", STRLEN-len);
  strncpy(cfg->fname_HT + len, ".HT.dat", STRLEN-len);
  strncpy(cfg->fname_HP + len, ".HP.dat", STRLEN-len);
  strncpy(cfg->fname_HE + len, ".HE.dat", STRLEN-len);
  strncpy(cfg->fname_MT + len, ".MT.dat", STRLEN-len);
  strncpy(cfg->fname_MP + len, ".MP.dat", STRLEN-len);
  strncpy(cfg->fname_Mx + len, ".Mx.dat", STRLEN-len);
  strncpy(cfg->fname_Hx + len, ".Hx.dat", STRLEN-len);
  strncpy(cfg->fname_Ex + len, ".Ex.dat", STRLEN-len);

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
  ret = init_dirac_str(dirac_str_file);
  assert(ret == E_OK);
  ret = init_geom(bc_kind, geom_file);
  assert(ret == E_OK);
  print_geom();
  ret = init_lat(init_kind, cfg->string_sep, init_file);
  assert(ret == E_OK);
  init_lat_apply_geom();
  measure_Gx(false);

  return ret;
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
  FILE *f_MT = fopen(cfg.fname_MT, "wb");
  FILE *f_MP = fopen(cfg.fname_MP, "wb");
  FILE *f_Mx = fopen(cfg.fname_Mx, "wb");
  FILE *f_Hx = fopen(cfg.fname_Hx, "wb");
  FILE *f_Ex = fopen(cfg.fname_Ex, "wb");
  if (f == NULL || f_meta == NULL ||
      f_HT == NULL || f_HP == NULL || f_HE == NULL ||
      f_MT == NULL || f_MP == NULL ||
      f_Mx == NULL || f_Hx == NULL || f_Ex == NULL)  {
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
  uint64_t* HT = malloc(len_Hi * sizeof(uint64_t));
  uint64_t* HP = malloc(len_Hi * sizeof(uint64_t));
  uint64_t* HE = malloc(len_Hi * sizeof(uint64_t));
  // Raw magnetization
  size_t len_Mi = cfg.n_meas;
  uint64_t* MT = malloc(len_Mi * sizeof(uint64_t));
  uint64_t* MP = malloc(len_Mi * sizeof(uint64_t));

  clock_t start = clock();
  double acc_tri = 0.0, acc_pet = 0.0;
  for (int i = 0; i < cfg.n_iter; ++i) {
    // progress
    if ((i+1) % 100 == 0) {
      clock_t ticks = clock() - start;
      double elapsed = ticks / (double)CLOCKS_PER_SEC;
      double expected = (elapsed * cfg.n_iter) / (i+1);
      double rate = (i+1) / elapsed;
      printf(
          "Iter %d / %d (%.2f / %.2fs | %.2f it/s)\n",
          i+1, cfg.n_iter, elapsed, expected, rate);
      printf(
          "... acc tri = %.2f, pet = %.2f\n",
          100.0*acc_tri / (i+1), 100.0*acc_pet / (i+1));
    }
    // FORNOW:
    // printf("Iter %d (pre), measure Gx:\n", i+1);
    // measure_Gx(true);
    // triangle sublattice
    sample_tri_bonds();
    acc_tri += update_tri_spins();
    // FORNOW:
    // printf("Iter %d (tri), measure Gx:\n", i+1);
    // measure_Gx(true);
    // petal sublattice
    sample_pet_bonds();
    acc_pet += update_pet_spins();
    // FORNOW:
    // printf("Iter %d (pet), measure Gx:\n", i+1);
    // measure_Gx(true);

    // accumulate Mx
    for (int t = 0; t < NT; ++t) {
      for (int i = 0; i < EROWS * ECOLS; ++i) {
        mag_accum[i] += lattice[t * EROWS * ECOLS + i];
      }
      mag_n++;
    }
    // accumulate Hx, Ex
    for (int t = 0; t < NT; ++t) {
      // triangles (HT)
      for (int i = 0; i < NROWS; ++i) {
        for (int j = 0; j < 2*NCOLS; ++j) {
          int ind = ind_tri(0, i, j);
          measure_HT_tri(
              t, i, j, &HT_plus_accum[ind], &HT_minus_accum[ind]);
        }
      }
      // even petals (HP, HE)
      for (int i = 0; i < NROWS; ++i) {
        for (int j = 0; j < 2*NCOLS; ++j) {
          int ind = ind_pet_even(0, i, j);
          measure_HP_HE_pet_even(
              t, i, j, &HP_plus_accum[ind], &HP_minus_accum[ind],
              &HE_plus_accum[ind], &HE_minus_accum[ind]);
          measure_Ex_pet_even(t, i, j, &Ex_accum[ind]);
        }
      }
      // odd petals (HP, HE)
      for (int i = 0; i < NROWS; ++i) {
        for (int j = 0; j < NCOLS; ++j) {
          int ind = ind_pet_odd(0, i, j);
          measure_HP_HE_pet_odd(
              t, i, j, &HP_plus_accum[ind], &HP_minus_accum[ind],
              &HE_plus_accum[ind], &HE_minus_accum[ind]);
          measure_Ex_pet_odd(t, i, j, &Ex_accum[ind]);
        }
      }
      H_n++;
      Ex_n++;
    }
    // measure
    if ((i+1) % cfg.meas_freq == 0) {
      const int ind = ((i+1) / cfg.meas_freq) - 1;
      measure_HT(&HT[2*ind], &HT[2*ind+1]);
      measure_HP_HE(&HP[2*ind], &HP[2*ind+1], &HE[2*ind], &HE[2*ind+1]);
      measure_MT_MP(&MT[ind], &MP[ind]);
    }
    // save
    if ((i+1) % cfg.save_freq == 0) {
      write_lattice(f);
    }
  }

  // FORNOW:
  printf("Final Gx state:\n");
  measure_Gx(true);

  fwrite(HT, sizeof(uint64_t), len_Hi, f_HT);
  fwrite(HP, sizeof(uint64_t), len_Hi, f_HP);
  fwrite(HE, sizeof(uint64_t), len_Hi, f_HE);
  fwrite(MT, sizeof(uint64_t), len_Mi, f_MT);
  fwrite(MP, sizeof(uint64_t), len_Mi, f_MP);

  double mag[EROWS * ECOLS];
  for (int i = 0; i < EROWS * ECOLS; ++i) {
    mag[i] = mag_accum[i] / (double)mag_n;
  }
  fwrite(mag, sizeof(double), EROWS*ECOLS, f_Mx);
  double Ex[EROWS * ECOLS];
  for (int i = 0; i < EROWS * ECOLS; ++i) {
    Ex[i] = Ex_accum[i] / (double)Ex_n;
  }
  fwrite(Ex, sizeof(double), EROWS*ECOLS, f_Ex);

  double HT_plus[EROWS * ECOLS], HT_minus[EROWS * ECOLS];
  double HP_plus[EROWS * ECOLS], HP_minus[EROWS * ECOLS];
  double HE_plus[EROWS * ECOLS], HE_minus[EROWS * ECOLS];
  for (int i = 0; i < EROWS * ECOLS; ++i) {
    HT_plus[i] = HT_plus_accum[i] / (double)H_n;
    HT_minus[i] = HT_minus_accum[i] / (double)H_n;
    HP_plus[i] = HP_plus_accum[i] / (double)H_n;
    HP_minus[i] = HP_minus_accum[i] / (double)H_n;
    HE_plus[i] = HE_plus_accum[i] / (double)H_n;
    HE_minus[i] = HE_minus_accum[i] / (double)H_n;
  }
  fwrite(HT_plus, sizeof(double), EROWS*ECOLS, f_Hx);
  fwrite(HT_minus, sizeof(double), EROWS*ECOLS, f_Hx);
  fwrite(HP_plus, sizeof(double), EROWS*ECOLS, f_Hx);
  fwrite(HP_minus, sizeof(double), EROWS*ECOLS, f_Hx);
  fwrite(HE_plus, sizeof(double), EROWS*ECOLS, f_Hx);
  fwrite(HE_minus, sizeof(double), EROWS*ECOLS, f_Hx);

  free(HT);
  free(HP);
  free(HE);
  free(MT);
  free(MP);
  fclose(f);
  fclose(f_HT);
  fclose(f_HP);
  fclose(f_HE);
  fclose(f_Mx);
  fclose(f_Hx);
}
