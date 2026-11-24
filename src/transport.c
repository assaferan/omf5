/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : transport.c
 *
 * Description: Hecke columns by carrying the classes and
 *              the isometries of the p-neighbors along
 *              isometries in Z[1/q] (issue #63).
 *
 * Conventions: an isometry s from q1 to q2 is S / d with
 * S^t q1 S = d^2 q2 (isometry_is_isom). As a map x -> S x / d on
 * column vectors it takes the lattice of q2 onto a lattice in the
 * space of q1. The line of a column vector v of L_i is moved by an
 * isometry g of L_i to that of g v.
 *
 * Let L = L_i, and c_v : L_j -> N_v an isometry onto the p-neighbor
 * of L at v. Let g be an isometry of the space of L, with g(L) the
 * q-neighbor of L at some line, and l_g the line of L at which
 * g^(-1)(L) is its q-neighbor. The p-neighbor of L at g v is g(N'),
 * where N' is N_v except at q, where it is g^(-1)(L). Then c_v^(-1)(N')
 * is the q-neighbor of L_j at the line w of c_v^(-1)(l_g), and with
 * t : L_{j'} -> c_v^(-1)(N') from the table of T_q at (j, w),
 * c_{gv} = g c_v t is an isometry from L_{j'} onto the p-neighbor at g v.
 * For an automorphism a of L, c_{av} = a c_v.
 *
 **********************************************************
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <flint/ulong_extras.h>

#include "arith.h"
#include "aut_grp.h"
#include "hash.h"
#include "isometry.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "spinor.h"
#include "transport.h"

// ---- a set of lines mod p, by open addressing on the index of the normalized vector ----

typedef struct {
  W64* keys;
  W64 mask;
  W64 count;
} line_set_struct;

#define LINE_SET_EMPTY ((W64)(-1))

static void line_set_init(line_set_struct* S, W64 min_size)
{
  W64 cap = 16;
  W64 i;
  while (cap < 2 * min_size)
    cap <<= 1;
  S->keys = (W64*)malloc(cap * sizeof(W64));
  for (i = 0; i < cap; i++)
    S->keys[i] = LINE_SET_EMPTY;
  S->mask = cap - 1;
  S->count = 0;
}

static void line_set_clear(line_set_struct* S)
{
  free(S->keys);
}

static W64 line_set_slot(const line_set_struct* S, W64 key)
{
  W64 h = (key * 0x9E3779B97F4A7C15ULL) & S->mask;
  while ((S->keys[h] != LINE_SET_EMPTY) && (S->keys[h] != key))
    h = (h + 1) & S->mask;
  return h;
}

// inserts the key, and returns whether it was not there. The set is kept at most half full.
static bool line_set_insert(line_set_struct* S, W64 key)
{
  W64 h, i, old_cap;
  W64* old;

  h = line_set_slot(S, key);
  if (S->keys[h] == key)
    return false;
  S->keys[h] = key;
  S->count++;
  if (2 * S->count > S->mask) {
    old = S->keys;
    old_cap = S->mask + 1;
    S->mask = 2 * old_cap - 1;
    S->keys = (W64*)malloc(2 * old_cap * sizeof(W64));
    for (i = 0; i <= S->mask; i++)
      S->keys[i] = LINE_SET_EMPTY;
    for (i = 0; i < old_cap; i++)
      if (old[i] != LINE_SET_EMPTY)
	S->keys[line_set_slot(S, old[i])] = old[i];
    free(old);
  }
  return true;
}

// the normalized vector of the line of v mod p, and its index in base p
static W64 line_key(vector_t v, Z64 p)
{
  W64 key;
  int k;

  vector_mod_p(v, p);
  normalize_mod_p(v, p);
  key = 0;
  for (k = QF_RANK - 1; k >= 0; k--)
    key = key * (W64)p + (W64)v[k];
  return key;
}

// ---- the set of the isotropic lines mod p reached, as a bitmap of 2 (p^3 + p^2 + p + 1) bits ----
//
// For a coordinate k with Q(e_k) = a != 0 mod p, Q(v) = a v_k^2 + b v_k + c, where b = sum_{i != k} G_ik v_i
// and c do not depend on v_k. An isotropic v is not a multiple of e_k, so its other four coordinates are a
// point y of P^3, and v_k is one of the two roots of this quadratic, whose sum is -b/a. The line of v is
// then the index of y among the p^3 + p^2 + p + 1 points of P^3, and one bit telling which root v_k is.
// When there is no such k (Q(e_k) = 0 mod p for each k), the lines are kept in the hash set above.

typedef struct {
  int k;               // the coordinate solved for, or -1 for the hash set
  Z64 p;
  Z64 a_inv;           // the inverse of Q(e_k) mod p
  Z64 g[QF_RANK];      // G_ik mod p
  W64 offs[QF_RANK - 1]; // the index of the first point of P^3 whose first nonzero coordinate is the t-th
  W64* bits;
  line_set_struct hash;
  Z64 g_full[QF_RANK][QF_RANK]; // G mod 2p, for the check
  bool check;           // with OMF5_TRANSPORT_CHECK, the lines are also kept in shadow, to check the bitmap
  line_set_struct shadow;
} seen_set_struct;

// the k for a Gram matrix mod p, or -1
static int seen_set_chart(const square_matrix_t G, Z64 p)
{
  int k;
  for (k = QF_RANK - 1; k >= 0; k--)
    if (((G[k][k] / 2) % p) != 0)
      return k;
  return -1;
}

static W64 seen_set_num_points(Z64 p)
{
  return (W64)p * p * p + (W64)p * p + p + 1;
}

static void seen_set_init(seen_set_struct* S, const square_matrix_t G, Z64 p, int k)
{
  Z64 x, y, a;
  W64 num_words;
  int i, t;

  S->p = p;
  S->k = k;
  S->bits = NULL;
  S->check = (k >= 0) && (getenv("OMF5_TRANSPORT_CHECK") != NULL);
  if (S->check) {
    int r, c;
    line_set_init(&(S->shadow), seen_set_num_points(p));
    for (r = 0; r < QF_RANK; r++)
      for (c = 0; c < QF_RANK; c++)
	S->g_full[r][c] = ((G[r][c] % (2 * p)) + 2 * p) % (2 * p);
  }
  if (k < 0) {
    line_set_init(&(S->hash), seen_set_num_points(p));
    return;
  }
  a = ((G[k][k] / 2) % p + p) % p;
  gcdext(a, p, &x, &y);
  S->a_inv = ((x % p) + p) % p;
  for (i = 0; i < QF_RANK; i++)
    S->g[i] = ((G[i][k] % p) + p) % p;
  S->offs[0] = 0;
  for (t = 1; t < QF_RANK - 1; t++) {
    W64 n = 1;
    for (i = t; i < QF_RANK - 1; i++)
      n *= (W64)p;
    S->offs[t] = S->offs[t - 1] + n;
  }
  num_words = (2 * seen_set_num_points(p) + 63) / 64;
  S->bits = (W64*)calloc(num_words, sizeof(W64));
  if (S->bits == NULL) {
    fprintf(stderr, "Error: out of memory for the lines of T_%ld in hecke_col_all_conductors_transport\n", (long)p);
    exit(EXIT_FAILURE);
  }
}

static void seen_set_clear(seen_set_struct* S)
{
  if (S->k < 0)
    line_set_clear(&(S->hash));
  else
    free(S->bits);
  if (S->check)
    line_set_clear(&(S->shadow));
}


// the bit of the line of v, for an isotropic v
static W64 seen_set_bit(const seen_set_struct* S, const vector_t v)
{
  Z64 p = S->p, y[QF_RANK - 1], x, inv, s, b, other;
  int i, n, t, k = S->k;
  W64 idx;

  for (i = 0, n = 0; i < QF_RANK; i++)
    if (i != k)
      y[n++] = ((v[i] % p) + p) % p;
  x = ((v[k] % p) + p) % p;
  for (t = 0; (t < QF_RANK - 1) && (y[t] == 0); t++);
  assert(t < QF_RANK - 1);
  gcdext(y[t], p, &inv, &s);
  inv = ((inv % p) + p) % p;
  b = 0;
  for (i = 0, n = 0; i < QF_RANK; i++)
    if (i != k) {
      y[n] = (y[n] * inv) % p;
      b = (b + S->g[i] * y[n]) % p;
      n++;
    }
  x = (x * inv) % p;
  // the other root: -b/a - x
  other = ((((p - b) % p) * S->a_inv) % p - x + p) % p;
  idx = 0;
  for (i = t + 1; i < QF_RANK - 1; i++)
    idx = idx * (W64)p + (W64)y[i];
  return 2 * (S->offs[t] + idx) + ((x > other) ? 1 : 0);
}

// inserts the line of v, and returns whether it was not there
static bool seen_set_insert(seen_set_struct* S, const vector_t v)
{
  W64 bit, m;
  vector_t w;

  if (S->k < 0) {
    vector_set(w, v);
    return line_set_insert(&(S->hash), line_key(w, S->p));
  }
  bit = seen_set_bit(S, v);
  if (S->check) {
    Z64 qv = 0;
    int r, c;
    vector_set(w, v);
    for (r = 0; r < QF_RANK; r++)
      for (c = 0; c < QF_RANK; c++)
	qv = (qv + (((S->g_full[r][c] * ((v[r] % S->p + S->p) % S->p)) % (2 * S->p)) * ((v[c] % S->p + S->p) % S->p))) % (2 * S->p);
    if (qv % (2 * S->p) != 0)
      fprintf(stderr, "transport: a line of the bitmap is not isotropic mod %ld\n", (long)S->p);
    if (line_set_insert(&(S->shadow), line_key(w, S->p)) != !((S->bits[bit >> 6] >> (bit & 63)) & 1))
      fprintf(stderr, "transport: the bitmap of the lines differs from the hash set, p = %ld\n", (long)S->p);
  }
  m = ((W64)1) << (bit & 63);
  if (S->bits[bit >> 6] & m)
    return false;
  S->bits[bit >> 6] |= m;
  return true;
}

// res = S v mod p, for an integral matrix S (the numerator of an isometry, whose denominator is prime to p)
static void mat_vec_mod_p(vector_t res, const square_matrix_t S, const vector_t v, Z64 p)
{
  int r, k;
  Z64 acc;

  for (r = 0; r < QF_RANK; r++) {
    acc = 0;
    for (k = 0; k < QF_RANK; k++)
      acc = (acc + ((S[r][k] % p + p) % p) * (v[k] % p)) % p;
    res[r] = acc;
  }
}

// divides the numerators and denominators of s and of its inverse by their gcd
static Z64 z64_gcd(Z64 a, Z64 b)
{
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  while (b != 0) {
    Z64 t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static void isometry_reduce(isometry_t s)
{
  Z64 g;
  int r, k;

  g = s->denom;
  for (r = 0; r < QF_RANK; r++)
    for (k = 0; k < QF_RANK; k++)
      g = z64_gcd(g, s->s[r][k]);
  if (g < 0) g = -g;
  if (g > 1) {
    for (r = 0; r < QF_RANK; r++)
      for (k = 0; k < QF_RANK; k++)
	s->s[r][k] /= g;
    s->denom /= g;
  }
  g = s->inv_denom;
  for (r = 0; r < QF_RANK; r++)
    for (k = 0; k < QF_RANK; k++)
      g = z64_gcd(g, s->s_inv[r][k]);
  if (g < 0) g = -g;
  if (g > 1) {
    for (r = 0; r < QF_RANK; r++)
      for (k = 0; k < QF_RANK; k++)
	s->s_inv[r][k] /= g;
    s->inv_denom /= g;
  }
}

// the line mod q of S x / d, for an integral x, when it is integral at q: d = q^e m with m prime to q
static void line_mod_q(vector_t res, const square_matrix_t S, Z64 d, const vector_t x, Z64 q)
{
  __int128 acc[QF_RANK];
  int r, k, e;
  Z64 m;

  for (r = 0; r < QF_RANK; r++) {
    acc[r] = 0;
    for (k = 0; k < QF_RANK; k++)
      acc[r] += (__int128)S[r][k] * x[k];
  }
  if (d < 0) d = -d;
  for (e = 0, m = d; m % q == 0; e++, m /= q)
    for (r = 0; r < QF_RANK; r++) {
      assert(acc[r] % q == 0);
      acc[r] /= q;
    }
  (void)e;
  for (r = 0; r < QF_RANK; r++)
    res[r] = (Z64)(((acc[r] % q) + q) % q);
}

// the line l_g of L at which g^(-1)(L) is its q-neighbor: q g^(-1) is integral, of rank 1 mod q
static void line_of_inverse(vector_t l, const isometry_t g, Z64 q)
{
  int r, k;
  Z64 d = g->inv_denom;
  bool found = false;

  for (k = 0; (k < QF_RANK) && !found; k++) {
    for (r = 0; r < QF_RANK; r++) {
      __int128 x = (__int128)g->s_inv[r][k] * q;
      assert(x % d == 0);
      x /= d;
      l[r] = (Z64)(((x % q) + q) % q);
      if (l[r] != 0)
	found = true;
    }
  }
  assert(found);
}

// ---- the p^2-neighbors (T_{p^2}), by their lattices ----
//
// A p^2-neighbor N of L is determined by N_p, and p N contains p^3 L (as in nbr_data_build_neighbor), so
// it is determined by the Z-span of the columns of a matrix B with p N = B Z^5 + p^3 L, and the key of N
// is the Hermite normal form of this span, with its entries in [0, p^3). An isometry g of L[1/q] which is
// integral at p maps it to the neighbor at g B (its denominator, a power of q, is a unit at p).

// the Hermite normal form H (rows, upper triangular, H[r][c] in [0, H[c][c]) for r < c) of the Z-span
// of the columns of B and of m Z^5
static void lattice_hnf_mod(Z64 H[QF_RANK][QF_RANK], const Z64 B[QF_RANK][QF_RANK], Z64 m)
{
  Z64 A[2 * QF_RANK][QF_RANK];
  int r, c, i, j;

  for (r = 0; r < QF_RANK; r++)
    for (c = 0; c < QF_RANK; c++) {
      A[r][c] = ((B[c][r] % m) + m) % m;
      A[QF_RANK + r][c] = (r == c) ? m : 0;
    }
  for (j = 0; j < QF_RANK; j++) {
    // gcd of the column j of the rows j, j+1, ... into the row j
    for (i = j + 1; i < 2 * QF_RANK; i++) {
      Z64 a = A[j][j], b = A[i][j], g, x, y;
      if (b == 0)
	continue;
      if (a == 0) {
	for (c = 0; c < QF_RANK; c++) {
	  Z64 t = A[j][c];
	  A[j][c] = A[i][c];
	  A[i][c] = t;
	}
	continue;
      }
      g = gcdext(a, b, &x, &y);
      if (g < 0) { g = -g; x = -x; y = -y; }
      for (c = 0; c < QF_RANK; c++) {
	__int128 u = (__int128)x * A[j][c] + (__int128)y * A[i][c];
	__int128 v = (__int128)(a / g) * A[i][c] - (__int128)(b / g) * A[j][c];
	A[j][c] = (Z64)(((u % m) + m) % m);
	A[i][c] = (Z64)(((v % m) + m) % m);
      }
      // the pivot is g exactly (g divides m, as m is a row of the span)
      A[j][j] = g;
      A[i][j] = 0;
    }
    if (A[j][j] == 0)
      A[j][j] = m;
    // reduce the rows above
    for (r = 0; r < j; r++) {
      Z64 t = A[r][j] / A[j][j];
      if (t != 0)
	for (c = j; c < QF_RANK; c++)
	  A[r][c] = (((A[r][c] - t * A[j][c]) % m) + m) % m;
      // the entries are in [0, m) and A[j][j] divides m, so A[r][j] is now in [0, A[j][j])
    }
  }
  for (r = 0; r < QF_RANK; r++)
    for (c = 0; c < QF_RANK; c++)
      H[r][c] = (c < r) ? 0 : A[r][c];
}

#define LATTICE_KEY_WORDS 8

typedef struct {
  W64* keys;     // LATTICE_KEY_WORDS * cap words
  bool* used;
  W64 mask;
  W64 count;
  int nw;        // the words of a key
  int bits;      // the bits of an entry
  Z64 m;
} lattice_set_struct;

static void lattice_set_init(lattice_set_struct* S, Z64 p, W64 min_size)
{
  W64 cap = 16, i;
  S->m = p * p * p;
  for (S->bits = 1; (((Z64)1) << S->bits) <= S->m; S->bits++);
  S->nw = (15 * S->bits + 63) / 64;
  assert(S->nw <= LATTICE_KEY_WORDS);
  while (cap < 2 * min_size)
    cap <<= 1;
  S->keys = (W64*)malloc(cap * S->nw * sizeof(W64));
  S->used = (bool*)malloc(cap * sizeof(bool));
  if ((S->keys == NULL) || (S->used == NULL)) {
    fprintf(stderr, "Error: out of memory for the neighbors of T_{%ld^2} in hecke_col_all_conductors_transport\n", (long)p);
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < cap; i++)
    S->used[i] = false;
  S->mask = cap - 1;
  S->count = 0;
}

static void lattice_set_clear(lattice_set_struct* S)
{
  free(S->keys);
  free(S->used);
}

// the key of the span of the columns of B and m Z^5: the 15 entries of its Hermite form, packed
static void lattice_key(W64 key[LATTICE_KEY_WORDS], const lattice_set_struct* S, const Z64 B[QF_RANK][QF_RANK])
{
  Z64 H[QF_RANK][QF_RANK];
  int r, c, w, pos;

  lattice_hnf_mod(H, B, S->m);
  for (w = 0; w < S->nw; w++)
    key[w] = 0;
  pos = 0;
  for (r = 0; r < QF_RANK; r++)
    for (c = r; c < QF_RANK; c++) {
      W64 x = (W64)H[r][c];
      int b;
      for (b = 0; b < S->bits; b++, pos++)
	if ((x >> b) & 1)
	  key[pos >> 6] |= ((W64)1) << (pos & 63);
    }
}

static W64 lattice_set_slot(const lattice_set_struct* S, const W64* key)
{
  W64 h = 0, i;
  int w;
  for (w = 0; w < S->nw; w++)
    h = (h ^ key[w]) * 0x9E3779B97F4A7C15ULL;
  // the product carries only upwards, so the high bits are folded into the low ones, which the mask keeps
  h ^= h >> 32;
  h &= S->mask;
  for (i = h; S->used[i]; i = (i + 1) & S->mask) {
    for (w = 0; (w < S->nw) && (S->keys[i * S->nw + w] == key[w]); w++);
    if (w == S->nw)
      return i;
  }
  return i;
}

// inserts the key, and returns whether it was not there. The set is kept at most half full.
static bool lattice_set_insert(lattice_set_struct* S, const W64* key)
{
  W64 h, i, old_cap;
  W64* old_keys;
  bool* old_used;
  int w;

  h = lattice_set_slot(S, key);
  if (S->used[h])
    return false;
  S->used[h] = true;
  for (w = 0; w < S->nw; w++)
    S->keys[h * S->nw + w] = key[w];
  S->count++;
  if (2 * S->count > S->mask) {
    old_keys = S->keys;
    old_used = S->used;
    old_cap = S->mask + 1;
    S->mask = 2 * old_cap - 1;
    S->keys = (W64*)malloc(2 * old_cap * S->nw * sizeof(W64));
    S->used = (bool*)malloc(2 * old_cap * sizeof(bool));
    if ((S->keys == NULL) || (S->used == NULL)) {
      fprintf(stderr, "Error: out of memory in hecke_col_all_conductors_transport\n");
      exit(EXIT_FAILURE);
    }
    for (i = 0; i <= S->mask; i++)
      S->used[i] = false;
    for (i = 0; i < old_cap; i++)
      if (old_used[i]) {
	h = lattice_set_slot(S, old_keys + i * S->nw);
	S->used[h] = true;
	for (w = 0; w < S->nw; w++)
	  S->keys[h * S->nw + w] = old_keys[i * S->nw + w];
      }
    free(old_keys);
    free(old_used);
  }
  return true;
}

// res = S B mod m, for an integral matrix S
static void mat_mat_mod(Z64 res[QF_RANK][QF_RANK], const square_matrix_t S, const Z64 B[QF_RANK][QF_RANK], Z64 m)
{
  int r, c, k;
  for (r = 0; r < QF_RANK; r++)
    for (c = 0; c < QF_RANK; c++) {
      __int128 acc = 0;
      for (k = 0; k < QF_RANK; k++)
	acc += (__int128)(((S[r][k] % m) + m) % m) * B[k][c];
      res[r][c] = (Z64)(acc % m);
    }
}

typedef struct {
  vector_t v;                  // the line, for T_p
  Z64 B[QF_RANK][QF_RANK];     // p times a basis of the neighbor mod p^3, for T_{p^2}
  int cls;
  W64 bits;
  isometry_t c;
} transport_item;

// The stack of lines to visit holds at most this many items (about 2 MB). A line reached when it is
// full is not marked as seen: it is reached again from another line, or found later as the first line of
// a new orbit, at the cost of one more isometry test. At D = 61, p = 97 this gives 48351 isometry tests
// for the 8 columns, with 10 MB, in less than the time of the orbits.
#define TRANSPORT_MAX_STACK (((slong)1) << 12)

static slong transport_memory_mb = TRANSPORT_DEFAULT_MEMORY_MB;

void hecke_transport_set_memory(slong mb)
{
  transport_memory_mb = mb;
}

// the memory of the lines of p, in bits: 2 a line with the bitmap, up to 256 with the hash set
static W64 transport_bits_needed(int p, bool bitmap)
{
  W64 lines = (W64)p * p * p + (W64)p * p + p + 1;
  return lines * (bitmap ? TRANSPORT_BITS_PER_LINE : TRANSPORT_BITS_PER_LINE_HASH);
}

bool hecke_transport_fits(int p)
{
  return transport_bits_needed(p, true) <= (((W64)transport_memory_mb) << 23);
}

int hecke_transport_max_prime(void)
{
  int p, best = 0;
  for (p = 2; hecke_transport_fits(p); p = (int)n_nextprime(p, 1))
    best = p;
  return best;
}

// the spinor bits of the isometry c : keys[i] -> keys[j], normalized to the frame of keys[0] as in
// process_isotropic_vector_nbr_data_all_conductors
static W64 transport_bits(const genus_t genus, int i, int j, const isometry_t c)
{
  isometry_t s_norm, s_inv;
  W64 bits;

  isometry_init(s_norm);
  isometry_mul(s_norm, genus->isoms[i], c);
  isometry_inv(s_inv, genus->isoms[j]);
  isometry_muleq_right(s_norm, s_inv);
  bits = spinor_norm_isom(genus->spinor, s_norm);
  isometry_clear(s_norm);
  isometry_clear(s_inv);
  return bits;
}

// the images of an item under the generator a, as in the comments of hecke_col_all_conductors_transport:
// c_{av} = a c_v for an automorphism, and c_{gv} = g c_v t for g, with t from the table at (cls, line of
// c_v^(-1) l_g mod q)
static void transport_child(transport_item* child, const transport_item* item, int a, int num_auts,
			    const isometry_struct* g, const vector_t* lgs, const W64* bits_gen, W64 mask,
			    const tq_table_t T)
{
  isometry_t tmp;
  slong idx;

  isometry_init(child->c);
  if (a < num_auts) {
    isometry_mul(child->c, g, item->c);
    child->cls = item->cls;
    child->bits = item->bits ^ bits_gen[a];
  }
  else {
    vector_t lw;
    line_mod_q(lw, item->c->s_inv, item->c->inv_denom, lgs[a - num_auts], T->q);
    idx = tq_table_pos(T, item->cls, lw);
    if (idx < 0) {
      fprintf(stderr, "Error: a carried line is not in the table of T_%ld\n", (long)T->q);
      exit(EXIT_FAILURE);
    }
    isometry_init(tmp);
    isometry_mul(tmp, g, item->c);
    isometry_mul(child->c, tmp, T->isom[item->cls][idx]);
    isometry_clear(tmp);
    child->cls = T->cls[item->cls][idx];
    child->bits = item->bits ^ bits_gen[a] ^ (T->spin[item->cls][idx] & mask);
  }
  isometry_reduce(child->c);
}

// with OMF5_TRANSPORT_CHECK: the carried isometry is an isometry onto the class, with denominators
// products of p and q
static void transport_check_child(transport_item* child, int a, int num_auts, int gen_idx, int p,
				  Z64 q, const genus_t genus)
{
  Z64 d1 = child->c->denom, d2 = child->c->inv_denom;

  if (!isometry_is_isom(child->c, genus->genus_reps->keys[gen_idx], genus->genus_reps->keys[child->cls]))
    fprintf(stderr, "transport: not an isometry, generator %d (%s), class %d -> %d\n", a,
	    (a < num_auts) ? "aut" : "g", gen_idx, child->cls);
  while (d1 % p == 0) d1 /= p;
  while (d1 % q == 0) d1 /= q;
  while (d2 % p == 0) d2 /= p;
  while (d2 % q == 0) d2 /= q;
  if ((d1 != 1 && d1 != -1) || (d2 != 1 && d2 != -1))
    fprintf(stderr, "transport: denominators %ld, %ld, generator %d (%s), class %d -> %d\n",
	    (long)child->c->denom, (long)child->c->inv_denom, a, (a < num_auts) ? "aut" : "g", gen_idx,
	    child->cls);
}

slong hecke_col_all_conductors_transport(slong* hist, int p, int k, int gen_idx, const genus_t genus,
					 const tq_table_t T)
{
  neighbor_manager_t nbr_man;
  nbr_data_t nd;
  aut_grp_t grp;
  square_matrix_t nbr;
  isometry_t *auts, *gs, s_nbr, hash_isom;
  W64 *bits_gen, mask;
  W64 key[LATTICE_KEY_WORDS];
  bool check, is_new, has_root;
  vector_t *lgs, u, w;
  Z64 B[QF_RANK][QF_RANK];
  slong num_gs, idx, lc, num_items, cap_items, num_roots, num_nbrs;
  int a, i, r, c, num_isom, chart;
  double theta_time, isom_time;
  Z64 q = T->q, m = (Z64)p * p * p;
  seen_set_struct seen;
  lattice_set_struct lseen;
  transport_item* stack;
  transport_item item, child;
  fmpz_mat_t nbr_fmpz, s_fmpz;

  assert(q != p);
  assert((k == 1) || (k == 2));
  assert(gen_idx < T->num_classes);
  chart = -1;
  if (k == 1) {
    // the lines are kept in a bitmap when the Gram matrix of the class has a diagonal entry prime to p,
    // and otherwise in a hash set, which takes more memory: then, if it does not fit, the column is left
    // to the orbits (-1)
    chart = seen_set_chart(genus->genus_reps->keys[gen_idx], p);
    if ((chart < 0) && (transport_bits_needed(p, false) > (((W64)transport_memory_mb) << 23)))
      return -1;
  }
  else {
    // the p^2-neighbors are kept in a hash set of their keys, of up to 4 slots of 8 LATTICE_KEY_WORDS
    // bytes and a flag for each of them; if it does not fit, the column is left to the neighbors (-1)
    nbr_data_init(nd, genus->genus_reps->keys[gen_idx], p, 2);
    num_nbrs = number_of_neighbors(nd, false);
    if ((W64)num_nbrs * 4 * (8 * LATTICE_KEY_WORDS + 1) > (((W64)transport_memory_mb) << 20)) {
      nbr_data_clear(nd);
      return -1;
    }
  }
  num_isom = 0;
  theta_time = isom_time = 0;

  // the generators: the automorphisms of L (grp->gens, as isometries), and the
  // isometries of the table which go back to L
  aut_grp_init_square_matrix(grp, genus->genus_reps->keys[gen_idx]);
  auts = (isometry_t*)malloc((grp->num_gens + 1) * sizeof(isometry_t));
  // the elements of nbr_process_init are the transposes of the products of grp->gens, and their
  // transposes are isometries, so the generators themselves are isometries
  for (a = 0; a < grp->num_gens; a++)
    isometry_init_set_square_matrix(auts[a], grp->gens[a], 1);
  check = (getenv("OMF5_TRANSPORT_CHECK") != NULL);
  mask = (((W64)1) << genus->spinor->num_primes) - 1;
  gs = (isometry_t*)malloc((T->num_lines[gen_idx] + 1) * sizeof(isometry_t));
  lgs = (vector_t*)malloc((T->num_lines[gen_idx] + 1) * sizeof(vector_t));
  num_gs = 0;
  // the spinor bits of each generator: the spinor norm is multiplicative, so the bits of a product of
  // isometries, normalized to the frame of keys[0], are the XOR of theirs
  bits_gen = (W64*)malloc((grp->num_gens + T->num_lines[gen_idx] + 1) * sizeof(W64));
  for (a = 0; a < grp->num_gens; a++)
    bits_gen[a] = transport_bits(genus, gen_idx, gen_idx, auts[a]);
  for (idx = 0; idx < T->num_lines[gen_idx]; idx++)
    if (T->cls[gen_idx][idx] == gen_idx) {
      isometry_init_set(gs[num_gs], T->isom[gen_idx][idx]);
      line_of_inverse(lgs[num_gs], gs[num_gs], q);
      bits_gen[grp->num_gens + num_gs] = T->spin[gen_idx][idx] & mask;
      num_gs++;
    }

  if (k == 1)
    seen_set_init(&seen, genus->genus_reps->keys[gen_idx], p, chart);
  else {
    lattice_set_init(&lseen, p, num_nbrs);
    fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
    fmpz_mat_init(s_fmpz, QF_RANK, QF_RANK);
  }
  cap_items = 1024;
  stack = (transport_item*)malloc(cap_items * sizeof(transport_item));

  lc = 0;
  num_roots = 0;
  // the neighbors, as the managers of the orbits (T_p) and of the neighbors (T_{p^2}) give them: the first
  // one of each orbit is tested for isometry, and the others are carried from it
  for (i = 0; i < ((k == 1) ? p : 1); i++) {
    if (k == 1)
      nbr_process_init(nbr_man, genus->genus_reps->keys[gen_idx], p, i, genus->isoms[gen_idx]);
    while ((k == 1) ? !nbr_process_has_ended(nbr_man) : !nbr_data_has_ended(nd)) {
      if (k == 1) {
	vector_set(u, nbr_man->iso_vec);
	is_new = seen_set_insert(&seen, u);
	has_root = is_new && nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man);
      }
      else {
	nbr_data_build_neighbor(nbr_fmpz, s_fmpz, nd);
	for (r = 0; r < QF_RANK; r++)
	  for (c = 0; c < QF_RANK; c++)
	    B[r][c] = (Z64)(((fmpz_get_si(fmpz_mat_entry(s_fmpz, r, c)) % m) + m) % m);
	lattice_key(key, &lseen, (const Z64 (*)[QF_RANK])B);
	is_new = lattice_set_insert(&lseen, key);
	has_root = is_new;
	if (has_root) {
	  square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
	  isometry_init_set_fmpz_mat(s_nbr, s_fmpz, p);
	}
      }
      if (has_root) {
	// a new orbit: its first neighbor is tested for isometry
	int j = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time,
					  &num_isom);
	num_roots++;
	if (k == 1)
	  vector_set(stack[0].v, u);
	else
	  memcpy(stack[0].B, B, sizeof(B));
	stack[0].cls = j;
	isometry_init(stack[0].c);
	isometry_mul(stack[0].c, s_nbr, hash_isom);
	stack[0].bits = transport_bits(genus, gen_idx, j, stack[0].c);
	isometry_clear(hash_isom);
	num_items = 1;
	while (num_items > 0) {
	  item = stack[--num_items];
	  // the entry of this neighbor
	  if (check && (item.bits != transport_bits(genus, gen_idx, item.cls, item.c)))
	    fprintf(stderr, "transport: carried spinor bits differ, class %d -> %d\n", gen_idx, item.cls);
	  hist[(((W64)item.cls) << (genus->spinor->num_primes)) | item.bits]++;
	  lc++;
	  // its images under the generators
	  for (a = 0; a < grp->num_gens + num_gs; a++) {
	    const isometry_struct* g = (a < grp->num_gens) ? auts[a] : gs[a - grp->num_gens];
	    // a neighbor reached when the stack is full is left for later: it is reached again from another
	    // one, or begins a new orbit
	    if (num_items == TRANSPORT_MAX_STACK)
	      continue;
	    if (k == 1) {
	      mat_vec_mod_p(w, g->s, item.v, p);
	      if (!seen_set_insert(&seen, w))
		continue;
	      vector_set(child.v, w);
	    }
	    else {
	      mat_mat_mod(child.B, g->s, (const Z64 (*)[QF_RANK])item.B, m);
	      lattice_key(key, &lseen, (const Z64 (*)[QF_RANK])child.B);
	      if (!lattice_set_insert(&lseen, key))
		continue;
	    }
	    transport_child(&child, &item, a, grp->num_gens, g, lgs, bits_gen, mask, T);
	    if (check)
	      transport_check_child(&child, a, grp->num_gens, gen_idx, p, q, genus);
	    if (num_items == cap_items) {
	      cap_items *= 2;
	      stack = (transport_item*)realloc(stack, cap_items * sizeof(transport_item));
	      if (stack == NULL) {
		fprintf(stderr, "Error: out of memory in hecke_col_all_conductors_transport\n");
		exit(EXIT_FAILURE);
	      }
	    }
	    stack[num_items++] = child;
	  }
	  isometry_clear(item.c);
	}
	square_matrix_clear(nbr);
	isometry_clear(s_nbr);
      }
      if (k == 1)
	nbr_process_advance(nbr_man);
      else
	nbr_data_get_next_neighbor(nd);
    }
    if (k == 1)
      nbr_process_clear(nbr_man);
  }

  if (check)
    fprintf(stderr, "transport stats: class %d, p = %d, k = %d, q = %ld, %ld generators of Aut, %ld g's, %ld roots, %ld neighbors\n",
	    gen_idx, p, k, (long)q, (long)grp->num_gens, (long)num_gs, (long)num_roots, (long)lc);
  free(stack);
  if (k == 1)
    seen_set_clear(&seen);
  else {
    lattice_set_clear(&lseen);
    fmpz_mat_clear(nbr_fmpz);
    fmpz_mat_clear(s_fmpz);
    nbr_data_clear(nd);
  }
  for (a = 0; a < grp->num_gens; a++)
    isometry_clear(auts[a]);
  free(auts);
  for (idx = 0; idx < num_gs; idx++)
    isometry_clear(gs[idx]);
  free(gs);
  free(lgs);
  free(bits_gen);
  aut_grp_clear(grp);

  return lc;
}
