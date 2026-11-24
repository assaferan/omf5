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

#include <flint/ulong_extras.h>

#include "arith.h"
#include "aut_grp.h"
#include "hash.h"
#include "isometry.h"
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

static bool line_set_contains(const line_set_struct* S, W64 key)
{
  return S->keys[line_set_slot(S, key)] == key;
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

static bool seen_set_contains(const seen_set_struct* S, const vector_t v)
{
  W64 bit;
  vector_t w;

  if (S->k < 0) {
    vector_set(w, v);
    return line_set_contains(&(S->hash), line_key(w, S->p));
  }
  bit = seen_set_bit(S, v);
  if (S->check) {
    vector_set(w, v);
    if (line_set_contains(&(S->shadow), line_key(w, S->p)) != (bool)((S->bits[bit >> 6] >> (bit & 63)) & 1))
      fprintf(stderr, "transport: the bitmap of the lines differs from the hash set, p = %ld\n", (long)S->p);
  }
  return (S->bits[bit >> 6] >> (bit & 63)) & 1;
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

typedef struct {
  vector_t v;
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

slong hecke_col_all_conductors_transport(slong* hist, int p, int gen_idx, const genus_t genus,
					 const tq_table_t T)
{
  neighbor_manager_t nbr_man;
  aut_grp_t grp;
  square_matrix_t nbr;
  isometry_t *auts, *gs, s_nbr, hash_isom, tmp;
  W64 *bits_gen, mask;
  bool check;
  vector_t *lgs, u, w;
  slong num_gs, idx, lc, num_items, cap_items, num_roots;
  int a, i, j2, num_isom;
  double theta_time, isom_time;
  Z64 q = T->q;
  int chart;
  seen_set_struct seen;
  transport_item* stack;
  transport_item item, child;

  assert(q != p);
  assert(gen_idx < T->num_classes);
  // the lines are kept in a bitmap when the Gram matrix of the class has a diagonal entry prime to p, and
  // otherwise in a hash set, which takes more memory: then, if it does not fit, the column is left to the
  // orbits (-1)
  chart = seen_set_chart(genus->genus_reps->keys[gen_idx], p);
  if ((chart < 0) && (transport_bits_needed(p, false) > (((W64)transport_memory_mb) << 23)))
    return -1;
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

  seen_set_init(&seen, genus->genus_reps->keys[gen_idx], p, chart);
  cap_items = 1024;
  stack = (transport_item*)malloc(cap_items * sizeof(transport_item));

  lc = 0;
  num_roots = 0;
  for (i = 0; i < p; i++) {
    nbr_process_init(nbr_man, genus->genus_reps->keys[gen_idx], p, i, genus->isoms[gen_idx]);
    while (!nbr_process_has_ended(nbr_man)) {
      vector_set(u, nbr_man->iso_vec);
      if (seen_set_insert(&seen, u)) {
	// a new orbit: its first neighbor is tested for isometry
	num_roots++;
	if (nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man)) {
	  int j = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time,
					    &num_isom);
	  num_items = 0;
	  vector_set(stack[0].v, u);
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
	      mat_vec_mod_p(w, g->s, item.v, p);
	      // a line reached when the stack is full is left for later: it is reached again from another
	      // line, or begins a new orbit
	      if (seen_set_contains(&seen, w) || (num_items == TRANSPORT_MAX_STACK))
		continue;
	      seen_set_insert(&seen, w);
	      vector_set(child.v, w);
	      isometry_init(child.c);
	      if (a < grp->num_gens) {
		// c_{av} = a c_v
		isometry_mul(child.c, g, item.c);
		child.cls = item.cls;
		child.bits = item.bits ^ bits_gen[a];
	      }
	      else {
		// c_{gv} = g c_v t, t from the table at (cls, line of c_v^(-1) l_g)
		vector_t lw;
		line_mod_q(lw, item.c->s_inv, item.c->inv_denom, lgs[a - grp->num_gens], q);
		idx = tq_table_pos(T, item.cls, lw);
		if (idx < 0) {
		  fprintf(stderr, "Error: a carried line is not in the table of T_%ld\n", (long)q);
		  exit(EXIT_FAILURE);
		}
		j2 = T->cls[item.cls][idx];
		isometry_init(tmp);
		isometry_mul(tmp, g, item.c);
		isometry_mul(child.c, tmp, T->isom[item.cls][idx]);
		isometry_clear(tmp);
		child.cls = j2;
		child.bits = item.bits ^ bits_gen[a] ^ (T->spin[item.cls][idx] & mask);
	      }
	      isometry_reduce(child.c);
	      if (check) {
		square_matrix_t nb;
		isometry_t sn;
		vector_t uu;
		bool okk;
		// the neighbor at w, built directly, must be the image of keys[cls] under c
		okk = isometry_is_isom(child.c, genus->genus_reps->keys[gen_idx],
				       genus->genus_reps->keys[child.cls]);
		vector_set(uu, w);
		(void)uu; (void)nb; (void)sn;
		if (!okk)
		  fprintf(stderr, "transport: not an isometry, generator %d (%s), class %d -> %d\n", a,
			  (a < grp->num_gens) ? "aut" : "g", item.cls, child.cls);
		{
		  // integrality away from p and q: the denominators of c and c^(-1) are products of p and q
		  Z64 d1 = child.c->denom, d2 = child.c->inv_denom;
		  while (d1 % p == 0) d1 /= p;
		  while (d1 % q == 0) d1 /= q;
		  while (d2 % p == 0) d2 /= p;
		  while (d2 % q == 0) d2 /= q;
		  if ((d1 != 1 && d1 != -1) || (d2 != 1 && d2 != -1))
		    fprintf(stderr, "transport: denominators %ld, %ld, generator %d (%s), class %d -> %d\n",
			    (long)child.c->denom, (long)child.c->inv_denom, a,
			    (a < grp->num_gens) ? "aut" : "g", item.cls, child.cls);
		}
	      }
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
	}
	square_matrix_clear(nbr);
	isometry_clear(s_nbr);
      }
      nbr_process_advance(nbr_man);
    }
    nbr_process_clear(nbr_man);
  }

  if (getenv("OMF5_TRANSPORT_CHECK") != NULL)
    fprintf(stderr, "transport stats: class %d, p = %d, q = %ld, %ld generators of Aut, %ld g's, %ld roots, %ld lines\n",
	    gen_idx, p, (long)q, (long)grp->num_gens, (long)num_gs, (long)num_roots, (long)lc);
  free(stack);
  seen_set_clear(&seen);
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
