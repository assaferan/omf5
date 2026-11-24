/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : tq_table.c
 *
 * Description: The table of T_q (see tq_table.h).
 *
 **********************************************************
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>

#include "hash.h"
#include "nbr_data.h"
#include "spinor.h"
#include "tq_table.h"

slong tq_table_line_index(const vector_t w, int q)
{
  vector_t v;
  slong idx;
  int i;

  vector_set(v, w);
  vector_mod_p(v, q);
  for (i = 0; (i < QF_RANK) && (v[i] == 0); i++);
  if (i == QF_RANK)
    return -1;
  normalize_mod_p(v, q);
  idx = 0;
  for (i = QF_RANK - 1; i >= 0; i--)
    idx = idx * q + v[i];
  return idx;
}

slong tq_table_pos(const tq_table_t T, slong j, const vector_t w)
{
  slong idx = tq_table_line_index(w, T->q);
  return (idx < 0) ? -1 : T->pos[j][idx];
}

void tq_table_init(tq_table_t T, const genus_t genus, int q)
{
  slong j, idx, h, n, max_lines;
  int i, k;
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_fmpz, nbr_isom;
  square_matrix_t nbr;
  isometry_t s_nbr, hash_isom, s_inv, s_norm;
  vector_t w;
  int num_isom;
  double theta_time, isom_time;

  assert(fmpz_fdiv_ui(genus->disc, q) != 0);

  h = genus->genus_reps->num_stored;
  T->q = q;
  T->num_classes = h;
  T->size = 1;
  for (k = 0; k < QF_RANK; k++)
    T->size *= q;
  T->pos = (int**)malloc(h * sizeof(int*));
  T->cls = (int**)malloc(h * sizeof(int*));
  T->isom = (isometry_t**)malloc(h * sizeof(isometry_t*));
  T->spin = (W64**)malloc(h * sizeof(W64*));
  T->line = (vector_t**)malloc(h * sizeof(vector_t*));
  T->num_lines = (slong*)malloc(h * sizeof(slong));
  // the isotropic lines of a form mod q which is nondegenerate or has a radical of dimension 1 (q does
  // not divide the discriminant, so the radical is 0, or a line at q = 2): at most q^3 + q^2 + q + 1
  max_lines = (slong)q * q * q + (slong)q * q + q + 1;

  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  num_isom = 0;
  theta_time = isom_time = 0;

  for (j = 0; j < h; j++) {
    T->pos[j] = (int*)malloc(T->size * sizeof(int));
    T->cls[j] = (int*)malloc(max_lines * sizeof(int));
    T->isom[j] = (isometry_t*)malloc(max_lines * sizeof(isometry_t));
    T->spin[j] = (W64*)malloc(max_lines * sizeof(W64));
    T->line[j] = (vector_t*)malloc(max_lines * sizeof(vector_t));
    for (idx = 0; idx < T->size; idx++)
      T->pos[j][idx] = -1;
    T->num_lines[j] = 0;

    nbr_data_init(nbr_man, genus->genus_reps->keys[j], q, 1);
    while (!nbr_data_has_ended(nbr_man)) {
      // the isotropic line, in the coordinates of keys[j] (X_skew is its lift mod q^2)
      for (k = 0; k < QF_RANK; k++)
	w[k] = (Z64)nmod_mat_entry(nbr_man->X_skew, 0, k);
      idx = tq_table_line_index(w, q);
      if ((idx < 0) || (T->pos[j][idx] != -1) || (T->num_lines[j] == max_lines)) {
	fprintf(stderr, "Error: the table of T_%d has a line which is 0 mod q, or twice, or too many\n", q);
	exit(EXIT_FAILURE);
      }
      n = T->num_lines[j];
      T->pos[j][idx] = n;
      vector_set(T->line[j][n], w);
      vector_mod_p(T->line[j][n], q);
      normalize_mod_p(T->line[j][n], q);

      // as in process_isotropic_vector_nbr_data_all_conductors
      nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
      square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
      isometry_init_set_fmpz_mat(s_nbr, nbr_isom, q);
      i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
      assert((0 <= i) && (i < h));

      // isom = s_nbr * hash_isom : keys[j] -> keys[i]
      isometry_init(T->isom[j][n]);
      isometry_mul(T->isom[j][n], s_nbr, hash_isom);
      T->cls[j][n] = i;

      // the spinor bits, as for spin_vals: of isoms[j] * isom * isoms[i]^(-1), an isometry of keys[0]
      isometry_init(s_norm);
      isometry_mul(s_norm, genus->isoms[j], T->isom[j][n]);
      isometry_inv(s_inv, genus->isoms[i]);
      isometry_muleq_right(s_norm, s_inv);
      T->spin[j][n] = (((W64)i) << (genus->spinor->num_primes)) | spinor_norm_isom(genus->spinor, s_norm);
      T->num_lines[j]++;

      isometry_clear(s_norm);
      isometry_clear(s_inv);
      isometry_clear(s_nbr);
      isometry_clear(hash_isom);
      square_matrix_clear(nbr);
      nbr_data_get_next_neighbor(nbr_man);
    }
    nbr_data_clear(nbr_man);
  }

  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);
  return;
}

void tq_table_clear(tq_table_t T)
{
  slong j, n;

  for (j = 0; j < T->num_classes; j++) {
    for (n = 0; n < T->num_lines[j]; n++)
      isometry_clear(T->isom[j][n]);
    free(T->pos[j]);
    free(T->cls[j]);
    free(T->isom[j]);
    free(T->spin[j]);
    free(T->line[j]);
  }
  free(T->pos);
  free(T->cls);
  free(T->isom);
  free(T->spin);
  free(T->line);
  free(T->num_lines);
  return;
}
