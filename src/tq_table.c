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
  normalize_mod_p(v, q);
  idx = 0;
  for (i = QF_RANK - 1; i >= 0; i--)
    idx = idx * q + v[i];
  return idx;
}

void tq_table_init(tq_table_t T, const genus_t genus, int q)
{
  slong j, idx, h;
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
  T->cls = (int**)malloc(h * sizeof(int*));
  T->isom = (isometry_t**)malloc(h * sizeof(isometry_t*));
  T->spin = (W64**)malloc(h * sizeof(W64*));
  T->num_lines = (slong*)malloc(h * sizeof(slong));

  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  num_isom = 0;
  theta_time = isom_time = 0;

  for (j = 0; j < h; j++) {
    T->cls[j] = (int*)malloc(T->size * sizeof(int));
    T->isom[j] = (isometry_t*)malloc(T->size * sizeof(isometry_t));
    T->spin[j] = (W64*)malloc(T->size * sizeof(W64));
    for (idx = 0; idx < T->size; idx++)
      T->cls[j][idx] = -1;
    T->num_lines[j] = 0;

    nbr_data_init(nbr_man, genus->genus_reps->keys[j], q, 1);
    while (!nbr_data_has_ended(nbr_man)) {
      // the isotropic line, in the coordinates of keys[j] (X_skew is its lift mod q^2)
      for (k = 0; k < QF_RANK; k++)
	w[k] = (Z64)nmod_mat_entry(nbr_man->X_skew, 0, k);
      idx = tq_table_line_index(w, q);
      assert(T->cls[j][idx] == -1);

      // as in process_isotropic_vector_nbr_data_all_conductors
      nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
      square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
      isometry_init_set_fmpz_mat(s_nbr, nbr_isom, q);
      i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
      assert((0 <= i) && (i < h));

      // isom = s_nbr * hash_isom : keys[j] -> keys[i]
      isometry_init(T->isom[j][idx]);
      isometry_mul(T->isom[j][idx], s_nbr, hash_isom);
      T->cls[j][idx] = i;

      // the spinor bits, as for spin_vals: of isoms[j] * isom * isoms[i]^(-1), an isometry of keys[0]
      isometry_init(s_norm);
      isometry_mul(s_norm, genus->isoms[j], T->isom[j][idx]);
      isometry_inv(s_inv, genus->isoms[i]);
      isometry_muleq_right(s_norm, s_inv);
      T->spin[j][idx] = (((W64)i) << (genus->spinor->num_primes)) | spinor_norm_isom(genus->spinor, s_norm);
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
  slong j, idx;

  for (j = 0; j < T->num_classes; j++) {
    for (idx = 0; idx < T->size; idx++)
      if (T->cls[j][idx] >= 0)
	isometry_clear(T->isom[j][idx]);
    free(T->cls[j]);
    free(T->isom[j]);
    free(T->spin[j]);
  }
  free(T->cls);
  free(T->isom);
  free(T->spin);
  free(T->num_lines);
  return;
}
