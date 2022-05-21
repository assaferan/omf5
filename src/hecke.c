#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int process_isotropic_vector(matrix_TYP* v, matrix_TYP* w, matrix_TYP* Q,
			     int p, matrix_TYP* b, int* T, int* th61)
{
  matrix_TYP* iso_vec;
  int iso_j;

  iso_j = 0;
  do {
    iso_vec = get_next_isotropic_vector(Q, p, v, w, b, &iso_j);
    if (iso_vec != NULL)
      T[q61_id(q61_nb(Q, p, iso_vec), th61)]++;
  } while ((iso_j != 0) && (iso_j != p));

  return 0;
}

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig)
{
  matrix_TYP *Q, *v, *b, *w_mat;
  int *th61;
  int *w;

  nbrs_data* init;

  /* printf("initializing q61 data\n"); */
  if (init_orig == NULL)
    init = q61_init(p, 2);
  else
    init = init_orig;

  Q = init->Q;
  v = init->v;
  th61 = init->th61;

  /* printf("initialized Q: \n"); */
  /* print_mat(Q); */
  /* printf("isotropic vector: "); */
  /* print_mat(v); */
  
  b = mat_mul(v, Q);
  modp_mat(b,p);

  /* printf("b = v*Q mod p: \n"); */
  /* print_mat(b); */

  w_mat = init_mat(1,5, "");
  w = w_mat->array.SZ[0];
  
  if (i % p == 0) {
    /* printf("Doing [v]..\n"); */
    /* do [v] */
    T[q61_id(q61_nb(Q, p, v), th61)]++;

    w[0] = w[1] = w[2] = w[3] = 0;
    w[4] = 1;

    process_isotropic_vector(v, w_mat, Q, p, b, T, th61);
    
  }

  w[0] = w[1] = w[2] = 0;
  w[3] = 1;
  w[4] = i;

  while (w[0] == 0) {
    process_isotropic_vector(v, w_mat, Q, p, b, T, th61);
    update_pivot(w, p, i);
  }

  return 0;
}
