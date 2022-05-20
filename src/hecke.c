#include "arith.h"
#include "hecke.h"
#include "neighbor.h"

int process_isotropic_vector(matrix_TYP* v, matrix_TYP* w_mat, matrix_TYP* Q,
			     int p, matrix_TYP* b, int* T, int* th61)
{
  matrix_TYP *t_mat, *n_mat, *tmp;
  int j, t, n, n_inv, dummy;
  
  n_mat = mat_mul(mat_mul(w_mat, Q), tr_pose(w_mat));
  n = n_mat->array.SZ[0][0]/2 % p;
  t_mat = mat_mul(w_mat, tr_pose(b));
  t = t_mat->array.SZ[0][0] % p;

  /* printf("v = "); */
  /* print_mat(v); */
  
  /* printf("w = "); */
  /* print_mat(w_mat); */

  /* printf("n = %d, t = %d\n", n, t); */
    
  if (n) {
    if (t) {
      /* do [n*v-t*w] */
      gcdext(n, p, &n_inv, &dummy);
      tmp = imat_add(v, w_mat, 1, -t*n_inv);
      /* printf("tmp = "); */
      /* print_mat(tmp); */
      /* This differs from the gp script */
      /* modp_mat(tmp, p); */
      for (j = 0; j < 5; j++)
	tmp->array.SZ[0][j] %= p;
      /* printf("tmp mod p = "); */
      /* print_mat(tmp); */
      T[q61_id(q61_nb(Q, p, tmp), th61)]++;
    }
  }
  else {
    /* do [w] */
    T[q61_id(q61_nb(Q, p, w_mat), th61)]++;
    if (!t) {
      for (j = 1; j < p; j++)
	/* do [v+j*w] */
	T[q61_id(q61_nb(Q, p, imat_add(v,w_mat,1,j)), th61)]++;
    }
  }

  return 0;
}

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig)
{
  matrix_TYP *Q, *v, *b, *w_mat;
  int *th61;
  int *w, j, k;

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

  process_isotropic_vector(v, w_mat, Q, p, b, T, th61);

  for (j = 0; j < p; j++) {
    w[0] = w[1] = 0;
    w[2] = 1;
    w[3] = i;
    w[4] = j;

    process_isotropic_vector(v, w_mat, Q, p, b, T, th61);

    for (k = 0; k < p; k++) {
      w[0] = 0;
      w[1] = 1;
      w[2] = i;
      w[3] = j;
      w[4] = k;

      process_isotropic_vector(v, w_mat, Q, p, b, T, th61);
      
    }
  }

  return 0;
}
