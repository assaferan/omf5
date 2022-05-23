#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

/*
int process_isotropic_vector(matrix_TYP* v, matrix_TYP* w, matrix_TYP* Q,
			     int p, matrix_TYP* b, int* T, int* th61)
*/
int process_isotropic_vector(neighbor_manager* nbr_man, int* T, int* th61)
{

  T[q61_id(q61_nb(nbr_man->Q, nbr_man->p, nbr_man->iso_vec), th61)]++;
  
  return 0;
}

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig)
{
  matrix_TYP *Q, *v; // , *b, *w_mat;
  int *th61;
  // int *w;

  nbrs_data* init;
  neighbor_manager nbr_man;
  
  /* printf("initializing q61 data\n"); */
  if (init_orig == NULL) {
    init = (nbrs_data*) malloc(sizeof(nbrs_data));
    init_nbrs_data(init, p, 2);
  }
  else
    init = init_orig;

  Q = init->Q;
  v = init->v;
  th61 = init->th61;

  /* printf("initialized Q: \n"); */
  /* print_mat(Q); */
  /* printf("isotropic vector: "); */
  /* print_mat(v); */

  init_nbr_process(&nbr_man, Q, p, i);

  // while (nbr_man.w->array.SZ[0][0] == 0) {
  while (!(has_ended(&nbr_man))) {
     process_isotropic_vector(&nbr_man, T, th61);
     advance_nbr_process(&nbr_man);
  }

  if (init_orig == NULL)
    free_nbrs_data(init);
 
  return 0;
}
