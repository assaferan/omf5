#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int process_isotropic_vector(neighbor_manager* nbr_man, int* T,
			     int* th61, hash_table* genus)
{

  T[indexof(genus, q61_nb(nbr_man), 0)]++;
  
  return 0;
}

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig, hash_table* genus)
{
  matrix_TYP *Q, *v;
  int *th61;

  nbrs_data* init;
  neighbor_manager nbr_man;

#ifdef DEBUG_LEVEL_FULL
  printf("initializing q61 data\n");
#endif // DEBUG_LEVEL_FULL
  if (init_orig == NULL) {
    init = (nbrs_data*) malloc(sizeof(nbrs_data));
    init_nbrs_data(init, p, 2);
  }
  else
    init = init_orig;

  Q = init->Q;
  v = init->v;
  th61 = init->th61;

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
  printf("isotropic vector: ");
  print_mat(v);
#endif // DEBUG_LEVEL_FULL

  init_nbr_process(&nbr_man, Q, p, i);

  while (!(has_ended(&nbr_man))) {
    process_isotropic_vector(&nbr_man, T, th61, genus);
    advance_nbr_process(&nbr_man);
  }

  free_nbr_process(&nbr_man);
  
  if (init_orig == NULL) {
    free_nbrs_data(init);
    free(init);
  }
 
  return 0;
}
