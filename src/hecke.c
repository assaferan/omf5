#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int process_isotropic_vector(neighbor_manager* nbr_man, int* T, hash_table* genus)
{

  T[indexof(genus, q61_nb(nbr_man), 0)]++;
  
  return 0;
}

int process_neighbour_chunk(int* T, int p, int i, int gen_idx, hash_table* genus)
{
  matrix_TYP *Q, *v;
 
  neighbor_manager nbr_man;

  Q = genus->keys[gen_idx];
  v = get_isotropic_vector(Q, p);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
  printf("isotropic vector: ");
  print_mat(v);
#endif // DEBUG_LEVEL_FULL

  init_nbr_process(&nbr_man, Q, p, i);

  while (!(has_ended(&nbr_man))) {
    process_isotropic_vector(&nbr_man, T, genus);
    advance_nbr_process(&nbr_man);
  }

  free_nbr_process(&nbr_man);
 
  return 0;
}

// assumes T is initialized to zeros
void hecke_col(int* T, int p, int gen_idx, hash_table* genus)
{
  int num;
  
  for (num = 0; num < p; num++) {
    process_neighbour_chunk(T, p, num, gen_idx, genus);
  }

  return;
}

matrix_TYP* hecke_matrix(hash_table* genus, int p)
{
  matrix_TYP* hecke;
  int gen_idx;

  hecke = init_mat(genus->num_stored, genus->num_stored, "");

  for (gen_idx = 0; gen_idx < genus->num_stored; gen_idx++)
    hecke_col(hecke->array.SZ[gen_idx], p, gen_idx, genus);
  
  return hecke;
}
