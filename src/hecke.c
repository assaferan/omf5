#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int process_isotropic_vector(neighbor_manager* nbr_man, int* T, hash_table* genus)
{

  int i;
  
  i = indexof(genus, build_nb(nbr_man), 0);
  
  T[i]++;
  
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

void get_hecke_ev(nf_elem_t e, hash_table* genus, eigenvalues* evs, int p, int ev_idx)
{
  int* a;
  int num, k, pivot;
  nf_elem_t prod;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(genus->num_stored * sizeof(int));
  for (num = 0; num < genus->num_stored; num++)
    a[num] = 0;

  pivot = 0;
  for (k = 0; k < evs->dim; k++)
    if (nf_elem_is_zero(evs->eigenvecs[ev_idx][k], evs->nfs[ev_idx])) {
      pivot++;
    }
  
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
  
  hecke_col(a, p, pivot, genus);

  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;

  nf_elem_zero(e, evs->nfs[ev_idx]);
  for (k = 0; k < evs->dim; k++) {
    nf_elem_set_si(prod, a[k], evs->nfs[ev_idx]);
    nf_elem_mul(prod, prod, evs->eigenvecs[ev_idx][k], evs->nfs[ev_idx]);
    nf_elem_add(e, e, prod, evs->nfs[ev_idx]);
  }
  nf_elem_div(e, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);
  free(a);
  
  return;
}
