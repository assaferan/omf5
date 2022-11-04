#include <assert.h>
#include <inttypes.h>

#ifdef DEBUG
#include <carat/symm.h>
#endif // DEBUG

#include "arith.h"
#include "matrix_tools.h"
#include "neighbor.h"
#include "typedefs.h"

/* Compute one p-neighbour for Q_orig corresponding to vector x 
 * Returns false if the vector is in the radical
*/

bool is_isotropic_vector(const vector_t x, const square_matrix_t Q, Z64 p)
{
  vector_t Qx;
  Z64 xQx;

  square_matrix_mul_vec_left(Qx, x, Q);
  xQx = scalar_product(x, Qx);

  return (xQx % p == 0);
}

bool nbr_process_build_nb(square_matrix_t Q, const neighbor_manager_t nbr_man)
{
  vector_t x, Qx;
  Z64 xQx;
  Z64 q, y;
  int row, col;
  Z64 a1, a2, a3, a4;
  
#ifdef DEBUG_LEVEL_FULL
  printf("Computing p-neighbor for Q: \n");
  square_matrix_print(nbr_man->Q);

  printf("Corresponding to isotropic vector ");
  vector_print(nbr_man->iso_vec);
#endif // DEBUG_LEVEL_FULL

  square_matrix_set(Q, nbr_man->Q);
  square_matrix_mul_vec_left(Qx, nbr_man->iso_vec, Q);

  /* From now until the end, Q is computed only upper triangular */

#ifdef DEBUG_LEVEL_FULL
  printf("Qx = \n");
  vector_print(Qx);
#endif // DEBUG_LEVEL_FULL
  
  vector_set(x, nbr_man->iso_vec);
  xQx = scalar_product(x, Qx);
  
  if (x[0] == 1) {
      /* M[,1] = x~ */
      
      Q[0][0] = xQx;
      for (col = 1; col < QF_RANK; col++)
	Q[0][col] = Qx[col];
#ifdef DEBUG_LEVEL_FULL
      printf("put x as first vector, now Q = \n");
      square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  }
  else {
    if (x[1] == 1) {
      /* M[,2] = M[,1] ; M[,1] = x~ */
      for (col = 2; col < QF_RANK; col++)
	Q[1][col] = Q[0][col];
      Q[1][1] = Q[0][0];
      for (col = 2; col < QF_RANK; col++)
	Q[0][col] = Qx[col];
      Q[0][1] = Qx[0];
      Q[0][0] = xQx;
    }
    else {
      if (x[2] == 1) {
	/* M[,3] = M[,1] ; M[,1] = x~ */
	for (col = 3; col < QF_RANK; col++)
	  Q[2][col] = Q[0][col];
	Q[2][2] = Q[0][0];
	Q[1][2] = Q[0][1];
	for (col = 3; col < QF_RANK; col++)
	  Q[0][col] = Qx[col];
	Q[0][2] = Qx[0];
	Q[0][1] = Qx[1];
	Q[0][0] = xQx;
      }
      else {
	if (x[3] == 1) {
	  /* M[,4] = M[,1] ; M[,1] = x~ */
	  Q[3][4] = Q[0][4];
	  Q[3][3] = Q[0][0];
	  Q[2][3] = Q[0][2];
	  Q[1][3] = Q[0][1];
	  Q[0][4] = Qx[4];
	  Q[0][3] = Qx[0];
	  Q[0][2] = Qx[2];
	  Q[0][1] = Qx[1];
	  Q[0][0] = xQx;
	}
	else {
	  assert(x[4]==1);
	  /* M[,N] = M[,1] ; M[,1] = x~ */
	  Q[4][4] = Q[0][0];
	  Q[3][4] = Q[0][3];
	  Q[2][4] = Q[0][2];
	  Q[1][4] = Q[0][1];
	  Q[0][4] = Qx[0];
	  Q[0][3] = Qx[3];
	  Q[0][2] = Qx[2];
	  Q[0][1] = Qx[1];
	  Q[0][0] = xQx;
	}
      }
    }
  }

  assert(!(Q[0][0] % nbr_man->p));
  
  if (!(Q[0][4] % nbr_man->p)) {
    if (Q[0][3] % nbr_man->p) {
      /* M[,4] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
      printf("swapping rows 4 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
      for (row = 0; row < 3; row++) {
	square_matrix_swap_elts(Q, row, 3, row, 4);
      }
      square_matrix_swap_elts(Q, 3, 3, 4, 4);
#ifdef DEBUG_LEVEL_FULL
      square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
    }
    else {
      if (Q[0][2] % nbr_man->p) {
	/* M[,3] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
	printf("swapping rows 3 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
	square_matrix_swap_elts(Q, 0, 2, 0, 4);
	square_matrix_swap_elts(Q, 1, 2, 1, 4);
	square_matrix_swap_elts(Q, 2, 2, 4, 4);
	square_matrix_swap_elts(Q, 2, 3, 3, 4);
#ifdef DEBUG_LEVEL_FULL
	square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
      }
      else {
	if (Q[0][1] % nbr_man->p) {
	  /* M[,2] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
	  printf("swapping rows 2 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
	  square_matrix_swap_elts(Q, 0, 1, 0, 4);
	  square_matrix_swap_elts(Q, 1, 1, 4, 4);
	  square_matrix_swap_elts(Q, 1, 2, 2, 4);
	  square_matrix_swap_elts(Q, 1, 3, 3, 4);
#ifdef DEBUG_LEVEL_FULL
	  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
	}
	else {
	  /* a singular zero, return Q as-is*/
	  /* Re-symmetrize */
	  square_matrix_resymmetrize(Q);
	  return false;
	}
      }
    }
  }
  
  /* M[,1] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[1] <-- R[1] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] /= nbr_man->p;
  gcdext(Q[0][4], nbr_man->p, &q, &y);
  a1 = -Q[0][0]*q/2 % nbr_man->p;
  if (a1 < 0)
    a1 += nbr_man->p;
#ifdef DEBUG_LEVEL_FULL
  printf("q = %" PRId64 ", a1 = %" PRId64 " (y = %"  PRId64 ")\n", q, a1, y);
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] += a1*(2*Q[0][4]+nbr_man->p*a1*Q[4][4]);
  Q[0][1] += nbr_man->p*a1*Q[1][4];
  Q[0][2] += nbr_man->p*a1*Q[2][4];
  Q[0][3] += nbr_man->p*a1*Q[3][4];
  Q[0][4] += nbr_man->p*a1*Q[4][4];
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,2] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[2] <-- R[2] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  a2 = -Q[0][1]*q % nbr_man->p;
  if (a2 < 0)
    a2 += nbr_man->p;
  
  Q[0][1] += a2*Q[0][4];
  Q[1][1] += a2*(2*Q[1][4]+a2*Q[4][4]);
  Q[1][2] += a2*Q[2][4];
  Q[1][3] += a2*Q[3][4];
  Q[1][4] += a2*Q[4][4];
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,3] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[3] <-- R[3] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
   
  a3 = -Q[0][2]*q % nbr_man->p;
  if (a3 < 0)
    a3 += nbr_man->p;
  
  Q[0][2] += a3*Q[0][4];
  Q[1][2] += a3*Q[1][4];
  Q[2][2] += a3*(2*Q[2][4]+a3*Q[4][4]);
  Q[2][3] += a3*Q[3][4];
  Q[2][4] += a3*Q[4][4];

#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,4] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[4] <-- R[4] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  a4 = -Q[0][3]*q % nbr_man->p;
  if (a4 < 0)
    a4 += nbr_man->p;
  
  Q[0][3] += a4*Q[0][4];
  Q[1][3] += a4*Q[1][4];
  Q[2][3] += a4*Q[2][4];
  Q[3][3] += a4*(2*Q[3][4]+a4*Q[4][4]);
  Q[3][4] += a4*Q[4][4];

#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,1] /= p */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[1] <-- R[1]/p, now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] /= nbr_man->p;
  Q[0][1] /= nbr_man->p;
  Q[0][2] /= nbr_man->p;
  Q[0][3] /= nbr_man->p;

#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,5] *= p */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[5] <-- R[5]/p, now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[1][4] *= nbr_man->p;
  Q[2][4] *= nbr_man->p;
  Q[3][4] *= nbr_man->p;
  Q[4][4] *= (nbr_man->p)*(nbr_man->p);

#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* Re-symmetrize */
  square_matrix_resymmetrize(Q);

#ifdef DEBUG_LEVEL_FULL
  printf("Resulting neighbor is: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  assert(square_matrix_is_positive_definite(Q));
#endif // DEBUG
  
  return true;
};

/* find an isotropic vector for Q mod p */
/* return a row vector 1xn */
bool get_isotropic_vector(vector_t x, const square_matrix_t Q, Z64 p)
{
  vector_t Qx;
  Z64 v1, v2, v3, v4, n;

  /* find one zero "v0" */
  /* TODO : We have a much better way to do this.... */
  x[0] = 1;
  for (v1 = 0; v1 < p; v1++) {
    x[1] = v1;
    for (v2 = 0; v2 < p; v2++) {
      x[2] = v2;
      for (v3 = 0; v3 < p; v3++) {
	x[3] = v3;
	for (v4 = 0; v4 < p; v4++) {
	  x[4] = v4;
	  square_matrix_mul_vec_left(Qx, x, Q);
	  n = (scalar_product(x, Qx) / 2) % p;
	  if (!n)
	    return true;
	}
      }
    }
  }
  
  //  printf("Didn't find an isotropic vector!?\n");
  return false;
}

/* update the pivot vector v */
void update_pivot(vector_t v, Z64 p, int i)
{
  int pivot, lsb;
  
  pivot = 0;
  while (v[pivot] == 0) {
    pivot++;
  }

  assert(pivot != 0);
  assert(pivot != QF_RANK);
  assert(v[pivot] == 1);

  lsb = QF_RANK-1;
  v[lsb]++;

  /* carry */
  while ((v[lsb] == p) && (pivot+1 < lsb)) {
    v[lsb] = 0;
    lsb--;
    v[lsb]++;
  }
  /* We have completed this round
   */
  if (lsb <= pivot+1) {
    if (pivot < QF_RANK-1)
      v[pivot+1] = 0;
    v[pivot] = i;
    v[pivot-1] = 1;
  }

  return;
}

/* get isotropic vector, correposnding to the pivot vector w. updates nbr_man->iso_vec */

bool get_next_isotropic_vector(neighbor_manager_t nbr_man)
{
  vector_t wQ;
 
  Z64 t, n, n_inv, dummy;

  if (nbr_man->iso_j) {
    vector_lin_comb(nbr_man->iso_vec, nbr_man->v, nbr_man->w, 1, nbr_man->iso_j);
    (nbr_man->iso_j)++;
    return true;
  }

  square_matrix_mul_vec_left(wQ, nbr_man->w, nbr_man->Q);
  n = (scalar_product(wQ, nbr_man->w) / 2) % nbr_man->p;
  t = scalar_product(nbr_man->b, nbr_man->w) % nbr_man->p;

#ifdef DEBUG_LEVEL_FULL
  printf("v = ");
  vector_print(nbr_man->v);
  
  printf("w = ");
  vector_print(nbr_man->w);

  printf("n = %" PRId64 ", t = %" PRId64 "\n", n, t);
#endif // DEBUG_LEVEL_FULL
    
  if (n) {
    if (t) {
      /* do [n*v-t*w] */
      gcdext(n, nbr_man->p, &n_inv, &dummy);
      vector_lin_comb(nbr_man->iso_vec, nbr_man->v, nbr_man->w, 1, -t*n_inv);
#ifdef DEBUG_LEVEL_FULL
      printf("iso_vec = ");
      vector_print(nbr_man->iso_vec);
#endif // DEBUG_LEVEL_FULL
      vector_mod_p(nbr_man->iso_vec, nbr_man->p);
#ifdef DEBUG_LEVEL_FULL
      printf("n*v-t*w = ");
      vector_print(nbr_man->iso_vec);
#endif // DEBUG_LEVEL_FULL
      return true;
    }
  }
  else {
    if (nbr_man->iso_j == 0)
      vector_set(nbr_man->iso_vec, nbr_man->w);
    else
      vector_lin_comb(nbr_man->iso_vec, nbr_man->v, nbr_man->w, 1, nbr_man->iso_j);
#ifdef DEBUG_LEVEL_FULL
    printf("v+%" PRId64 "*w = ", nbr_man->iso_j);
    vector_print(nbr_man->iso_vec);
#endif // DEBUG_LEVEL_FULL
    if (!t)
      (nbr_man->iso_j)++;
    return true;
  }

  return false;
}

void nbr_process_init(neighbor_manager_t nbr_man, const square_matrix_t Q, Z64 p, int i)
{
  aut_grp_t grp;
  bool found;
  
  square_matrix_set(nbr_man->Q, Q);
  nbr_man->p = p;
  get_isotropic_vector(nbr_man->v, Q, p);
  square_matrix_mul_vec_left(nbr_man->b, nbr_man->v, Q);
  vector_mod_p(nbr_man->b,p);

  nbr_man->i = i;
  nbr_man->iso_j = 0;
  vector_zero(nbr_man->w);
  if (i == 0)
    vector_set(nbr_man->iso_vec, nbr_man->v);
  else {
    nbr_man->w[0] = nbr_man->w[1] = nbr_man->w[2] = 0;
    nbr_man->w[3] = 1;
    nbr_man->w[4] = i;
    found = false;
    while ((nbr_man->w[0] == 0) && (!found)) {
      found = get_next_isotropic_vector(nbr_man);
      if (!found)
	update_pivot(nbr_man->w, nbr_man->p, nbr_man->i);
    }
  }

  nbr_man->first_iter = (i == 0);

  // It could be that there are no isotropic vectors in this batch
  //  assert(is_isotropic_vector(nbr_man->iso_vec, nbr_man->Q, nbr_man->p));

  aut_grp_init_square_matrix(grp, nbr_man->Q);

  nbr_man->auts = (square_matrix_t*)malloc((grp->order / 2) * sizeof(square_matrix_t));
  aut_grp_get_special_elements(nbr_man->auts,grp);
  nbr_man->num_auts = grp->order / 2;
  
  aut_grp_clear(grp);
  return;
}

void nbr_process_advance(neighbor_manager_t nbr_man)
{
  bool found;
  
  // This is the first iteration
  if (nbr_man->first_iter) {
    nbr_man->w[0] = nbr_man->w[1] = nbr_man->w[2] = nbr_man->w[3] = 0;
    nbr_man->w[4] = 1;
    nbr_man->first_iter = false;
  }
  else {
    
    // finished looping over v+jw
    if (nbr_man->iso_j == nbr_man->p) {
      nbr_man->iso_j = 0; 
    }
    
    // update to next pivot
    if (nbr_man->iso_j == 0) {
      update_pivot(nbr_man->w, nbr_man->p, nbr_man->i);
    }
  }

  found = false;
  while ((nbr_man->w[0] == 0) && (!found)) {
    found = get_next_isotropic_vector(nbr_man);

    if (!found) {
      // finished looping over v+jw
      if (nbr_man->iso_j == nbr_man->p) {
	nbr_man->iso_j = 0; 
      }
	
      // update to next pivot
      if (nbr_man->iso_j == 0) {
	update_pivot(nbr_man->w, nbr_man->p, nbr_man->i);
      }
    }
  }

  //  assert(is_isotropic_vector(nbr_man->iso_vec, nbr_man->Q, nbr_man->p));

  return;
}

bool nbr_process_has_ended(const neighbor_manager_t nbr_man)
{
  return (nbr_man->w[0] != 0);
}

void nbr_process_clear(neighbor_manager_t nbr_man)
{
  slong i;
  
  for (i = 0; i < nbr_man->num_auts; i++)
    square_matrix_clear(nbr_man->auts[i]);

  free(nbr_man->auts);
  return;
}

// right now we copy and paste and just modify in order to not slow down the fast case
bool nbr_process_build_nb_and_isom(square_matrix_t Q, isometry_t s, const neighbor_manager_t nbr_man)
{
  vector_t x, Qx;
  Z64 xQx;
  Z64 q, y;
  int row, col;
  Z64 a1, a2, a3, a4;
#ifdef DEBUG
  square_matrix_t full_Q;
  isometry_t s_inv;
  
#ifdef DEBUG_LEVEL_FULL
  printf("Computing p-neighbor for Q: \n");
  square_matrix_print(nbr_man->Q);

  printf("Corresponding to isotropic vector ");
  vector_print(nbr_man->iso_vec);
#endif // DEBUG_LEVEL_FULL
#endif // DEBUG
  
  isometry_init(s);
  square_matrix_set(Q, nbr_man->Q);
  square_matrix_mul_vec_left(Qx, nbr_man->iso_vec, Q);

  /* From now until the end, Q is computed only upper triangular */

#ifdef DEBUG_LEVEL_FULL
  printf("Qx = \n");
  vector_print(Qx);
#endif // DEBUG_LEVEL_FULL
  
  vector_set(x, nbr_man->iso_vec);
  xQx = scalar_product(x, Qx);
  
  if (x[0] == 1) {
      /* M[,1] = x~ */
      
      Q[0][0] = xQx;
      for (col = 1; col < QF_RANK; col++)
	Q[0][col] = Qx[col];

      isometry_replace_vec(s, 0, x);
      
#ifdef DEBUG
      square_matrix_set(full_Q, Q);
      square_matrix_resymmetrize(full_Q);
      assert(isometry_is_isom(s, nbr_man->Q, full_Q));
      isometry_inv(s_inv, s);
      assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#ifdef DEBUG_LEVEL_FULL
      printf("put x as first vector, now Q = \n");
      square_matrix_print(Q);
      
      printf("full matrix is: \n");
      square_matrix_print(full_Q);
#endif // DEBUG_LEVEL_FULL
#endif // DEBUG
  }
  else {
    if (x[1] == 1) {
      /* M[,2] = M[,1] ; M[,1] = x~ */
      for (col = 2; col < QF_RANK; col++)
	Q[1][col] = Q[0][col];
      Q[1][1] = Q[0][0];
      for (col = 2; col < QF_RANK; col++)
	Q[0][col] = Qx[col];
      Q[0][1] = Qx[0];
      Q[0][0] = xQx;

      isometry_replace_vec(s, 1, x);
#ifdef DEBUG
      square_matrix_set(full_Q, Q);
      square_matrix_resymmetrize(full_Q);
      assert(isometry_is_isom(s, nbr_man->Q, full_Q));
      isometry_inv(s_inv, s);
      assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
    }
    else {
      if (x[2] == 1) {
	/* M[,3] = M[,1] ; M[,1] = x~ */
	for (col = 3; col < QF_RANK; col++)
	  Q[2][col] = Q[0][col];
	Q[2][2] = Q[0][0];
	Q[1][2] = Q[0][1];
	for (col = 3; col < QF_RANK; col++)
	  Q[0][col] = Qx[col];
	Q[0][2] = Qx[0];
	Q[0][1] = Qx[1];
	Q[0][0] = xQx;
	
	isometry_replace_vec(s, 2, x);
      
#ifdef DEBUG
	square_matrix_set(full_Q, Q);
	square_matrix_resymmetrize(full_Q);
	assert(isometry_is_isom(s, nbr_man->Q, full_Q));
	isometry_inv(s_inv, s);
	assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
      }
      else {
	if (x[3] == 1) {
	  /* M[,4] = M[,1] ; M[,1] = x~ */
	  Q[3][4] = Q[0][4];
	  Q[3][3] = Q[0][0];
	  Q[2][3] = Q[0][2];
	  Q[1][3] = Q[0][1];
	  Q[0][4] = Qx[4];
	  Q[0][3] = Qx[0];
	  Q[0][2] = Qx[2];
	  Q[0][1] = Qx[1];
	  Q[0][0] = xQx;

	  isometry_replace_vec(s, 3, x);
      
#ifdef DEBUG
	square_matrix_set(full_Q, Q);
	square_matrix_resymmetrize(full_Q);
	assert(isometry_is_isom(s, nbr_man->Q, full_Q));
	isometry_inv(s_inv, s);
	assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
	}
	else {
	  assert(x[4]==1);
	  /* M[,5] = M[,1] ; M[,1] = x~ */
	  Q[4][4] = Q[0][0];
	  Q[3][4] = Q[0][3];
	  Q[2][4] = Q[0][2];
	  Q[1][4] = Q[0][1];
	  Q[0][4] = Qx[0];
	  Q[0][3] = Qx[3];
	  Q[0][2] = Qx[2];
	  Q[0][1] = Qx[1];
	  Q[0][0] = xQx;

	  isometry_replace_vec(s, 4, x);
      
#ifdef DEBUG
	square_matrix_set(full_Q, Q);
	square_matrix_resymmetrize(full_Q);
	assert(isometry_is_isom(s, nbr_man->Q, full_Q));
	isometry_inv(s_inv, s);
	assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
	}
      }
    }
  }

  assert(!(Q[0][0] % nbr_man->p));
  
  if (!(Q[0][4] % nbr_man->p)) {
    if (Q[0][3] % nbr_man->p) {
      /* M[,4] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
      printf("swapping rows 4 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
      for (row = 0; row < 3; row++) {
	square_matrix_swap_elts(Q, row, 3, row, 4);
      }
      square_matrix_swap_elts(Q, 3, 3, 4, 4);
      isometry_swap_vecs(s, 3, 4);
#ifdef DEBUG
      square_matrix_set(full_Q, Q);
      square_matrix_resymmetrize(full_Q);
      assert(isometry_is_isom(s, nbr_man->Q, full_Q));
      isometry_inv(s_inv, s);
      assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
#ifdef DEBUG_LEVEL_FULL
      square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
    }
    else {
      if (Q[0][2] % nbr_man->p) {
	/* M[,3] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
	printf("swapping rows 3 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
	square_matrix_swap_elts(Q, 0, 2, 0, 4);
	square_matrix_swap_elts(Q, 1, 2, 1, 4);
	square_matrix_swap_elts(Q, 2, 2, 4, 4);
	square_matrix_swap_elts(Q, 2, 3, 3, 4);
	isometry_swap_vecs(s, 2, 4);
#ifdef DEBUG
	square_matrix_set(full_Q, Q);
	square_matrix_resymmetrize(full_Q);
	assert(isometry_is_isom(s, nbr_man->Q, full_Q));
	isometry_inv(s_inv, s);
	assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
#ifdef DEBUG_LEVEL_FULL
	square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
      }
      else {
	if (Q[0][1] % nbr_man->p) {
	  /* M[,2] <--> M[,5] */
#ifdef DEBUG_LEVEL_FULL
	  printf("swapping rows 2 and 5, now Q = \n");
#endif // DEBUG_LEVEL_FULL
	  square_matrix_swap_elts(Q, 0, 1, 0, 4);
	  square_matrix_swap_elts(Q, 1, 1, 4, 4);
	  square_matrix_swap_elts(Q, 1, 2, 2, 4);
	  square_matrix_swap_elts(Q, 1, 3, 3, 4);
	  isometry_swap_vecs(s, 1, 4);
#ifdef DEBUG
	  square_matrix_set(full_Q, Q);
	  square_matrix_resymmetrize(full_Q);
	  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
	  isometry_inv(s_inv, s);
	  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
#ifdef DEBUG_LEVEL_FULL
	  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
	}
	else {
	  /* a singular zero, return Q as-is*/
	  /* Re-symmetrize */
	  square_matrix_resymmetrize(Q);
	  return false;
	}
      }
    }
  }
  
  /* M[,1] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[1] <-- R[1] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] /= nbr_man->p;
  gcdext(Q[0][4], nbr_man->p, &q, &y);
  a1 = -Q[0][0]*q/2 % nbr_man->p;
  if (a1 < 0)
    a1 += nbr_man->p;

#ifdef DEBUG_LEVEL_FULL
  printf("q = %"  PRId64 ", a1 = %"  PRId64 "(y = %" PRId64 ")\n", q, a1, y);
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] += a1*(2*Q[0][4]+nbr_man->p*a1*Q[4][4]);
  Q[0][1] += nbr_man->p*a1*Q[1][4];
  Q[0][2] += nbr_man->p*a1*Q[2][4];
  Q[0][3] += nbr_man->p*a1*Q[3][4];
  Q[0][4] += nbr_man->p*a1*Q[4][4];

  isometry_add_vec(s, 0, nbr_man->p*a1, 4);

#ifdef DEBUG
  square_matrix_set(full_Q, Q);
  square_matrix_resymmetrize(full_Q);
  full_Q[0][0] *= nbr_man->p;
  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
  isometry_inv(s_inv, s);
  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,2] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[2] <-- R[2] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  a2 = -Q[0][1]*q % nbr_man->p;
  if (a2 < 0)
    a2 += nbr_man->p;
  
  Q[0][1] += a2*Q[0][4];
  Q[1][1] += a2*(2*Q[1][4]+a2*Q[4][4]);
  Q[1][2] += a2*Q[2][4];
  Q[1][3] += a2*Q[3][4];
  Q[1][4] += a2*Q[4][4];

  isometry_add_vec(s, 1, a2, 4);

#ifdef DEBUG
  square_matrix_set(full_Q, Q);
  square_matrix_resymmetrize(full_Q);
  full_Q[0][0] *= nbr_man->p;
  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
  isometry_inv(s_inv, s);
  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,3] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[3] <-- R[3] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
   
  a3 = -Q[0][2]*q % nbr_man->p;
  if (a3 < 0)
    a3 += nbr_man->p;
  
  Q[0][2] += a3*Q[0][4];
  Q[1][2] += a3*Q[1][4];
  Q[2][2] += a3*(2*Q[2][4]+a3*Q[4][4]);
  Q[2][3] += a3*Q[3][4];
  Q[2][4] += a3*Q[4][4];

  isometry_add_vec(s, 2, a3, 4);

#ifdef DEBUG
  square_matrix_set(full_Q, Q);
  square_matrix_resymmetrize(full_Q);
  full_Q[0][0] *= nbr_man->p;
  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
  isometry_inv(s_inv, s);
  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,4] += a*M[,5]  */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[4] <-- R[4] + aR[5], now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  a4 = -Q[0][3]*q % nbr_man->p;
  if (a4 < 0)
    a4 += nbr_man->p;
  
  Q[0][3] += a4*Q[0][4];
  Q[1][3] += a4*Q[1][4];
  Q[2][3] += a4*Q[2][4];
  Q[3][3] += a4*(2*Q[3][4]+a4*Q[4][4]);
  Q[3][4] += a4*Q[4][4];

  isometry_add_vec(s, 3, a4, 4);

#ifdef DEBUG
  square_matrix_set(full_Q, Q);
  square_matrix_resymmetrize(full_Q);
  full_Q[0][0] *= nbr_man->p;
  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
  isometry_inv(s_inv, s);
  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,1] /= p */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[1] <-- R[1]/p, now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[0][0] /= nbr_man->p;
  Q[0][1] /= nbr_man->p;
  Q[0][2] /= nbr_man->p;
  Q[0][3] /= nbr_man->p;

  isometry_vec_scalar_div(s,0,nbr_man->p);

#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* M[,5] *= p */
#ifdef DEBUG_LEVEL_FULL
  printf("Doing R[5] <-- R[5]/p, now Q = \n");
#endif // DEBUG_LEVEL_FULL
  
  Q[1][4] *= nbr_man->p;
  Q[2][4] *= nbr_man->p;
  Q[3][4] *= nbr_man->p;
  Q[4][4] *= (nbr_man->p)*(nbr_man->p);

  isometry_vec_scalar_mul(s,4,nbr_man->p);

#ifdef DEBUG
  square_matrix_set(full_Q, Q);
  square_matrix_resymmetrize(full_Q);
  assert(isometry_is_isom(s, nbr_man->Q, full_Q));
  isometry_inv(s_inv, s);
  assert(isometry_is_isom(s_inv, full_Q, nbr_man->Q));
#endif // DEBUG
  
#ifdef DEBUG_LEVEL_FULL
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL
  
  /* Re-symmetrize */
  square_matrix_resymmetrize(Q);

#ifdef DEBUG_LEVEL_FULL
  printf("Resulting neighbor is: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  assert(square_matrix_is_positive_definite(Q));
#endif // DEBUG
  
  return true;
};
