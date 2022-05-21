#include "arith.h"
#include "matrix_tools.h"
#include "neighbor.h"

/* Compute one p-neighbour for Q_orig corresponding to vector x 
 * On error, return NULL.
*/
matrix_TYP* q61_nb(matrix_TYP* Q_orig, int p, matrix_TYP* x_mat)
{
  matrix_TYP *Q_mat, *Qx_mat, *xQx_mat;
  int q, y, row, col, *x, **Q, *Qx, xQx;
  int a1, a2, a3, a4;

  /* printf("Computing p-neighbor for Q_orig: \n"); */
  /* print_mat(Q_orig); */

  /* printf("Corresponding to isotropic vector "); */
  /* print_mat(x_mat); */
  
  Q_mat = copy_mat(Q_orig);
  Qx_mat = mat_mul(x_mat, Q_mat);
  /* From now until the end, Q is computed only upper triangular */

  /* printf("Qx = \n"); */
  /* print_mat(Qx_mat); */
  
  x = x_mat->array.SZ[0];
  Q = Q_mat->array.SZ;
  Qx = Qx_mat->array.SZ[0];

  xQx_mat = mat_mul(Qx_mat, tr_pose(x_mat));

  xQx = xQx_mat->array.SZ[0][0];
  
  if (x[0] == 1) {
      /* M[,1] = x~ */
      
      Q[0][0] = xQx;
      for (col = 1; col < 5; col++)
	Q[0][col] = Qx[col];

      /* printf("put x as first vector, now Q = \n"); */
      /* print_mat(Q_mat); */
  }
  else {
    if (x[1] == 1) {
      /* M[,2] = M[,1] ; M[,1] = x~ */
      for (col = 2; col < 5; col++)
	Q[1][col] = Q[0][col];
      Q[1][1] = Q[0][0];
      for (col = 2; col < 5; col++)
	Q[0][col] = Qx[col];
      Q[0][1] = Qx[0];
      Q[0][0] = xQx;
    }
    else {
      if (x[2] == 1) {
	/* M[,3] = M[,1] ; M[,1] = x~ */
	for (col = 3; col < 5; col++)
	  Q[2][col] = Q[0][col];
	Q[2][2] = Q[0][0];
	Q[1][2] = Q[0][1];
	for (col = 3; col < 5; col++)
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
	  if (x[4] == 1) {
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
	  }
	  else {
	    printf("BAD\n");
	    return NULL;
	  }
	}
      }
    }
  }
 
  if (Q[0][0] % p) {
    printf("NOT 0 mod p\n");
    return NULL;
  }
  
  if (!(Q[0][4] % p)) {
    if (Q[0][3] % p) {
      /* M[,4] <--> M[,5] */
      /* printf("swapping rows 4 and 5, now Q = \n"); */
      for (row = 0; row < 3; row++) {
	swap(Q, row, 3, row, 4);
      }
      swap(Q, 3, 3, 4, 4);

      /* print_mat(Q_mat); */
    }
    else {
      if (Q[0][2] % p) {
	/* M[,3] <--> M[,5] */
	/* printf("swapping rows 3 and 5, now Q = \n"); */
	swap(Q, 0, 2, 0, 4);
	swap(Q, 1, 2, 1, 4);
	swap(Q, 2, 2, 4, 4);
	swap(Q, 2, 3, 3, 4);

	/* print_mat(Q_mat); */
      }
      else {
	if (Q[0][1] % p) {
	  /* M[,2] <--> M[,5] */
	  /* printf("swapping rows 2 and 5, now Q = \n"); */
	  swap(Q, 0, 1, 0, 4);
	  swap(Q, 1, 1, 4, 4);
	  swap(Q, 1, 2, 2, 4);
	  swap(Q, 1, 3, 3, 4);
	  
	  /* print_mat(Q_mat); */
	}
	else {
	  /* a singular zero, return Q as-is*/
	  /* Re-symmetrize */
	  resymmetrize(Q);
	  printf("???\n");
	  return Q_mat;
	}
      }
    }
  }
  
  /* M[,1] += a*M[,5]  */
  /* printf("Doing R[1] <-- R[1] + aR[5], now Q = \n"); */
  
  Q[0][0] /= p;
  gcdext(Q[0][4], p, &q, &y);
  a1 = -Q[0][0]*q/2 % p;
  if (a1 < 0)
    a1 += p;

  /* printf("q = %d, a1 = %d (y = %d)\n", q, a1, y); */
  
  Q[0][0] += a1*(2*Q[0][4]+p*a1*Q[4][4]);
  Q[0][1] += p*a1*Q[1][4];
  Q[0][2] += p*a1*Q[2][4];
  Q[0][3] += p*a1*Q[3][4];
  Q[0][4] += p*a1*Q[4][4];

  /* print_mat(Q_mat); */
  
  /* M[,2] += a*M[,5]  */
  /* printf("Doing R[2] <-- R[2] + aR[5], now Q = \n"); */
  
  a2 = -Q[0][1]*q % p;
  if (a2 < 0)
    a2 += p;
  
  Q[0][1] += a2*Q[0][4];
  Q[1][1] += a2*(2*Q[1][4]+a2*Q[4][4]);
  Q[1][2] += a2*Q[2][4];
  Q[1][3] += a2*Q[3][4];
  Q[1][4] += a2*Q[4][4];

  /* print_mat(Q_mat); */
  
  /* M[,3] += a*M[,5]  */
  /* printf("Doing R[3] <-- R[3] + aR[5], now Q = \n"); */
   
  a3 = -Q[0][2]*q % p;
  if (a3 < 0)
    a3 += p;
  
  Q[0][2] += a3*Q[0][4];
  Q[1][2] += a3*Q[1][4];
  Q[2][2] += a3*(2*Q[2][4]+a3*Q[4][4]);
  Q[2][3] += a3*Q[3][4];
  Q[2][4] += a3*Q[4][4];

  /* print_mat(Q_mat); */
  
  /* M[,4] += a*M[,5]  */
  /* printf("Doing R[4] <-- R[4] + aR[5], now Q = \n"); */
  
  a4 = -Q[0][3]*q % p;
  if (a4 < 0)
    a4 += p;
  
  Q[0][3] += a4*Q[0][4];
  Q[1][3] += a4*Q[1][4];
  Q[2][3] += a4*Q[2][4];
  Q[3][3] += a4*(2*Q[3][4]+a4*Q[4][4]);
  Q[3][4] += a4*Q[4][4];

  /* print_mat(Q_mat); */
  
  /* M[,1] /= p */
  /* printf("Doing R[1] <-- R[1]/p, now Q = \n"); */
  
  Q[0][0] /= p;
  Q[0][1] /= p;
  Q[0][2] /= p;
  Q[0][3] /= p;

  /* print_mat(Q_mat); */
  
  /* M[,5] *= p */
  /* printf("Doing R[5] <-- R[5]/p, now Q = \n"); */
  
  Q[1][4] *= p;
  Q[2][4] *= p;
  Q[3][4] *= p;
  Q[4][4] *= p*p;

  /* print_mat(Q_mat); */
  
  /* Re-symmetrize */
  resymmetrize(Q);

  /* printf("Resulting neighbor is: \n"); */
  /* print_mat(Q_mat); */

  /* is_definite = definite_test(Q_mat); */

  /* if (!is_definite) */
  /*   printf("Neighbor is not positive definite!!!\n"); */
  
  return Q_mat;
};

/* find an isotropic vector for Q mod p */
/* return a row vector 1x5 */
matrix_TYP* get_isotropic_vector(matrix_TYP* Q, int p)
{
  matrix_TYP *v, *Qx, *n_mat;
  int* x;
  int v1, v2, v3, v4, n;
  
  v = init_mat(1,5,"");
  x = v->array.SZ[0];

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
	  Qx = mat_mul(v, Q);
	  n_mat = mat_mul(v, tr_pose(Qx));
	  n = n_mat->array.SZ[0][0] / 2 % p;
	  if (!n)
	    return v;
	}
      }
    }
  }
  
  printf("Didn't find an isotropic vector!?\n");
  return v;
}

/* update the pivot vector v */
void update_pivot(int* v, int p, int i)
{
  int pivot, lsb;
  
  pivot = 0;
  while (v[pivot] == 0) {
    pivot++;
  }

  if (pivot == 0)
    printf("Error! Got a vector with pivot 0! shouldn't have gotten here\n");
  
  if (pivot == 5) {
    printf("Error! Got the zero vector!\n");
  }
  if (v[pivot] != 1)
    printf("Error! value at pivot is not 1!\n");

  lsb = 4;
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
    if (pivot < 4)
      v[pivot+1] = 0;
    v[pivot] = i;
    v[pivot-1] = 1;
  }

  return;
}

/* get isotropic vector , correposnding to the pivot vector w. */
/* If none is found returns NULL */
matrix_TYP* get_next_isotropic_vector(matrix_TYP*Q, int p,
				      matrix_TYP* v, matrix_TYP* w_mat,
				      matrix_TYP* b, int* iso_j)
{
  matrix_TYP *t_mat, *n_mat, *iso_vec;
  int j, t, n, n_inv, dummy;


  if (*iso_j) {
    iso_vec = imat_add(v,w_mat,1,*iso_j);
    (*iso_j)++;
    return iso_vec;
  }
  
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
      iso_vec = imat_add(v, w_mat, 1, -t*n_inv);
      /* printf("tmp = "); */
      /* print_mat(tmp); */
      /* This differs from the gp script */
      /* modp_mat(tmp, p); */
      for (j = 0; j < 5; j++)
	iso_vec->array.SZ[0][j] %= p;
      /* printf("n*v-t*w = "); */
      /* print_mat(iso_vec); */
      return iso_vec;
    }
  }
  else {
    if (*iso_j == 0)
      iso_vec = w_mat;
    else
      iso_vec = imat_add(v,w_mat,1,*iso_j);
    /* printf("v+%d*w = ", *iso_j); */
    /* print_mat(iso_vec); */
    if (!t)
      (*iso_j)++;
    return iso_vec;
  }

  return NULL;
}

void init_nbr_process(neighbor_manager* nbr_man, matrix_TYP* Q, int p)
{ 
  nbr_man->Q = Q;
  nbr_man->p = p;
  nbr_man->v = get_isotropic_vector(Q,p);
  nbr_man->b = mat_mul(nbr_man->v, Q);
  modp_mat(nbr_man->b,p);
  nbr_man->w = init_mat(1,5, "");
  nbr_man->i = 0;
  nbr_man->iso_j = 0;
  nbr_man->iso_vec = nbr_man->v;

  return;
}

void advance_nbr_process(neighbor_manager* nbr_man)
{
  int* w = nbr_man->w->array.SZ[0];
  
  // This is the first iteration
  if (nbr_man->iso_vec == nbr_man->v) {
    w[0] = w[1] = w[2] = w[3] = 0;
    w[4] = 1;
  }
  else {
    
    // finished looping over v+jw
    if (nbr_man->iso_j == nbr_man->p) {
      nbr_man->iso_j = 0; 
    }
    
    // update to next pivot
    if (nbr_man->iso_j == 0) {
      update_pivot(w, nbr_man->p, nbr_man->i);
    }
  }

  nbr_man->iso_vec = NULL;
  while ((nbr_man->iso_vec == NULL) && (nbr_man->i < nbr_man->p) ) {
    
    while ((w[0] == 0) && (nbr_man->iso_vec == NULL)) {
      nbr_man->iso_vec = get_next_isotropic_vector(nbr_man->Q, nbr_man->p,
						   nbr_man->v, nbr_man->w,
						   nbr_man->b,
						   &(nbr_man->iso_j));

      if (nbr_man->iso_vec == NULL) {
	// finished looping over v+jw
	if (nbr_man->iso_j == nbr_man->p) {
	  nbr_man->iso_j = 0; 
	}
	
	// update to next pivot
	if (nbr_man->iso_j == 0) {
	  update_pivot(w, nbr_man->p, nbr_man->i);
	}
      }
    }

    if (w[0] != 0) {
      (nbr_man->i)++;
      w[0] = w[1] = w[2] = 0;
      w[3] = 1;
      w[4] = nbr_man->i;
    }
  }

  return;
}
