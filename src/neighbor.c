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
