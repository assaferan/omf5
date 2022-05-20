#include "carat/symm.h"

#include "genus.h"
#include "matrix_tools.h"

/* initialize the neighbors data */
nbrs_data* q61_init(int p, int k)
{
  nbrs_data* dtm;
  matrix_TYP *Qx, *n_mat;
  
  int* x;
  int s, i, norm, n, v1, v2, v3, v4;
  int num_short[3];

  matrix_TYP* genus[8];

  /* Set of representatives for the quadratic forms */
  
  int genus_coeffs[8][15] =
    {
     {2,1,2,1,0,2,0,1,0,16,0,0,0,0,2},
     {2,1,2,0,0,2,1,1,1,6,-1,-1,0,-3,6},
     {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6},
     {2,-1,2,-1,1,2,1,-1,-1,6,-1,1,0,2,8},
     {2,1,2,1,0,4,-1,-1,0,4,0,0,1,0,4},
     {2,1,2,-1,0,2,0,0,0,6,-1,0,1,-1,6},
     {2,-1,2,0,0,2,1,0,-1,4,1,0,0,0,8},
     {2,0,2,0,0,2,0,0,1,4,0,1,1,2,6}
    };

  /* printf("in q61_init. Initializing genus representatives to:\n"); */
  
  for (i = 0; i < 8; i++) {
    genus[i] = init_sym_matrix(genus_coeffs[i]);
    /* print_mat(genus[i]); */
  }
  
  dtm = (nbrs_data*) malloc(sizeof(nbrs_data));

  /* Hardcode detect class */
  // !! This depends on the size of our hash
  // In general replace by some hash
  dtm->th61 = (int*) malloc(20*sizeof(int));

  /* printf("Computing theta series...\n"); */
  for (i = 0; i < 8; i++) {
    /* printf("For the lattice: \n"); */
    /* print_mat(genus[i]); */
    for (norm = 2; norm <= 6; norm += 2) {
      /* !! TODO !! Make sure this does what I think it does  */
      /* short_vectors(Q61_mat_i, norm, norm, 0, 1, &(num_short[norm/2])); */
      short_vectors(genus[i], norm, norm, 0, 1, &(num_short[norm/2 - 1]));
    }
    /* printf("We get: %d q + %d q^2 + %d q^3 + O(q^4).\n", */
    /* 	   num_short[0], num_short[1], num_short[2]); */
    // Gonzalo chose this one since it is enough to distinuish the
    // representatives
    s = num_short[0] - num_short[1] + num_short[2];

    /* printf("s = %d\n", s); */
    dtm->th61[s] = i;
  }

  dtm->Q = init_sym_matrix(genus_coeffs[k]);
  dtm->v = init_mat(1,5,"");
  x = dtm->v->array.SZ[0];

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
	  Qx = mat_mul(dtm->v, dtm->Q);
	  n_mat = mat_mul(dtm->v, tr_pose(Qx));
	  n = n_mat->array.SZ[0][0] / 2 % p;
	  if (!n)
	    return dtm;
	}
      }
    }
  }

  printf("Didn't find an isotropic vector!?\n");
  return dtm;
}


/* identify the genus representative of Q (returns the index) */
int q61_id(matrix_TYP* Q, int* th61)
{
  int norm, x;
  int num_short[3];
  
  for (norm = 2; norm <= 6; norm += 2) {
    // !! TODO !! Make sure this does what I think it does 
    short_vectors(Q, norm, norm, 0, 1, &(num_short[norm/2 - 1]));
  }
  x = num_short[0] - num_short[1] + num_short[2];

  /* printf("id = %d\n", th61[x]); */
  
  return th61[x];
}
