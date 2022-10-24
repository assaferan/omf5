#include "isometry.h"
#include "matrix_tools.h"

void isometry_init(isometry_t isom)
{
  // by default we initialize isometries to one
  square_matrix_one(isom->s);
  square_matrix_one(isom->s_inv);
  isom->denom = 1;
  isom->inv_denom = 1;

  return;
}

void isometry_init_set(isometry_t dest, const isometry_t src)
{
  square_matrix_set(dest->s, src->s);
  square_matrix_set(dest->s_inv, src->s_inv);
  dest->denom = src->denom;
  dest->inv_denom = src->inv_denom;

  return;
}

void isometry_init_set_square_matrix(isometry_t isom, const square_matrix_t s, int denom)
{
  // maybe we want to save time copying - check later
  square_matrix_set(isom->s,s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, s, isom->denom);

  return;
}

void isometry_init_set_fmpz_mat(isometry_t isom, const fmpz_mat_t s, int denom)
{

  square_matrix_set_fmpz_mat(isom->s, s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, s, isom->denom);

  return;
}

void isometry_mul(isometry_t prod, const isometry_t s1, const isometry_t s2)
{

  square_matrix_mul(prod->s, s1->s, s2->s);
  square_matrix_mul(prod->s_inv, s2->s_inv, s1->s_inv);
  prod->denom = (s1->denom) * (s2->denom);
  prod->inv_denom = (s1->inv_denom) * (s2->inv_denom);

  return;
}

void isometry_mul_mat_left(isometry_t prod, const square_matrix_t s1, int s1_denom, const isometry_t s2)
{
  
  square_matrix_mul(prod->s, s1, s2->s);
  prod->denom = s1_denom * (s2->denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);

  return;
}

void isometry_mul_mat_right(isometry_t prod, const isometry_t s1, const square_matrix_t s2, int s2_denom)
{
  
  square_matrix_mul(prod->s, s1->s, s2);
  prod->denom = s1->denom * (s2_denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);
  
  return;
}

void isometry_inv(isometry_t inv, const isometry_t isom)
{
  
  square_matrix_set(inv->s, isom->s_inv);
  square_matrix_set(inv->s_inv, isom->s);
  inv->denom = isom->inv_denom;
  inv->inv_denom = isom->denom;

  return;
}

void isometry_clear(isometry_t isom)
{
  // did not allocate any memory
  return;
}

// !! TODO - we can make this operation faster by doing it all at once
void isometry_transform_gram(square_matrix_t gtQg, const isometry_t g, const square_matrix_t Q)
{
  square_matrix_t gt, Qg;
  
  square_matrix_transpose(gt, g->s);
  square_matrix_mul(Qg, Q, g->s);
  square_matrix_mul(gtQg, gt, Qg);
  square_matrix_div_scalar(gtQg, gtQg, (g->denom)*(g->denom));

  return;
}

bool isometry_is_isom(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  square_matrix_t gram;

#ifdef DEBUG_LEVEL_FULL
  printf("s = \n");
  square_matrix_print(s->s);
  printf("\n");
  printf("q1 = \n");
  square_matrix_print(q1);
  printf("\n");
  printf("q2 = \n");
  square_matrix_print(q2);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  isometry_transform_gram(gram, s, q1);
  
  return square_matrix_is_equal(gram,q2);
}

/* // !! TODO - switch to use transform */
/* bool is_isometry(matrix_TYP* s, matrix_TYP* q1, matrix_TYP* q2, int denom) */
/* { */
/*   matrix_TYP *s_t, *scaled_q2; */
/*   matrix_TYP *q1_s, *s_t_q1_s; */
/*   // matrix_TYP *q1_s_t, *s_q1_s_t; */
/*   bool ret; */

/* #ifdef DEBUG_LEVEL_FULL */
/*   printf("s = \n"); */
/*   print_mat(s); */
/*   printf("\n"); */
/*   printf("q1 = \n"); */
/*   print_mat(q1); */
/*   printf("\n"); */
/*   printf("q2 = \n"); */
/*   print_mat(q2); */
/*   printf("\n"); */
/*   printf("denom = %d\n", denom); */
/* #endif // DEBUG_LEVEL_FULL */
  
/*   s_t = tr_pose(s); */
  
/*   q1_s = mat_mul(q1, s); */
/*   s_t_q1_s = mat_mul(s_t, q1_s); */

/*   // q1_s_t = mat_mul(q1, s_t); */
/*   // s_q1_s_t = mat_mul(s, q1_s_t); */
  
/* #ifdef DEBUG_LEVEL_FULL */
/*   printf("st_q1_s = \n"); */
/*   print_mat(s_t_q1_s); */
/*   // printf("s_q1_st = \n"); */
/*   // print_mat(s_q1_s_t); */
/*   printf("\n"); */
/* #endif // DEBUG_LEVEL_FULL */
  
/*   scaled_q2 = copy_mat(q2); */
/*   iscal_mul(scaled_q2, denom*denom); */

/*   ret = (cmp_mat(s_t_q1_s,scaled_q2) == 0); */

/*   // ret = (cmp_mat(s_q1_s_t,scaled_q2) == 0); */
  
/*   free(scaled_q2); */
/*   free_mat(s_t); */
/*   free_mat(s_t_q1_s); */
/*   free_mat(q1_s); */
/*   // free_mat(s_q1_s_t); */
/*   // free_mat(q1_s_t); */
/*   return ret; */
/* } */
