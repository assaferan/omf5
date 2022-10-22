#include "isometry.h"
#include "matrix_tools.h"

void isometry_init(isometry_t isom, matrix_TYP* s, slong denom)
{
  // maybe we want to save time copying - check later
  isom->s = copy_mat(s);
  isom->s->kgv *= denom;
  isom->s_inv = mat_inv(s);

  return;
}

void isometry_init_fmpz_mat(isometry_t isom, fmpz_mat_t s, slong denom)
{
  
  matrix_TYP_init_set_fmpz_mat(&(isom->s), s);
  isom->s->kgv *= denom;
  isom->s_inv = mat_inv(isom->s);

  return;
}

void isometry_mul(isometry_t prod, isometry_t s1, isometry_t s2)
{

  prod->s = mat_mul(s1->s, s2->s);
  prod->s_inv = mat_mul(s2->s_inv, s1->s_inv);

  return;
}

void isometry_mul_mat_left(isometry_t prod, matrix_TYP* s1, isometry_t s2)
{

  prod->s = mat_mul(s1, s2->s);
  prod->s_inv = mat_mul(s2->s_inv, mat_inv(s1));

  return;
}

void isometry_mul_mat_right(isometry_t prod, isometry_t s1, matrix_TYP* s2)
{

  prod->s = mat_mul(s1->s, s2);
  prod->s_inv = mat_mul(mat_inv(s2), s1->s_inv);

  return;
}


void isometry_inv(isometry_t inv, isometry_t isom)
{
  inv->s = copy_mat(isom->s_inv);
  inv->s_inv = copy_mat(isom->s);

  return;
}

void isometry_clear(isometry_t isom)
{
  free_mat(isom->s);
  free_mat(isom->s_inv);

  return;
}

// !! TODO - this is leaky. Should fix that, but make sure timing is not severly comrpomised
matrix_TYP* isometry_transform_gram(isometry_t g, matrix_TYP* Q)
{
  return mat_mul(tr_pose(g->s), mat_mul(Q, g->s));
}

bool isometry_is_isom(isometry_t s, matrix_TYP* q1, matrix_TYP* q2)
{
  bool ret;

#ifdef DEBUG_LEVEL_FULL
  printf("s = \n");
  print_mat(s->s);
  printf("\n");
  printf("q1 = \n");
  print_mat(q1);
  printf("\n");
  printf("q2 = \n");
  print_mat(q2);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  ret = (cmp_mat(isometry_transform_gram(s,q1),q2) == 0);
  
  return ret;
}

// !! TODO - switch to use transform
bool is_isometry(matrix_TYP* s, matrix_TYP* q1, matrix_TYP* q2, int denom)
{
  matrix_TYP *s_t, *scaled_q2;
  matrix_TYP *q1_s, *s_t_q1_s;
  // matrix_TYP *q1_s_t, *s_q1_s_t;
  bool ret;

#ifdef DEBUG_LEVEL_FULL
  printf("s = \n");
  print_mat(s);
  printf("\n");
  printf("q1 = \n");
  print_mat(q1);
  printf("\n");
  printf("q2 = \n");
  print_mat(q2);
  printf("\n");
  printf("denom = %d\n", denom);
#endif // DEBUG_LEVEL_FULL
  
  s_t = tr_pose(s);
  
  q1_s = mat_mul(q1, s);
  s_t_q1_s = mat_mul(s_t, q1_s);

  // q1_s_t = mat_mul(q1, s_t);
  // s_q1_s_t = mat_mul(s, q1_s_t);
  
#ifdef DEBUG_LEVEL_FULL
  printf("st_q1_s = \n");
  print_mat(s_t_q1_s);
  // printf("s_q1_st = \n");
  // print_mat(s_q1_s_t);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  scaled_q2 = copy_mat(q2);
  iscal_mul(scaled_q2, denom*denom);

  ret = (cmp_mat(s_t_q1_s,scaled_q2) == 0);

  // ret = (cmp_mat(s_q1_s_t,scaled_q2) == 0);
  
  free(scaled_q2);
  free_mat(s_t);
  free_mat(s_t_q1_s);
  free_mat(q1_s);
  // free_mat(s_q1_s_t);
  // free_mat(q1_s_t);
  return ret;
}
