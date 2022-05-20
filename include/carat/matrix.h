#ifndef CARAT_MATRIX_H
#define CARAT_MATRIX_H

#include "typedef.h"

/*--------------------------------------------------------------*\
| FILE add_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *imat_add (matrix_TYP *L_mat, matrix_TYP *R_mat,
     int Lc, int Rc);
extern matrix_TYP *imat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat,
     int Lc, int Rc);
extern matrix_TYP *rmat_add (matrix_TYP *L_mat, matrix_TYP *R_mat,
     rational L_coeff, rational R_coeff);
extern matrix_TYP *rmat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat,
     rational L_coeff, rational R_coeff);
extern matrix_TYP *pmat_add (matrix_TYP *L_mat, matrix_TYP *R_mat,
     int L_coeff, int R_coeff);
extern matrix_TYP *pmat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat,
     int L_coeff, int R_coeff);
extern matrix_TYP *mat_add (matrix_TYP *L_mat, matrix_TYP *R_mat,
     rational L_coeff, rational R_coeff);
extern matrix_TYP *mat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat,
     rational L_coeff, rational R_coeff);

/*--------------------------------------------------------------*\
| FILE col_row_ops_mat.c 
\*--------------------------------------------------------------*/
extern void row_per(matrix_TYP *M, int i, int j);
extern void col_per(matrix_TYP *M, int i, int j);
extern void row_add(matrix_TYP *M, int i, int j, int fac);
extern void col_add(matrix_TYP *M, int i, int j, int fac);
extern void row_mul(matrix_TYP *M, int i, int fac);
extern void col_mul(matrix_TYP *M, int i, int fac);

/*--------------------------------------------------------------*\
| FILE comp_mat.c 
\*--------------------------------------------------------------*/
extern int cmp_mat(matrix_TYP *A, matrix_TYP *B);

/*--------------------------------------------------------------*\
| FILE divide_by_gcd.c 
\*--------------------------------------------------------------*/
extern int divide_by_gcd(matrix_TYP *A);

/*--------------------------------------------------------------*\
| FILE construct_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *init_mat(int rows, int cols, const char *option);
extern matrix_TYP *copy_mat( matrix_TYP *old );
extern void free_mat (matrix_TYP *mat);
extern void Check_mat(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE elt_div_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *elt_div(matrix_TYP *Mat);

/*--------------------------------------------------------------*\
| FILE find_max_entry_mat.c 
\*--------------------------------------------------------------*/
extern int find_max_entry(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE gauss_mat.c 
\*--------------------------------------------------------------*/
extern int tgauss(matrix_TYP *mat);
extern int row_gauss(matrix_TYP *M);
extern matrix_TYP *ggauss(matrix_TYP *mat);

/************************************************************************\
|  FILE: hnf_mat.c
\************************************************************************/
extern int row_hnf_mat_trf(matrix_TYP *Mat, matrix_TYP *T);
extern int row_hnf_mat_simultaneous(matrix_TYP *Mat, matrix_TYP *T);
extern int row_hnf_mat(matrix_TYP *Mat);
extern int col_hnf_mat_trf(matrix_TYP *Mat, matrix_TYP *T);
extern int col_hnf_mat_simultaneous(matrix_TYP *Mat, matrix_TYP *T);
extern int col_hnf_mat(matrix_TYP *Mat);

/*--------------------------------------------------------------*\
| FILE inv_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *pmat_inv (matrix_TYP *mat);
extern matrix_TYP *mat_inv(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE kernel_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP kernel_mat(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE kgv2rat_mat.c 
\*--------------------------------------------------------------*/
extern int kgv2rat( matrix_TYP *mat );
extern int rat2kgv( matrix_TYP *mat );

/*--------------------------------------------------------------*\
| FILE kron_mat.c
\*--------------------------------------------------------------*/
extern matrix_TYP *kron_mat(matrix_TYP *A, matrix_TYP *B);

/*--------------------------------------------------------------*\
| FILE modp_mat.c
\*--------------------------------------------------------------*/
extern void modp_mat(matrix_TYP *M, int prime);

/*--------------------------------------------------------------*\
| FILE mul_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *pmat_mul (matrix_TYP *L_mat, matrix_TYP * R_mat);
extern matrix_TYP *mat_mul (matrix_TYP *L_mat, matrix_TYP * R_mat);
extern matrix_TYP *mat_muleq(matrix_TYP *L_mat, matrix_TYP *R_mat);
extern matrix_TYP *mat_kon(matrix_TYP *L_mat, matrix_TYP *M_mat,
     matrix_TYP *R_mat);

/*--------------------------------------------------------------*\
| FILE null_mat.c 
\*--------------------------------------------------------------*/
extern int null_mat(matrix_TYP *mat);
extern int save_null_mat(matrix_TYP *mat);
extern int quick_null_mat(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE p_mat_det.c 
\*--------------------------------------------------------------*/
extern int p_mat_det(matrix_TYP *M, int prime);

/*--------------------------------------------------------------*\
| FILE p_gauss_mat.c 
\*--------------------------------------------------------------*/
extern int p_gauss (matrix_TYP *L_mat);

/*--------------------------------------------------------------*\
| FILE p_lse_solve.c 
\*--------------------------------------------------------------*/
extern matrix_TYP **p_lse_solve(matrix_TYP *A, matrix_TYP *B, int *anz, int p);

/*--------------------------------------------------------------*\
| FILE p_solve_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP **p_solve (int *anz, matrix_TYP **L_mat, matrix_TYP **R_mat,
                      int option);

/*--------------------------------------------------------------*\
| FILE real_mat.c 
\*--------------------------------------------------------------*/
extern void real_mat(matrix_TYP *mat, int rows, int cols);

/*--------------------------------------------------------------*\
| FILE red_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *mat_red(matrix_TYP *Mat);
extern void dec_mat(matrix_TYP *Mat, matrix_TYP *Trf);

/*--------------------------------------------------------------*\
| FILE scal_pr_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *scal_pr ( matrix_TYP *vectors, matrix_TYP *form,
     boolean truth );

/*--------------------------------------------------------------*\
| FILE solve_mat.c 
\*--------------------------------------------------------------*/
extern int Trf_gauss(matrix_TYP *M, matrix_TYP *Trf);
extern matrix_TYP *solve_mat(matrix_TYP *M);

/*--------------------------------------------------------------*\
| FILE tools_mat.c 
\*--------------------------------------------------------------*/
extern boolean iset_entry(matrix_TYP *mat, int r, int c, int v);
extern boolean rset_entry(matrix_TYP *mat, int r, int c, rational v);
extern void iscal_mul(matrix_TYP *mat, int v);
extern void rscal_mul(matrix_TYP *mat, rational v);
extern boolean kill_row(matrix_TYP *mat, int row);
extern boolean kill_col(matrix_TYP *mat, int col);
extern boolean ins_row(matrix_TYP *mat, int row);
extern boolean ins_col(matrix_TYP *mat,  int col);
extern boolean imul_row(matrix_TYP *mat, int row, int v);
extern boolean rmul_row(matrix_TYP *mat, int row, rational v);
extern boolean imul_col(matrix_TYP *mat, int col, int v);
extern boolean rmul_col( matrix_TYP *mat, int col, rational v);
extern boolean iadd_row(matrix_TYP *mat,int t_row, int d_row, int v);
extern boolean radd_row(matrix_TYP *mat, int t_row, int d_row, rational v);
extern boolean iadd_col(matrix_TYP *mat, int t_col, int d_col, int v);
extern boolean radd_col(matrix_TYP *mat, int t_col, int d_col, rational v);
extern void normal_rows(matrix_TYP *mat);
extern void normal_cols(matrix_TYP *mat);
extern matrix_TYP *mat_to_line(matrix_TYP **gen, int num);
extern matrix_TYP **line_to_mat(matrix_TYP *mat, int row, int col);
extern int normal_mat(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE tr_pose_mat.c 
\*--------------------------------------------------------------*/
extern matrix_TYP *tr_pose(matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE trace_mat.c 
\*--------------------------------------------------------------*/
extern int trace( matrix_TYP *mat);

/*--------------------------------------------------------------*\
| FILE unity_mat.c
\*--------------------------------------------------------------*/
extern matrix_TYP *einheitsmatrix( int n);

#endif
