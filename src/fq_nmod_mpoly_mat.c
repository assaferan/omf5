/*************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_mpoly_mat.c
 *
 * Description: Matrices of multivariable polynomials over
 *              finite fields.
 *
 *************************************************************
 */

// Self dependencies

#include "fq_nmod_mpoly_mat.h"

// initialize a matrix
void fq_nmod_mpoly_mat_init(fq_nmod_mpoly_mat_t mat, slong rows, slong cols, const fq_nmod_mpoly_ctx_t R)
{
  slong i;

  if (rows != 0)
    mat->rows = (fq_nmod_mpoly_struct **) flint_malloc(rows * sizeof(fq_nmod_mpoly_struct *));
  else
    mat->rows = NULL;

  if (rows != 0 && cols != 0)
    {
      mat->entries = (fq_nmod_mpoly_struct *) flint_calloc(flint_mul_sizes(rows, cols), sizeof(fq_nmod_mpoly_struct));

      /* Set denominators */
      for (i = 0; i < rows * cols; i++)
	fq_nmod_mpoly_init(&(mat->entries[i]), R);

      for (i = 0; i < rows; i++)
	mat->rows[i] = mat->entries + i * cols;
    }
  else
    {
      mat->entries = NULL;
      if (rows != 0)
        {
	  for (i = 0; i < rows; i++)
	    mat->rows[i] = NULL;
        }
    }

  mat->r = rows;
  mat->c = cols;
}

// clear the memeory allocated for the matrix
void fq_nmod_mpoly_mat_clear(fq_nmod_mpoly_mat_t mat, const fq_nmod_mpoly_ctx_t R)
{
  if (mat->entries)
    {
      slong i;
      
      for (i = 0; i < mat->r * mat->c; i++)
	fq_nmod_mpoly_clear(mat->entries + i, R);

      flint_free(mat->entries);
      flint_free(mat->rows);
    } else if (mat->r != 0)
    flint_free(mat->rows);
}

// print the matrix
void fq_nmod_mpoly_mat_print(const fq_nmod_mpoly_mat_t mat, const char** var_names, const fq_nmod_mpoly_ctx_t R)
{
    slong i, j;

    flint_printf("<%wd x %wd matrix over R>\n", mat->r, mat->c);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < mat->c; j++)
        {
	  fq_nmod_mpoly_print_pretty(fq_nmod_mpoly_mat_entry(mat, i, j), var_names, R);
            if (j + 1 < mat->c)
                flint_printf(", ");
        }
        flint_printf("]\n");
    }
    flint_printf("\n");
}
