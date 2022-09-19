#include "fq_nmod_poly_mat.h"

void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t F)
{
  slong i;

  if (rows != 0)
    mat->rows = (fq_nmod_poly_struct **) flint_malloc(rows * sizeof(fq_nmod_poly_struct *));
  else
    mat->rows = NULL;

  if (rows != 0 && cols != 0)
    {
      mat->entries = (fq_nmod_poly_struct *) flint_calloc(flint_mul_sizes(rows, cols), sizeof(fq_nmod_poly_struct));

      /* Set denominators */
      for (i = 0; i < rows * cols; i++)
	fq_nmod_poly_init(&(mat->entries[i]), F);

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

void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t F)
{
  if (mat->entries)
    {
      slong i;
      
      for (i = 0; i < mat->r * mat->c; i++)
	fq_nmod_poly_clear(mat->entries + i, F);

      flint_free(mat->entries);
      flint_free(mat->rows);
    } else if (mat->r != 0)
    flint_free(mat->rows);
}
