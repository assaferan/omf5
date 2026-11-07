#ifndef __FLINT_COMPAT_H__
#define __FLINT_COMPAT_H__

/*****************************************************************************
 *
 * Package     : omf5
 * Filename    : flint_compat.h
 *
 * Description : Compatibility layer between flint 2.x (with antic) and
 *               flint >= 3.0, which absorbed antic and changed a few
 *               signatures.
 *
 *****************************************************************************
 */

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#if __FLINT_RELEASE >= 30000

// nf and nf_elem were merged from antic into flint
#include <flint/nf.h>
#include <flint/nf_elem.h>

// fq_nmod_ctx_prime returns a ulong
#define FQ_NMOD_CTX_PRIME_UI(F) (fq_nmod_ctx_prime(F))

// fq_nmod_mat_rref takes separate output and input (aliasing allowed)
#define FQ_NMOD_MAT_RREF(A, F) (fq_nmod_mat_rref((A), (A), (F)))

#else // flint 2.x

#include <antic/nf.h>
#include <antic/nf_elem.h>

// fq_nmod_ctx_prime returns a const fmpz *
#define FQ_NMOD_CTX_PRIME_UI(F) (fmpz_get_ui(fq_nmod_ctx_prime(F)))

// fq_nmod_mat_rref works in place
#define FQ_NMOD_MAT_RREF(A, F) (fq_nmod_mat_rref((A), (F)))

#endif // __FLINT_RELEASE

#endif // __FLINT_COMPAT_H__
