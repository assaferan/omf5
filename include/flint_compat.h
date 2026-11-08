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

// nf and nf_elem were merged from antic into flint in 3.0
#if __FLINT_RELEASE >= 30000
#include <flint/nf.h>
#include <flint/nf_elem.h>
#else
#include <antic/nf.h>
#include <antic/nf_elem.h>
#endif

// fq_nmod_ctx_prime returns a const fmpz * in flint 2.x and early 3.x,
// and a ulong in later versions. Dispatch on the type actually returned,
// so that only the matching conversion is compiled.
static inline ulong _omf5_prime_ui_from_ui(ulong p) { return p; }
static inline ulong _omf5_prime_ui_from_fmpz(const fmpz * p) { return fmpz_get_ui(p); }

#define FQ_NMOD_CTX_PRIME_UI(F)                          \
  (_Generic((fq_nmod_ctx_prime(F)),                      \
            ulong: _omf5_prime_ui_from_ui,               \
            default: _omf5_prime_ui_from_fmpz)(fq_nmod_ctx_prime(F)))

// fq_nmod_mat_rref takes separate output and input matrices since flint 3.1
// (aliasing allowed), and works in place before that.
#if __FLINT_RELEASE >= 30100
#define FQ_NMOD_MAT_RREF(A, F) (fq_nmod_mat_rref((A), (A), (F)))
#else
#define FQ_NMOD_MAT_RREF(A, F) (fq_nmod_mat_rref((A), (F)))
#endif

#endif // __FLINT_COMPAT_H__
