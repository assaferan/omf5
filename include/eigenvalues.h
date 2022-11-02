#ifndef __EIGENVALUES_H__
#define __EIGENVALUES_H__

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include <carat/typedef.h>

#include <flint/fmpz.h>

#include "genus.h"
#include "typedefs.h"

typedef struct {
  nf_t* nfs;
  nf_elem_t* eigenvals; // right now only stores the eigenvalue of the Hecke operator used to construct it. rethink that
  nf_elem_t** eigenvecs; 
  slong num;
  slong dim;
  bool* is_lift;
} eigenvalues;

typedef eigenvalues eigenvalues_t[1];

typedef struct {
  eigenvalues_t* evs;
  slong num;
} all_eigen;

void eigenvalues_init(eigenvalues_t evs, slong num, slong dim);

void eigenvalues_init_set_mat(eigenvalues_t evs, matrix_TYP* mat);

void eigenvalues_set_lifts(eigenvalues_t evs, slong prec, slong c, const genus_t genus);

void eigenvalues_clear(eigenvalues_t evs);

#endif // __EIGENVALUES_H__
