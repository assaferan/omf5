#ifndef __IO_H__
#define __IO_H__

#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

bool parse_matrix(Z64* Q_coeffs, const char* mat_str);
size_t read_genus(square_matrix_t** genus, const char* fname, size_t disc);
size_t read_genus_and_isom(square_matrix_t** p_genus,
			   isometry_t** p_isom,
			   const char* fname, size_t disc);

#endif // _IO_H__
