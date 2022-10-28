#ifndef __IO_H__
#define __IO_H__

#include "square_matrix.h"
#include "typedefs.h"

bool parse_matrix(const char* mat_str, int* Q_coeffs);
size_t read_genus(square_matrix_t** genus, const char* fname, size_t disc);

#endif // _IO_H__
