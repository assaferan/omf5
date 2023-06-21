/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : io.h
 *
 * Description: Functions for handling input and output
 *
 *****************************************************************
 */

#ifndef __IO_H__
#define __IO_H__

// Self dependencies

#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

/***************************************************************************
 *
 * Function: parse_matrix
 *
 * Description: Given a string of entries describing a symmetric 5x5 matrix,
 *              return an array of 64-bit integers containing the entries.
 *              It is assumed that there are exactly 15 entries,
 *              describing the lower triangular part of the matrix.
 *              Assumes the reuired target memory has been allocated.
 *
 * Arguments:
 *     + mat_str (const char*) - the input string
 *
 * Returns:
 *     + (bool) - True if the data was formatted as expected
 *     + Q_coeffs (Z64*) - an array containing the integers from the input
 *
 ***************************************************************************
 */

bool parse_matrix(Z64* Q_coeffs, const char* mat_str);

/***************************************************************************
 *
 * Function: read_genus
 *
 * Description: Reads a list of genus representatives from a file.
 *              The file is assumed to be in the form of a python array,
 *              where the N-th entry is an array of arrays of integers,
 *              each corresponding to the gram matrix of a genus
 *              representative for the space of discrminant N.
 *
 * Arguments:
 *     + fname (const char*) - the name of the file containing the genus
 *     + disc (size_t) - the discriminant of the genus to read
 *
 * Returns:
 *     + (size_t) - number of genus representatives read.
 *     + genus (square_matrix_t**) - a pointer to an array of the gram
 *                                   matrices of genus representatives
 *
 ***************************************************************************
 */

size_t read_genus(square_matrix_t** genus, const char* fname, size_t disc);

/***************************************************************************
 *
 * Function: read_genus_and_isom
 *
 * Description: Reads a list of genus representatives from a file.
 *              The file is assumed to be in the form of a python array,
 *              where the N-th entry is an array of pairs of matrices,
 *              each corresponding to the gram matrix of a genus
 *              representative for the space of discrminant N,
 *              and an isometry (over Q) from the first representative. 
 *
 * Arguments:
 *     + fname (const char*) - the name of the file containing the genus
 *     + disc (size_t) - the discriminant of the genus to read
 *
 * Returns:
 *     + (size_t) - number of genus representatives read.
 *     + p_genus (square_matrix_t**) - a pointer to an array of the gram
 *                                     matrices of genus representatives
 *     + p_isom (isometry_t**) - a pointer to an array of the isometries
 *
 ***************************************************************************
 */

size_t read_genus_and_isom(square_matrix_t** p_genus,
			   isometry_t** p_isom,
			   const char* fname, size_t disc);

#endif // _IO_H__
