omf5 - Orthogonal Modular Forms of Rank 5
=========================================

A library for computing Hecke matrices for positive definite rational quinary quadratic forms.

The code in this repository is based on the work of many, among them:

1. Gonzalo Tornaria - His fast implementation for discriminant 61 is the basis for this code.
2. Jeffery P. Hein - his repository of ternary_birch contains the case of rank 3 and we have widely used ideas from there.
3. Matthew Greenberg and John Voight - their magma code for p-neighbors is the algorithm on which this code is based. (https://github.com/assaferan/ModFrmAlg)

We also use the CARAT package for working with lattices

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contribution)

## Requirements

- [autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html)
- [gcc](https://gcc.gnu.org/)

* This code has only been tested so far on:
- macOS Monterey 12.3.1, compiled with gcc 4.2.1 on Intel i9 8-Core Processor.

The main (and develop) branch require the following libraries:

- [CARAT](https://github.com/lbfm-rwth/carat)
- [flint](https://github.com/wbhart/flint2)
- [antic](https://github.com/wbhart/antic)

Installation of all the requirements on macOS can be done as follows, from a folder containing the above 3 libraries.

    brew update && brew install autoconf automake libtool
    cd flint2 && ./configure && make && make check && make install && cd ..
    cd antic && ./configure && make && make install && make check && cd ..
    cd carat && ./autogen.sh && ./configure && make && make install && cp libfunctions.a /usr/local/lib && cd ..

## Installation

To build as a C library on a standard Linux/Mac OS system with autotools:

    ./autogen.sh
    ./configure --prefix=DIR
    make
    make install

Or all at once:

   ./autogen.sh && ./configure && make && make install

## Usage

The executable is src/omf5.
It runs with command-line arguments as follows

src/omf5 [-tests] [-quad=Q] [-format=f] [-prec=L] [-form_idx=idx]

where the  arguments are
[Q] - the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f, currently one of two.
[f] - either 'GG' or 'A', the former for the format in Rama-Toranria webpage, the latter for the Magma format, as in the ModFrmAlg package.
[L] - the preicision up to which to compute the hecke eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L).
[idx]- the index of the form in the decomposition to eigenvectors.
If either L or i is not supplied, only decomposes the space, and finds Hecke eigenvectors.
If the flag -tests is supplied, additionally runs standard tests.

## Contribution

If you want to help develop this project, please create your own fork on Github and submit a pull request. I will do my best to integrate any additional useful features as necessary. Alternatively, submit a patch to me via email at assaferan@gmail.com.
