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

## Installation

To build as a C library on a standard Linux/Mac OS system with autotools:

    ./autogen.sh
    ./configure --prefix=DIR
    make
    make install

## Usage



## Contribution

If you want to help develop this project, please create your own fork on Github and submit a pull request. I will do my best to integrate any additional useful features as necessary. Alternatively, submit a patch to me via email at assaferan@gmail.com.
