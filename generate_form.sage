import sys
from sage.all_cmdline import *
from sage.rings.integer import Integer

load("generate_dprt.sage")

def print_with_spaces(coeffs):
    for coeff in coeffs:
        print(coeff, end=' ')
    print()
dminus = Integer(sys.argv[1])
dplus = Integer(sys.argv[2])
qf = generate_form_dprt(dminus, dplus)
print_with_spaces(qf.coefficients())



