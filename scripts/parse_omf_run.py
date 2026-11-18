#! /usr/local/bin/sage -python

import sys
from parse_omf import parse_omf5, DEFAULT_GENUS

# parse_omf_run.py N [GENUS]: GENUS is the genus file of the run (omf5's -genus),
# by default data/qf5db.sage of this omf5 checkout
N = int(sys.argv[1])
genus = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_GENUS
print(N)
parse_omf5(3,0,N,"data_nl_200/",genus=genus)
    
