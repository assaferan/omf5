#! /usr/local/bin/sage -python

import sys
from parse_omf import parse_omf5

N = int(sys.argv[1])
print(N)
parse_omf5(3,0,N,"data_nl_200/")
    
