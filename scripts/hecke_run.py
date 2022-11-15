#! /usr/local/bin/sage -python

import os
from sage.all import prime_range

os.system("export LD_LIBRARY_PATH=/home/assaferan/lib")

MAX_DISC = 1000
MAX_P_ROW = 100
MAX_P_FULL = 10

#runnning over all discriminants up to 1000
discs = range(1,MAX_DISC)

# full hecke matrix T_p for p up to 10
ps = {d : [p for p in prime_range(MAX_P_FULL) if (d % p != 0)] for d in discs}
cmds = [["./src/omf5 -genus=data/qf5.db -disc=%d -hecke -p=%d > data/hecke_mat_%d_%d.dat 2> logs/hecke_mat_%d_%d.log &" %(d, p, d, p, d, p) for p in ps[d] ] for d in discs]
cmds = reduce(lambda x,y : x+y, cmds)
for cmd in cmds:
    print "Executing %s" %(cmd)
    os.system(cmd)

# first row of hecke matrix T_p for p up to 100
ps = {d : [p for p in prime_range(MAX_P_ROW) if (d % p != 0)] for d in discs}
cmds = [["./src/omf5 -genus=data/qf5.db -disc=%d -hecke -p=%d -row > data/hecke_row_%d_%d.dat 2> logs/hecke_row_%d_%d.log &" %(d, p, d, p, d, p) for p in ps[d] ] for d in discs]
cmds = reduce(lambda x,y : x+y, cmds)
for cmd in cmds:
    print "Executing %s" %(cmd)
    os.system(cmd)
    
exit()

    
