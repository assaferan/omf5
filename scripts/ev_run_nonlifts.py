#! /usr/local/bin/sage -python

import os

os.system("export LD_LIBRARY_PATH=/home/assaferan/lib")

MIN_DISC = 1
MAX_DISC = 200
MAX_PREC = 200
k = 3
j = 0

#runnning over all discriminants up to 1000
discs = range(MIN_DISC,MAX_DISC)

nonlift_f = open("scripts/nonlift_idxs.dat")
r = nonlift_f.read()
nonlift_f.close()
nonlift_idxs = eval(r)

base_fname = "hecke_ev_%d_%d_nl_%d_" %(k, j, MAX_PREC)
cmds = ["./src/omf5 -genus=data/qf5db.sage -isom -disc=%d -prec=%d -nonlifts -idxs=%s > data/%s%d.dat 2> logs/%s%d.log &" %(d,MAX_PREC,str(nonlift_idxs[d])[1:-1].replace(" ",""),base_fname,d,base_fname,d) for d in discs if len(nonlift_idxs[d]) > 0]

for cmd in cmds:
    print "Executing %s" %(cmd)
    os.system(cmd)
    
exit()

    
