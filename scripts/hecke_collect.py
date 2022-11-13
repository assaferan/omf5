#! /usr/bin/python

fnames = ["data/hecke_"+str(i) + ".dat" for i in range(1,1000)]
hecke_strs = [open(f).read() for f in fnames]
hecke_str = "[{}," + ",".join(hecke_strs) + "]"
# not really needed, but we want to check validity
hecke = eval(hecke_str)
out_fname = "data/hecke_all.dat"
f = open(out_fname, "w")
f.write(hecke_str)
f.close()
