#! /usr/bin/python

import os
import sys
import commands
import subprocess

assert(len(sys.argv) in [3,5])

quad = commands.getoutput("sage generate_form.sage " + str(sys.argv[1]) + " " +str(sys.argv[2]))
quad = quad.replace(" ", ",")[:-1]
params = " -quad=" + quad + " -format=GG"
if (len(sys.argv) == 5):
    params += " -form_idx="+sys.argv[3]
    params += " -prec="+sys.argv[4]
os.system("./bin/omf5" + params)


