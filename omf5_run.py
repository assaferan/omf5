#! /usr/bin/python

import os
import sys
import commands
import subprocess

assert(len(sys.argv) in [2,4])

quad = commands.getoutput("sage generate_form.sage " + str(sys.argv[1]))
quad = quad.replace(" ", ",")[:-2]
params = " -quad=" + quad + " -format=GG"
if (len(sys.argv) == 4):
    params += " -form_idx="+sys.argv[2]
    params += " -prec="+sys.argv[3]
os.system("./bin/omf5" + params)

