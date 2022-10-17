#! /usr/bin/python

import os
import sys
import commands

assert(len(sys.argv) == 3)

quad = commands.getoutput("sage generate_form.sage " + str(sys.argv[1]) + " " + str(sys.argv[2]))
quad = quad.replace(" ", ",")[:-1]
os.system("./bin/omf5 -quad=" + quad + " -format=GG")

