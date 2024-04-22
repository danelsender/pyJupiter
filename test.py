import pyJupiter as jup

import os

OUT = '/home/ipa/ipaszulagyi/users/delsender/adiab_1e-2/ref_tests'

os.chdir(OUT)

N = 120

jup.output(N).plotmidslice(field='gasdensity',polar=True,grid=False)