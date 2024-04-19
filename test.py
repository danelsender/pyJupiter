import pyJupiter as jup

import os

OUT = '/home/ipa/ipaszulagyi/users/delsender/adiab_1e-3/adiab_3D_g_1e-3'

os.chdir(OUT)

N = 100

jup.output(N).hill_sphere_mass(field='gasdensity', planetmass=0.001)