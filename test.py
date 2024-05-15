import pyJupiter as jup

import os

OUT = '/home/ipa/ipaszulagyi/users/delsender/adiab_1e-2/adiab_3D_g_1e-2'

os.chdir(OUT)

N = 130

#jup.output(N).plotmidslice(field='gasdensity',polar=True,grid=False)
jup.output(N).hill_sphere_mass()
