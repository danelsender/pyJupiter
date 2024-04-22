import pyJupiter as jup

import os 


OUT = '/home/ipa/ipaszulagyi/users/delsender/adiab_1e-2/adiab_3D_g_1e-2/'
os.chdir(OUT)

N = 110
i = N

min = 1e+10
max = 0

for i in range(N):
    temppmin, tempmax = jup.output(i+1).get_minmax(field='gasdensity')

    if temppmin < min:
        min = temppmin
    if tempmax > max:
        max = tempmax

print(f"The retreived min and max are:\n i = {i} \n Minimum = {min}\n Maximum = {max}")

globmin, globmax = jup.get(N).global_minmax(field='gasdensity')

print(f"The retreived min and max are:\n i = {i} \n Minimum = {globmin}\n Maximum = {globmax}")