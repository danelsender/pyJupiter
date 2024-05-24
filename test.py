import numpy as np

M0 = 1.989e+33
R0 = 1.496e+13 * 5.2

SIGMA0 = M0/R0**2

Sig = 200

print(f'Sigma0 has units: {SIGMA0}')
print(f'R0 has units: {R0:.3e}')
print(f'M0 has units: {M0}')
print(f'Surface density of {Sig} g cm$^{-2}$ is {Sig / SIGMA0:.3e} in code units')
