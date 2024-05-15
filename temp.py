import numpy as np


# with open('test.txt','r') as f:
#     lines = f.readlines()

# dim = lines[0]

# xdim = int(dim.split(" ")[0])

# print(xdim)

arr = np.arange(0,60,1)

_zdim = 3
_ydim = 4
_xdim = 5

arr_reshape = np.reshape(arr,(_zdim,_ydim,_xdim))

print(arr)

print(arr_reshape)