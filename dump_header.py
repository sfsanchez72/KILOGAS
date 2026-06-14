#!/usr/bin/env python
# coding: utf-8

# In[2]:


from astropy.io import fits
import numpy as np
from typing import Tuple
import math
import sys


# %%
nargs=len(sys.argv)
if (nargs==2):
    file=sys.argv[1]
    ext=0
else:
    if (nargs==3):
        file=sys.argv[1]
        ext=sys.argv[2]
    else:   
        print('USE: new_diga_Pipe3D.py FILE [extension]')
        quit()
# %%
#DIR = 'data/'
hdu=fits.open(f'{file}')

hdu.info()

print('##########################')
print(f'Header of extension {ext}:')
print('##########################')
for key in hdu[ext].header.keys():
    print(f'{key}: {hdu[ext].header[key]}')



# %%
