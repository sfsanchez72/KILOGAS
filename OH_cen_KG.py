#!/usr/bin/env python
# coding: utf-8

# In[2]:


from astropy.io import fits
import numpy as np
from typing import Tuple
import math


# In[39]:


def create_circular_mask(h: int, w: int, center: Tuple[int, int], radius: int) -> np.ndarray:
    """Create a boolean mask for a circle centered at `center` with given `radius`."""
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y - center[1])**2)
    return dist_from_center <= radius

def extract_values_from_cube(fits_path: str, x: int, y: int, radius: int) -> np.ndarray:
    """Extract average values (or quadratic average for errors) from FITS datacube slice-by-slice."""
    
    with fits.open(fits_path) as hdul:
        data = hdul[0].data  # shape: (z, y, x)
        header = hdul[0].header
    desc_entries = extract_desc_entries(header)
    
    if data.ndim != 3:
        raise ValueError("Input FITS file must contain a 3D datacube.")

    nz, ny, nx = data.shape

    if not (0 <= x < nx and 0 <= y < ny):
        x = int(nx/2.0)
        y = int(ny/2.0)
#        raise ValueError("Provided (x, y) coordinates are outside the image bounds.")

#    print(x,y)
    mask = create_circular_mask(ny, nx, center=(x, y), radius=radius)

    result = []

    for z in range(nz):
        desc_key = f'DESC{z}'
        desc = header.get(desc_key, '')

        slice_data = data[z, :, :]
        values = slice_data[mask]

        if values.size == 0:
            result.append(np.nan)
            continue

        if desc.startswith('e_'):
            quadratic_mean = math.sqrt(np.mean(values ** 2))
            result.append(quadratic_mean)
        else:
            simple_mean = np.mean(values)
            result.append(simple_mean)

    return np.array(result),desc_entries

def extract_desc_entries(header):
    """Extract DESC# entries from FITS header as an array sorted by slice index."""
    desc_entries = []
    z_indices = []

    for key in header.keys():
        if key.startswith('DESC') and key[4:].isdigit():
            z_indices.append(int(key[4:]))

    z_indices.sort()

    for z in z_indices:
        desc_key = f'DESC{z}'
        desc_entries.append(header.get(desc_key, ''))

    return desc_entries


# In[64]:


HP_dir = 'data'
x_pixel = -1
y_pixel = -1
radius = 2.5
I = 0
for tab_pe_now in tab_pe:
    kgas_id = tab_pe_now['KGAS_ID']
    cubename = tab_pe_now['cubename']
    fits_file = f'{HP_dir}/{cubename}.OH.cube.fits.gz'
    try:
        extracted_values,keys = extract_values_from_cube(fits_file, x_pixel, y_pixel, radius)    
        if (I==0):
            tab_keys = []
            tab_keys.append('KGAS_ID')
            tab_keys.append('cubename')
            for key in keys:
                tab_keys.append(key)
            tab_keys = np.array(tab_keys)
            tab_OH_cen = Table(names=tab_keys)
            tab_OH_cen['KGAS_ID'] = tab_OH_cen['KGAS_ID'].astype('str')
            tab_OH_cen['cubename'] = tab_OH_cen['cubename'].astype('str')
        tab_values = []
        tab_values.append(kgas_id)
        tab_values.append(cubename)
        for vals in extracted_values:
            tab_values.append(vals)         
        tab_values = np.array(tab_values)
        tab_OH_cen.add_row(tab_values)
        I = I+1
    except:
        print(f'file {fits_file} not found')
    


# In[65]:


tab_OH_cen


# In[68]:


tab_OH_cen.write('tables/tab_OH_cen.ecsv',overwrite=True,delimiter=',')


# In[ ]:




