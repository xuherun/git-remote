# -*- coding: utf-8 -*-
# @Time    : 2022/8/25 14:44
# @Author  : Herun Xu
# @File description: This is a wind vector average mapping for Vietnam.

import matplotlib.pyplot as plt
import cartopy.mpl.ticker as cticker
import numpy as np
import xarray as xr
with xr.open_dataset('ERA52001-2021-10muv-monthly-means.nc') as data:
    u = data.u10
    v = data.v10

umeans = np.mean(u,axis=0)
vmeans = np.mean(v)
print(umeans)
print(v)


