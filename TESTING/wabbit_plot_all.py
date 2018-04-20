#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 15:43:42 2017

@author: engels
"""

import glob
import wabbit_tools
import os.path
import matplotlib.pyplot as plt

files = glob.glob('./B_SWIRL/adaptive1_swirl-tough_0_5.0e-5/phi_*.h5', recursive=True)
#files = glob.glob('./B_SWIRL/adaptive1*/phi_000005000000.h5', recursive=True)
#files.extend(glob.glob('./B_SWIRL/adaptive1*/phi_000002500000.h5', recursive=True))
#files.extend(glob.glob('./B_SWIRL/adaptive1*/phi_000010000000.h5', recursive=True))
#files = glob.glob('../swirl/meso/swirl_1*/phi_*.h5', recursive=True)

plt.close('all')
for file, i in zip(files, range(len(files))):
    f2 = file.replace('h5','png')

#    wabbit_tools.plot_wabbit_file( file, savepng=True, cmap='rainbow', caxis=[0,1] )
    print(f2)

    if not os.path.isfile( f2 ):
        wabbit_tools.plot_wabbit_file( file, savepng=True, cmap='rainbow', caxis=[0,1] )
        wabbit_tools.plot_wabbit_file( file, savepng=True, gridonly=True )
    else:
        print('already taken care of')
    print("[%i %i]" % (i, len(files)) )


