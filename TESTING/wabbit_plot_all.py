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


wabbit_tools.plot_wabbit_file( 'A_CONV/dx1_nonequi_4th-4th-4th_0/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/dx2_nonequi_2nd-2nd-4th_0/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )

wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_0_1.00000000e-07/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_0_1.00000000e-07/phi_000000500000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_0_1.00000000e-07/phi_000001000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )

wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_17_1.51177507e-05/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_17_1.51177507e-05/phi_000000500000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt1_17_1.51177507e-05/phi_000001000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )

wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_0_2.06913808e-05/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_0_2.06913808e-05/phi_000000500000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_0_2.06913808e-05/phi_000001000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )


wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_11_1.12883789e-03/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_11_1.12883789e-03/phi_000000500000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'A_CONV/adapt2_11_1.12883789e-03/phi_000001000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )

wabbit_tools.plot_wabbit_file( 'B_SWIRL/dt1_equi_4th-4th-4th_1.0e-3/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )
wabbit_tools.plot_wabbit_file( 'B_SWIRL/dt5_nonequi_4th-4th-4th_2.5e-3/phi_000000000000.h5', savepng=True, cmap='rainbow', caxis=[0,1] )

#
#
files = glob.glob('B_SWIRL/dt5_nonequi_4th-4th-4th_1.0e-3/phi_*.h5', recursive=True)
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
#        wabbit_tools.plot_wabbit_file( file, savepng=True, gridonly=True )
    else:
        print('already taken care of')
    print("[%i %i]" % (i, len(files)) )


