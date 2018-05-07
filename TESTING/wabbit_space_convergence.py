#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 18:17:56 2017

@author: engels
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import wabbit_tools

def do_test2(rootdir, title):
    dirsx = glob.glob( rootdir+'*' )
    dirsx.sort()
    print(dirsx)

    Nblocks=[]
    e=[]

    if len(dirsx)==0:
        raise ValueError('no data')

    for d in dirsx:
        if len(glob.glob(d+'/*.h5')) != 0:
            err = wabbit_tools.wabbit_error( d )#, norm=np.inf)
            N, Bs = wabbit_tools.fetch_Nblocks_dir(d, return_Bs=True)
            Nblocks.append( np.sqrt(N*Bs) )
            e.append( err )

    Nblocks, e = zip(*sorted(zip(Nblocks, e)))

    plt.loglog( Nblocks, e, label=title+" [%2.2f]" % (wabbit_tools.convergence_order(Nblocks,e)), marker='o')


plt.close('all')
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.serif"] = 'Times'

do_test2('A_CONV/dx1','nonequi (4th-4th-4th)')
do_test2('A_CONV/dx2','nonequi (2nd-2nd-4th)')
do_test2('A_CONV/dx3','equi (4th-4th-4th)')
do_test2('A_CONV/dx4','equi (2nd-2nd-4th)')

plt.xlim([1e0, 1e3])
plt.title('Convection test')
plt.grid()
plt.legend() #loc='upper left', bbox_to_anchor=(0,1.02,1,0.2), prop={'size': 6})
plt.xlabel('$\\sqrt{N_{b} \cdot B_{s}}')
plt.ylabel('$||\phi-\phi_{\\mathrm{ex}}||_2/||\phi_{\\mathrm{ex}}||_2$')
plt.gcf().subplots_adjust(top=0.76)
plt.gcf().savefig('convection-space-nblocks.pdf')
