#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 18:17:56 2017

@author: engels
"""

import wabbit_tools
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import glob

def do_test(rootdir, name):
    dirsx = glob.glob(rootdir+'*')

    EPS=[]
    err=[]
    Nblocks=[]

    for d in dirsx:
        EPS.append( wabbit_tools.fetch_eps_dir(d) )
        err.append( wabbit_tools.wabbit_error(d) )
        Nblocks.append( wabbit_tools.fetch_Nblocks_dir(d) )

    # sort the lists (by eps)
    EPS, err, Nblocks = zip(*sorted(zip(EPS, err, Nblocks)))

    name = name +" [%2.2f]" % (wabbit_tools.convergence_order(EPS,err))

    plt.figure(1)
    plt.loglog( EPS, err, label=name, marker='o')

    plt.figure(2)
    plt.loglog( EPS, Nblocks, label=name, marker='o')


#%% convection
plt.close('all')
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.serif"] = 'Times'

do_test("A_CONV/adapt1","adaptive 4th-4th-4th")
do_test("A_CONV/adapt2","adaptive 2nd-2nd-4th")


plt.figure(1)
plt.title('Convection test')
plt.grid()
plt.xlabel('$\\varepsilon$')
plt.ylabel('$||\phi-\phi_{\\mathrm{ex}}||_2/||\phi_{\\mathrm{ex}}||_2$')
plt.legend()#loc='upper left', bbox_to_anchor=(0,1.02,1,0.2), prop={'size': 6})
plt.gcf().subplots_adjust(top=0.82)
plt.gcf().savefig('conv-adaptive-eps.pdf')

plt.figure(2)
plt.title('Convection test')
plt.grid()
plt.xlabel('$\\varepsilon$')
plt.ylabel('$\max\\left({N_{b}(t)}\\right)$')
plt.legend()#loc='upper left', bbox_to_anchor=(0,1.02,1,0.2), prop={'size': 6})
plt.gcf().subplots_adjust(top=0.82)
plt.gcf().savefig('conv-adaptive-nblocks.pdf')


#%% swirl
plt.close('all')


do_test("B_SWIRL/adaptive1_swirl-tough","SWIRL-tough-4th-4th-4th (T=10.0, CFL=1.0)")


do_test("B_SWIRL/adapt_jmax4_swirl-easy","4th-4th-4th (T=4.0, Jmax=4)")
do_test("B_SWIRL/adapt_jmax5_swirl-easy","4th-4th-4th (T=4.0, Jmax=5)")
do_test("B_SWIRL/adapt_jmax6_swirl-easy","4th-4th-4th (T=4.0, Jmax=6)")
do_test("B_SWIRL/adapt_jmax7_swirl-easy","4th-4th-4th (T=4.0, Jmax=7)")
do_test("B_SWIRL/adaptive2_swirl-easy","4th-4th-4th (T=4.0, Jmax=$\infty$)")

plt.figure(1)
plt.title('Swirl test')
plt.grid()
plt.xlabel('$\\varepsilon$')
plt.ylabel('$||\phi-\phi_{\\mathrm{ex}}||_2/||\phi_{\\mathrm{ex}}||_2$')
plt.legend()#loc='upper left', bbox_to_anchor=(0,1.02,1,0.2), prop={'size': 6})
plt.gcf().subplots_adjust(top=0.82)
plt.gcf().savefig('swirl-adaptive-eps.pdf')

plt.figure(2)
plt.title('Swirl test')
plt.grid()
plt.xlabel('$\\varepsilon$')
plt.ylabel('$\max\\left({N_{b}(t)}\\right)$')
plt.legend()#loc='upper left', bbox_to_anchor=(0,1.02,1,0.2), prop={'size': 6})
plt.gcf().subplots_adjust(top=0.82)
plt.gcf().savefig('swirl-adaptive-nblocks.pdf')