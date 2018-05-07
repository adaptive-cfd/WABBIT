#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 18:17:56 2017

@author: engels
"""
import matplotlib.pyplot as plt
import glob
import wabbit_tools

def do_test( basedir, name,  blob_width=0.01):
    dirs = glob.glob( basedir+'*')

    if dirs != []:
        dt=[]
        err=[]
        for dir in dirs:
            if len(glob.glob(dir+'/*.h5')) != 0:
                dt.append( wabbit_tools.fetch_dt_dir(dir) )
                err.append( wabbit_tools.wabbit_error(dir) )
            else:
                print('There are no h5 files in '+dir)

        print(dirs)

        dt, err = zip(*sorted(zip(dt, err)))

        dt = list(dt)
        err = list(err)
        print(dt)
        plt.loglog(dt, err, label=name+" [%2.2f] " % (wabbit_tools.convergence_order(dt,err)), marker='o')

plt.close('all')
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.serif"] = 'Times'

do_test("B_SWIRL/dt1_equi_4th-4th-4th_","equi-4th-4th-4th")
do_test("B_SWIRL/dt2_equi_2nd-2nd-4th_","equi-2nd-2nd-4th")
do_test("B_SWIRL/dt3_equi_4th-4th-3rd_","equi-4th-4th-3rd")
do_test("B_SWIRL/dt5_nonequi_4th-4th-4th","nonequi-4th-4th-4th")

plt.grid()
plt.title('Swirl test')
plt.xlim([1e-5,1e-2])
plt.legend()
plt.xlabel('$\Delta t$')
plt.ylabel('$||\phi-\phi_{\\mathrm{ex}}||_2/||\phi_{\\mathrm{ex}}||_2$')
plt.savefig('time_convergence.pdf', format='pdf')
plt.savefig('time_convergence.png', format='png')
