#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 23:01:03 2021

@author: engels
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#---------------------------------------------------------
# gridfile = 'grid_3levelsbs5.h5'
gridfile = 'grid_2levelsBs13.h5'
# gridfile = 'grid_3levelsbs9.h5'
# gridfile = 'grid_equiBs13.h5'
#---------------------------------------------------------

with open( gridfile+'.'+'info.txt', 'r') as f:
    data = f.readlines()
header = data[0]
print(header)

fig = plt.figure(constrained_layout=False)
fig.set_size_inches( [15, 25] )
gs = GridSpec(3, 2, figure=fig, bottom=0.02, hspace=0.276, top=0.95)
ax1 = fig.add_subplot(gs[0, 0])    
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1:3, :])
fig.suptitle(header)


# plot the grid (FYI)
d = np.loadtxt(gridfile+'.'+'operator_grid_points.txt')

ax1.scatter( d[:,2], d[:,3], s=5 )
ax1.axis('equal')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('grid')



#%%
d = np.loadtxt(gridfile+'.'+'operator_matrix.txt')
N = int(np.max(d[:,0:1+1]))
D = np.zeros([N,N])

for i in range(d.shape[0]):
    # python zero based indexing
    ix, iy = int(d[i,0])-1, int(d[i,1])-1
    D[ix,iy] = d[i,2]
   
   
# ax3.spy(D)
# ax3.grid(True)
# ax3.set_title( header+' sync full')

v = np.linalg.eigvals(D)

ax2.plot( np.real(v), np.imag(v), 'o', mfc='none')
ax2.grid(True)
ax2.set_xlabel('$\\lambda_R$')
ax2.set_ylabel('$\\lambda_I$')
ax2.set_title('Operator eigenvalues (first derivative)')
    
# print(np.linalg.matrix_rank(D))
    
# zero lines removal
D2 = D
D2 = D[~np.all(D == 0, axis=1)].copy()
D2 = D2.transpose()
D2 = D2[~np.all(D2 == 0, axis=1)]#.copy()
D2 = D2.transpose()

ax3.spy(D2)
ax3.set_title('Operator (empty rows&cols removed)')

if '2nd' in header:
    order = '2nd'
elif '4th' in header:
    order = '4th'
elif 'CDF44'in header:
    order = 'CDF44'
    
if 'coarseWins=F' in header:
    method = 'fineWins'
else:
    method = 'coarseWins'
    
outfile = gridfile + '.' + method + '_' + order +'.pdf'
plt.savefig(outfile)



# #%%

# du = np.loadtxt(gridfile+'.'+'u_function.txt')
# du_dx = np.loadtxt(gridfile+'.'+'u_derivative.txt')

# n = int(np.max(du[:,0]))

# u = np.zeros([n])
# u_dx = np.zeros([n])

# for i in range(du.shape[0]):
#     j = int(du[i,0]) - 1 
#     u[j] = du[i,1]
    
# for i in range(du_dx.shape[0]):
#     j = int(du_dx[i,0]) - 1 
#     u_dx[j] = du_dx[i,1]


# error = np.max(np.abs( u_dx - np.matmul(D.transpose(),u) ))
# print(error)