#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:41:48 2017

@author: engels
"""
def read_wabbit_hdf5(file, return_treecode=False):

    import h5py
    import numpy as np

    fid = h5py.File(file,'r')
    b = fid['coords_origin'][:]
    x0 = np.array(b, dtype=float)

    b = fid['coords_spacing'][:]
    dx = np.array(b, dtype=float)

    b = fid['blocks'][:]
    data = np.array(b, dtype=float)

    b = fid['block_treecode'][:]
    treecode = np.array(b, dtype=float)

    # get the dataset handle
    dset_id = fid.get('blocks')

    # from the dset handle, read the attributes
    time = dset_id.attrs.get('time')
    iteration = dset_id.attrs.get('iteration')
    box = dset_id.attrs.get('domain-size')


    fid.close()

    N = data.shape[0]
    Bs = data.shape[1]

    if return_treecode:
        return data, x0, dx, time, N, Bs, treecode
    else:
        return data, x0, dx, time, N, Bs


def wabbit_error(dir, show=False, norm=2):
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import glob

    files = glob.glob(dir+'/phi_*.h5')
    files.sort()

    data, x0, dx, tmax, N, Bs = read_wabbit_hdf5(files[-1])


    if show:
        plt.close('all')
        plt.figure()
        for i in range(N):
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])
            block = data[i,:,:].copy().transpose()
            y = plt.pcolormesh( Y, X, block, cmap='rainbow' )
            y.set_clim(0,1.0)
            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False ))
        plt.colorbar()
        plt.axis('equal')

    # compute error:
    err = []
    exc = []
    for i in range(N):
        [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0] - 0.75, np.arange(Bs)*dx[i,1]+x0[i,1] -0.5 )
        [X1,Y1] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1] )

        # periodization
        X[X<-0.5] = X[X<-0.5] + 1.0
        Y[Y<-0.5] = Y[Y<-0.5] + 1.0

        exact = np.exp( -((X)**2 + (Y)**2) / 0.01  )
        block = data[i,:,:].copy().transpose()

        err = np.hstack( (err,np.ndarray.flatten(exact - block)) )
        exc = np.hstack( (exc,np.ndarray.flatten(exact)) )

    if show:
        plt.figure()
        c1 = []
        c2 = []
        h = []

        for i in range(N):
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0] - 0.75, np.arange(Bs)*dx[i,1]+x0[i,1] -0.5 )
            [X1,Y1] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1] )

            # periodization
            X[X<-0.5] = X[X<-0.5] + 1.0
            Y[Y<-0.5] = Y[Y<-0.5] + 1.0

            exact = np.exp( -((X)**2 + (Y)**2) / 0.01  )
            block = data[i,:,:].copy().transpose()

            err = np.hstack( (err,np.ndarray.flatten(exact - block)) )
            exc = np.hstack( (exc,np.ndarray.flatten(exact)) )

            y = plt.pcolormesh( Y1, X1, exact-block, cmap='rainbow' )
            a = y.get_clim()

            h.append(y)
            c1.append(a[0])
            c2.append(a[1])

            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False ))

        for i in range(N):
            h[i].set_clim( (min(c1),max(c2))  )

        plt.colorbar()
        plt.axis('equal')


    return( np.linalg.norm(err, ord=norm) / np.linalg.norm(exc, ord=norm) )


def key_parameters(dir):
    import configparser
    import glob

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    print(inifile[0])
    config = configparser.ConfigParser()
    config.read(inifile[0])


    adapt_mesh=config.get('Blocks','adapt_mesh',fallback='0')
    adapt_inicond=config.get('Blocks','adapt_inicond',fallback='0')
    eps= config.get('Blocks','eps',fallback='0')

    namestring = adapt_mesh+adapt_inicond+eps

    print(namestring)


def fetch_dt_dir(dir):
    import glob

    if dir[-1] == '/':
        dir = dir
    else:
        dir=dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if len(inifile) > 1:
        print('ERROR MORE THAN ONE INI FILE')

    ## Open the file with read only permit
    f = open(inifile[0], "r")
    ## use readlines to read all lines in the file
    ## The variable "lines" is a list containing all lines
    lines = f.readlines()
    ## close the file after reading the lines.
    f.close()

    dt = 0.0
    for line in lines:
        if line.find('dt_fixed=') != -1:
            line = line.replace('dt_fixed=','')
            line = line.replace(';','')
            dt = float(line)
    return(dt)



def fetch_Nblocks_dir(dir, return_Bs=False):
    import glob
    import numpy as np
    import os.path

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    files = glob.glob(dir+'/phi_*.h5')
    files.sort()
    data, x0, dx, tmax, N, Bs = read_wabbit_hdf5(files[-1])

    if os.path.isfile(dir+'timesteps_info.t'):
        d = np.loadtxt(dir+'timesteps_info.t')
        N = np.max(d[:,2])
    else:
        raise ValueError('timesteps_info.t not found in dir.'+dir)

    if (return_Bs):
        return(N,Bs)
    else:
        return(N)


def fetch_eps_dir(dir):
    import glob
    import configparser

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    print(inifile[0])
    config = configparser.ConfigParser()
    config.read(inifile[0])


    eps=config.get('Blocks','eps',fallback='0')
    eps = eps.replace(';','')

    return( float(eps) )



def convergence_order(N, err):
    import numpy as np

    if len(N) != len(err):
        raise ValueError('Convergence order args do not have same length')

    A = np.ones([len(err), 2])
    B = np.ones([len(err), 1])
    # ERR = A*N + B
    for i in range( len(N) ) :
        A[i,0] = np.log(N[i])
        B[i] = np.log(err[i])

    x, residuals, rank, singval  = np.linalg.lstsq(A,B)

    return x[0]



def plot_wabbit_dir(d, savepng=False):
    import glob

    files = glob.glob(d+'/phi_*.h5')
    files.sort()

    for file in files:
        plot_wabbit_file(file, savepng)


# given a treecode tc, return its level
def treecode_level( tc ):
    for j in range(len(tc)):
            if (tc[j]==-1):
                break

    level = j - 1 + 1 # note one-based level indexing.
    return(level)


# for a treecode list, return max and min level found
def get_max_min_level( treecode ):
    import numpy as np

    min_level = 99
    max_level = -99
    N = treecode.shape[0]
    for i in range(N):
        tc = treecode[i,:].copy()
        level = treecode_level(tc)
        min_level = min([min_level,level])
        max_level = max([max_level,level])

    return min_level, max_level



def plot_wabbit_file( file, savepng=False, cmap='rainbow', caxis=None, title=True, mark_blocks=True, gridonly=False, contour=False ):
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt



    data, x0, dx, time, N, Bs, treecode = read_wabbit_hdf5(file, return_treecode=True)

    min_level, max_level = get_max_min_level( treecode )
#    print(min_level, max_level)
    h =[]
    c1 = []
    c2 = []
    plt.gcf().clf()

    for i in range(N):
        if not gridonly:
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])
            block = data[i,:,:].copy().transpose()
            if not contour:
                hplot = plt.pcolormesh( Y, X, block, cmap=cmap, shading='flat' )
            else:
                hplot = plt.contour( Y, X, block, [0.1, 0.2, 0.5, 0.75] )

            # unfortunately, each patch of pcolor has its own colorbar, so we have to take care
            # that they all use the same.
            h.append(hplot)
            a = hplot.get_clim()
            c1.append(a[0])
            c2.append(a[1])

        if mark_blocks and not gridonly:
            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False, edgecolor='k' ))

        if gridonly:
            level = treecode_level( treecode[i,:] )
            color = 0.9 - 0.75*(level-min_level)/(max_level-min_level)
            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], facecolor=[color,color,color], edgecolor='k' ))


    if not gridonly:
        # unfortunately, each patch of pcolor has its own colorbar, so we have to take care
        # that they all use the same.
        if caxis is None:
            # automatic colorbar, using min and max throughout all patches
            for hplots in h:
                hplots.set_clim( (min(c1),max(c2))  )
        else:
            # set fixed (user defined) colorbar for all patches
            for hplots in h:
                hplots.set_clim( (min(caxis),max(caxis))  )
        plt.colorbar()

    if title:
        plt.title( "t=%f Nb=%i Bs=%i" % (time,N,Bs) )
#    plt.axis('equal')
    plt.axis('tight')
    plt.gcf().canvas.draw()


    if savepng:
        plt.savefig( file.replace('h5','png') )


def to_dense_grid( fname, ofile):
    import numpy as np
    import insect_tools
    import matplotlib.pyplot as plt
    data, x0, dx, time, N, Bs, treecode = read_wabbit_hdf5(fname, return_treecode=True)

    jmin, jmax = get_max_min_level( treecode )

    if jmin != jmax:
        print("ERROR! not an equidistant grid yet...")

    # note skipping of redundant points, hence the -1
    ny = int( np.sqrt(N)*(Bs-1))
    nx = int( np.sqrt(N)*(Bs-1))
    # all spacings should be the same - it does not matter which one we use.
    ddx = dx[0,0]


    print(np.max(data))

    print("Number of blocks %i" % (N))
    print("Dense field resolution %i x %i" % (nx, ny) )
    print("Spacing %e domain %e" % (ddx, ddx*nx))

    field = np.zeros([nx,ny])

    for i in range(N):
        ix0=int(x0[i,0]/dx[i,0])
        iy0=int(x0[i,1]/dx[i,1])
        field[ ix0:ix0+Bs-1, iy0:iy0+Bs-1 ] = data[i,0:-1,0:-1]

    plt.figure()
    plt.pcolormesh(field)
    plt.axis('equal')
    plt.colorbar()

    insect_tools.write_flusi_HDF5( ofile, time, [ddx*nx, ddx*nx], field)





#plot_wabbit_file('B_SWIRL/adaptive1_swirl-tough_0_5.0e-5/phi_000000000000.h5', cmap='gray', contour=True, mark_blocks=True)
#plot_wabbit_file('B_SWIRL/adaptive1_swirl-tough_0_5.0e-5/phi_000010000000.h5', cmap='gray', contour=True, mark_blocks=True)
#plot_wabbit_dir('C_SWIRL_new/dx1_equi_4th-4th-4th_4/', savepng=True)
#plot_wabbit_dir('A_CONV/adapt3_4_3.45510729e-06/', savepng=True)
#err = wabbit_error('../exact/')
#print(err)

#to_dense_grid('../shear-mario/fixed/Uy_000003000000.h5', '../shear-mario/fixed/denseUy_000003000000.h5')
to_dense_grid('../shear-mario/1e-2/Uy_000004000000.h5', '../shear-mario/1e-2/denseUy_000004000000.h5')
to_dense_grid('../shear-mario/1e-2/Ux_000004000000.h5', '../shear-mario/1e-2/denseUx_000004000000.h5')
