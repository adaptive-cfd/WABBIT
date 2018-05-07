#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:41:48 2017

@author: engels
"""
def read_wabbit_hdf5(file):
    """ Read a wabbit-type HDF5 of block-structured data.

    Return time, x0, dx, box, data, treecode.
    Get number of blocks and blocksize as

    N, Bs = data.shape[0], data.shape[1]
    """
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

    jmin, jmax = get_max_min_level( treecode )
    N = data.shape[0]
    Bs = data.shape[1]

    print("~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Reading file %s" % (file) )
    print("Time=%e it=%i N=%i Bs=%i Jmin=%i Jmax=%i" % (time, iteration, N, Bs, jmin, jmax) )
    print("~~~~~~~~~~~~~~~~~~~~~~~~~")

    return time, x0, dx, box, data, treecode


def wabbit_error(dir, show=False, norm=2):
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import glob

    files = glob.glob(dir+'/phi_*.h5')
    files.sort()

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5(files[-1])
    # get number of blocks and blocksize
    N, Bs = data.shape[0], data.shape[1]


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

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5(files[-1])

    # get number of blocks and blocksize
    N, Bs = data.shape[0], data.shape[1]

    if os.path.isfile(dir+'timesteps_info.t'):
        d = np.loadtxt(dir+'timesteps_info.t')
        if d.shape[1]==5:
            # old data has 5 cols
            N = np.max(d[:,2])
        else:
            # new data we added the col for cpu time (...seufz)
            N = np.max(d[:,3])
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

def fetch_jmax_dir(dir):
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


    eps=config.get('Blocks','max_treelevel',fallback='0')
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



def plot_wabbit_file( file, savepng=False, savepdf=False, cmap='rainbow', caxis=None, title=True, mark_blocks=True,
                     gridonly=False, contour=False, ax=None, fig=None, ticks=True, colorbar=True, dpi=300 ):
    """ Read a (2D) wabbit file and plot it as a pseudocolor plot.

    Keyword arguments:
        * savepng directly store a image file
        * cmap colormap for data
        * caxis manually specify glbal min / max for color plot
    """
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt


    import h5py

    # read procs table, if we want to draw the grid only
    if gridonly:
        fid = h5py.File(file,'r')
        b = fid['procs'][:]
        procs = np.array(b, dtype=float)
        fid.close()

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5( file )
    # get number of blocks and blocksize
    N, Bs = data.shape[0], data.shape[1]

    # we need these lists to modify the colorscale, as each block usually gets its own
    # and we would rather like to have a global one.
    h = []
    c1 = []
    c2 = []

    if fig is None:
        fig = plt.gcf()
        fig.clf()
    if ax is None:
        ax = fig.gca()

    # clear axes
    ax.cla()

    # if only the grid is plotted, we use grayscale for the blocks, and for
    # proper scaling we need to know the max/min level in the grid
    jmin, jmax = get_max_min_level( treecode )

    for i in range(N):
        if not gridonly:
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])
            block = data[i,:,:].copy().transpose()
            if not contour:
                hplot = ax.pcolormesh( Y, X, block, cmap=cmap, shading='flat' )
            else:
                hplot = ax.contour( Y, X, block, [0.1, 0.2, 0.5, 0.75] )

            # use rasterization for the patch we just draw
            hplot.set_rasterized(True)

            # unfortunately, each patch of pcolor has its own colorbar, so we have to take care
            # that they all use the same.
            h.append(hplot)
            a = hplot.get_clim()
            c1.append(a[0])
            c2.append(a[1])

        if mark_blocks and not gridonly:
            ax.add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False, edgecolor='k' ))

        if gridonly:
            level = treecode_level( treecode[i,:] )
            color = 0.9 - 0.75*(level-jmin)/(jmax-jmin)
            color = plt.cm.jet( procs[i]/np.max(procs) )
#            ax.add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], facecolor=[color,color,color], edgecolor='k' ))
            ax.add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], facecolor=color, edgecolor='k' ))


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

        if colorbar:
            plt.colorbar(h[0], ax=ax)

    if title:
        plt.title( "t=%f Nb=%i Bs=%i" % (time,N,Bs) )


    if not ticks:
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='left',      # ticks along the bottom edge are off
        top='right',         # ticks along the top edge are off
        labelleft='off') # labels along the bottom edge are off

    plt.axis('tight')
    plt.axes().set_aspect('equal')
    plt.gcf().canvas.draw()

    if not gridonly:
        if savepng:
            plt.savefig( file.replace('h5','png'), dpi=dpi, transparent=True )

        if savepdf:
            plt.savefig( file.replace('h5','pdf') )
    else:
        if savepng:
            plt.savefig( file.replace('.h5','-grid.png'), dpi=dpi, transparent=True )

        if savepdf:
            plt.savefig( file.replace('.h5','-grid.pdf') )

#%%
def to_dense_grid( fname_in, fname_out):
    """ Convert a WABBIT grid to a full dense grid in a single matrix.

    We asssume here that interpolation has already been performed, i.e. all
    blocks are on the same (finest) level.
    """
    import numpy as np
    import insect_tools
    import matplotlib.pyplot as plt

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5( fname_in )

    # convert blocks to complete matrix
    field, box = dense_matrix(  x0, dx, data, treecode )

    plt.figure()
    plt.pcolormesh(field)
    plt.axis('equal')
    plt.colorbar()

    # write data to FLUSI-type hdf file
    insect_tools.write_flusi_HDF5( fname_out, time, box, field)

#%%
def dense_matrix(  x0, dx, data, treecode ):
    import numpy as np
    """ Convert a WABBIT grid to a full dense grid in a single matrix.

    We asssume here that interpolation has already been performed, i.e. all
    blocks are on the same (finest) level.

    returns the full matrix and the domain size
    """
    # number of blocks
    N = data.shape[0]
    # size of each block
    Bs = data.shape[1]

    # check if all blocks are on the same level or not
    jmin, jmax = get_max_min_level( treecode )
    if jmin != jmax:
        print("ERROR! not an equidistant grid yet...")

    # note skipping of redundant points, hence the -1
    ny = int( np.sqrt(N)*(Bs-1) )
    nx = int( np.sqrt(N)*(Bs-1) )
    # all spacings should be the same - it does not matter which one we use.
    ddx = dx[0,0]

    print("Number of blocks %i" % (N))
    print("Dense field resolution %i x %i" % (nx, ny) )
    print("Spacing %e domain %e" % (ddx, ddx*nx))

    # allocate target field
    field = np.zeros([nx,ny])

    for i in range(N):
        # get starting index of block
        ix0 = int(x0[i,0]/dx[i,0])
        iy0 = int(x0[i,1]/dx[i,1])

        # copy block content to data field
        field[ ix0:ix0+Bs-1, iy0:iy0+Bs-1 ] = data[i,0:-1,0:-1]

    # domain size
    box = [dx[0,0]*nx, dx[0,1]*ny]

    return(field, box)

#plot_wabbit_file('B_SWIRL/adaptive1_swirl-tough_0_5.0e-5/phi_000000000000.h5', cmap='gray', contour=True, mark_blocks=True)
#plot_wabbit_file('B_SWIRL/adaptive1_swirl-tough_0_5.0e-5/phi_000010000000.h5', cmap='gray', contour=True, mark_blocks=True)
#plot_wabbit_dir('C_SWIRL_new/dx1_equi_4th-4th-4th_4/', savepng=True)
#plot_wabbit_dir('A_CONV/adapt3_4_3.45510729e-06/', savepng=True)
#err = wabbit_error('../exact/')
#print(err)

#to_dense_grid('../shear-mario/fixed/Uy_000003000000.h5', '../shear-mario/fixed/denseUy_000003000000.h5')
#to_dense_grid('../shear-mario/1e-2/Uy_000004000000.h5', '../shear-mario/1e-2/denseUy_000004000000.h5')
#to_dense_grid('../shear-mario/1e-2/Ux_000004000000.h5', '../shear-mario/1e-2/denseUx_000004000000.h5')
