#!/usr/bin/python
import h5py
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import floatpy.readers.padeops_reader as por

def visualize(reader, step, savefig):
    reader.step = step
    zslice = 0
    fs = 8
    vel = reader.readData(('u', 'v', 'w'))
    colormap = ['inferno','afmhot','viridis']
    
    fig, axarr = plt.subplots(1,3, figsize=(8,25), dpi=200)
    for icol in range(0,3):
        im = axarr[icol].imshow( np.transpose(vel[icol][:,:,0]),
                        cmap=colormap[icol], origin='lower', interpolation='none')#, aspect=1.)\n",
        divider1 = make_axes_locatable(axarr[icol])
        #axarr[icol].set_ylim([])\n",
        cax = divider1.append_axes('right', size='2%', pad=0.1)
        cbar = plt.colorbar(im, cax)
        cbar.ax.tick_params(labelsize=8)
        plt.subplots_adjust(wspace=0.7,hspace=0.1)

    if savefig:
        savefig('slices.png')
    else: 
        plt.show()


# Check input args
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: ')
        print(' python {} <fname>'.format(sys.argv[0]) )
        sys.exit()

    fname = sys.argv[1]
    directory = './'
    if '/' in fname: 
        directory = list(filter(None,fname.split('shearlayer')))[0]

    tID = list(filter(None,fname.split('_')))[-1]
    tID = list(filter(None,tID.split('.h5')))[-1]
    tID = int(tID)

    reader = por.PadeopsReader(directory+'/shearlayer_', periodic_dimensions=(False,False,True))
    print('Grid size: {}'.format(reader.domain_size))
    
    # Only get a single x-y slice
    zslice = 0
    reader.sub_domain = (0,0,zslice), (reader.domain_size[0]-1, reader.domain_size[1]-1, zslice)
    x, y, z = reader.readCoordinates()
    steps = sorted(reader.steps)
   

    visualize(reader, step=tID, savefig=False)
