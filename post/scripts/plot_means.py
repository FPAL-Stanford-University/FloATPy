#!/bin/python
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import get_namelist as nml
from common import *

xdir = 0
zdir = 2

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_list] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        tID_list = map(int, sys.argv[2].strip('[]').split(',')) 

    # setup the inputs object, get grid info
    dirname = os.path.dirname(filename_prefix)
    Nx,Ny,Nz,Lx,Ly,Lz = nml.read_grid_params(dirname,verbose=True)
    Ny = int(Ny)
    y = np.linspace(Ly/2.,-Ly/2.,Ny)
    inp = nml.inputs(dirname,verbose=True)
    du = inp.du
   
    fig,ax = plt.subplots(1,2,figsize=(6,3),dpi=100)
    alpha = 0.4 
    # Read the stats at each time step
    for tID,count in zip(tID_list,range(len(tID_list))):
        # Get dtheta
        time = nml.read_time(dirname,tID)
        eta = y/get_dtheta(dirname,time)

        # Read utilde
        utilde = np.fromfile(dirname + '/shearlayer_utilde_%04d.dat'%tID, sep=' ')
        
        # Read cbar
        cbar = np.fromfile(dirname + '/shearlayer_cbar_%04d.dat'%tID, sep=' ')

        # Read ru,rv,re 
        ru = np.fromfile(dirname + '/shearlayer_mean_ru_%04d.dat'%tID, sep=' ')
        rv = np.fromfile(dirname + '/shearlayer_mean_rv_%04d.dat'%tID, sep=' ')
        re = np.fromfile(dirname + '/shearlayer_mean_re_%04d.dat'%tID, sep=' ')

        # Read Rij and plot on ax[1]
        dat = np.fromfile(dirname + '/shearlayer_Rij_%04d.dat'%tID, sep=' ')
        Rij = np.reshape(dat,[Ny,6])

        ax[0].plot(eta,utilde,alpha=alpha,color='C0')
        ax[0].plot(eta,cbar,alpha=alpha,color='C1')
        ax[0].plot(eta,ru,alpha=alpha,color='C2')
        ax[0].plot(eta,rv,alpha=alpha,color='C3')
        ax[0].plot(eta,re,alpha=alpha,color='C4')
        for i in [0,1,3,5]:
            ax[1].plot(Rij[:,i])
   
    for a in ax:
        a.grid(True)
        a.set_xlabel(r'$\eta$')
    ax[0].legend(['utilde','cbar','ru','rv','re'])
    ax[1].legend(['$R_{11}$','$R_{12}$','$R_{22}$','$R_{33}$'])
    plt.show()
    plt.savefig(dirname+'/img/means.png',dpi=200,bbox_inches='tight')


