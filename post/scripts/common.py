#!/bin/python
import numpy as np
import get_namelist as nml
import statistics as stats

myblue=[51./255, 105./255, 169./255]
myred=[200./255, 90./255, 90./255]
fs = 12
xdir = 0
ydir = 1
zdir = 2

# 1D integration in y using midpoint rule:
def integrate_y(y,f):
    ny = np.size(y)
    dy = abs(y[0]-y[1])
    I = 0.
    for i in range(1,ny-1):
        I += dy*0.5*(f[i-1]+f[i+1])
    return I

# Get decorrelation length
def get_lscale(y,Lvv,thresh=0.20):
    Ny = np.size(Lvv)
    imax0 = np.argmax(Lvv)
    # find inflection point
    dLdx = (Lvv[:-1:]-Lvv[1::])/abs(y[0]-y[1])
    asign = np.sign(dLdx)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    imax1,imax2=None,None
                                    
    for i in range(1,Ny/4):
        if signchange[imax0-i]==1:
            imax1 = imax0-i
            break
    for i in range(1,Ny/4):
        if signchange[imax0+i]==1:
            imax2 = imax0+i
            break
    if (imax1 is None): 
        #print('Error finding inflection points 1')
        imax1 = 0
    if (imax2 is None): 
        #print('Error finding inflection points 2')
        imax2 = Ny-1
    i1 = np.argmin(abs(Lvv[imax1:imax0]-thresh)) + imax1
    i2 = np.argmin(abs(Lvv[imax0:imax2]-thresh)) + imax0 
    L = abs(y[i1]-y[i2])    
    return L,i1,i2

# Get 99% thickness. y must be from Ly/2,-Ly/2
def get_L99(y,utilde):
    if (y[0]<y[-1]): print('y must be from Ly/2,-Ly/2')
    utilde = np.squeeze(utilde)
    du = abs(utilde[-1]-utilde[0])
    utop = 0.99*du/2.
    ubot = -0.99*du/2.
    Ny = np.size(y)
    ibot = np.argmin(abs(utilde[:Ny/2]-ubot)[::])
    itop = np.argmin(abs(utilde[Ny/2:]-utop)[::])+Ny/2
    L99 = abs(y[itop]-y[ibot])
    return L99, itop, ibot

# get centerline from utilde=0 
def get_centerline(directory,y,tID):
    if (y[0]<y[-1]): print('y must be from Ly/2,-Ly/2')
    utilde = np.fromfile(directory + '/shearlayer_utilde_%04d.dat'%tID,sep=' ')
    ic = np.argmin(abs(utilde))
    return ic,y[ic]

# Thickness and growth rates
def get_dtheta(directory,time):
    filename_prefix = directory+'/shearlayer_'
    tlist,dtheta,rate = growth_rates(filename_prefix)
    idx = np.argmin(abs(tlist-time))
    return dtheta[idx]

def growth_rates(filename_prefix):
    fname = filename_prefix+'growth.dat'
    dat = np.fromfile(fname,dtype=float, count=-1, sep=' ')
    n = np.size(dat)
    nstats=3
    dat = np.reshape(dat,[n/nstats,nstats])
    time = dat[:,0]
    dtheta = dat[:,1]
    rate = dat[:,2]
    return time,dtheta,rate

def which_tID(tvec,t):
    tID = np.argmin(np.abs(tvec-t))
    return tID

def smooth_modes(f,nmodes):
    fhat = np.fft.fft(f)
    fhat[nmodes+1:-(nmodes+1):] = 0
    f = np.fft.ifft(fhat)
    return f

def merge_dicts(old1,old2):
    new = {}
    for i in range(3):
        key = old1.keys()[i]
        v1 = old1.values()[i][0]
        v2 = old2.values()[i][0]
        e1 = old1.values()[i][1]
        e2 = old2.values()[i][1]
        err = (e1**2+e2**2)**0.5
        new[key] = [(v1+v2)/2., err]
    return new

def grid_res(x,y,z):
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]
    return dx,dy,dz
   
# Transpose from 3d to y
def transpose2y(numSet,q3d): 
    bufy = np.empty(numSet.chunk_y_size, 
            dtype=np.dtype(np.float64, align=True), order='F')
    numSet.grid_partition.transpose_3d_to_y(q3d,bufy)
    q3d = None
    return bufy

def get_qpp(reader,avg,varname):
    r,q = reader.readData(('rho',varname))
    qtilde = stats.favre_average(avg,r,q)
    return q-qtilde[None,:,None]

def get_qp(reader,avg,varname):
    q = reader.readData(varname)
    qbar = stats.reynolds_average(avg,q[0])
    return q-qbar[None,:,None]
