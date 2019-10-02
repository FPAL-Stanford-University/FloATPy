#!/gpfs/mira-home/kmatsuno/floatpy_env/bin/python
import numpy as np
import os
import sys
import get_namelist as nml
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: "
        print "  python {} <prefix> [tID_start(default=0)] ".format(sys.argv[0])
        sys.exit()
    filename_prefix = sys.argv[1]
    start_index = 0
    if len(sys.argv) > 2:
        start_index = int(sys.argv[2])
    
    # setup the inputs object
    dirname = os.path.dirname(filename_prefix)
    verbose=True
    inp = nml.inputs(dirname,verbose)
    du = inp.du
    
    # Get Re_theta
    dat = np.fromfile(filename_prefix+'growth.dat',sep=' ')
    n = np.size(dat)
    dat = np.reshape(dat,[n/3,3])
    time = dat[:,0]
    dtheta = dat[:,1]
    Re_theta = inp.r_ref*inp.du*dtheta/inp.mu_ref

    # Get Re_omega
    dat = np.fromfile(filename_prefix+'domega.dat',sep=' ')
    n = np.size(dat)
    dat = np.reshape(dat,[n/2,2])
    time = dat[:,0]
    domega = dat[:,1]
    Re_omega = inp.r_ref*inp.du*domega/inp.mu_ref

    print('t\tdomega\t\tRe_omega')
    for i in range(len(time)):
        print('%i\t%2d\t%2d\t%2d\t%2d'%(time[i],dtheta[i],domega[i],Re_omega[i],Re_theta[i]))
