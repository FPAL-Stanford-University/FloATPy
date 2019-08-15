#!/bin/python
import sys
import os

# A class for accessing vars given a namelist
class inputs:
    def __init__(self,dirname,verbose=False):
        self.Mc,self.Re,self.rr = read_namelist(dirname,verbose)
        
        gam = 1.4
        p_ref = 1.
        T_ref = 1.
        r_ref = 1.
        mu_ref = 1./self.Re
        Rgas1 = (1.+self.rr)/2.			
        Rgas2 = (1.+self.rr)/(2.*self.rr) 
        c1 = (gam*p_ref/(r_ref/Rgas1))**0.5
        c2 = (gam*p_ref/(r_ref/Rgas2))**0.5
        
        self.du = self.Mc*(c1+c2)
        self.gam = gam  
        self.p_ref = p_ref 
        self.T_ref = T_ref 
        self.r_ref = r_ref 
        self.mu_ref = mu_ref

    def grid_params(self,dirname,verbose=False):
        return read_grid_params(dirname,verbose)
    
    def dimensionalize(this,rho_ref,p_ref,T_ref,mu_ref,Rgas,verbose=True):
        Rgas1 = (1.+self.rr)/2.			
        Rgas2 = (1.+self.rr)/(2.*rr) 
        c1 = (gam*p_ref/(r_ref/Rgas1))**0.5
        c2 = (gam*p_ref/(r_ref/Rgas2))**0.5
        u1 = Mc*c1
        u2 = Mc*c2
        du = Mc*(c1+c2)
        dtheta0 = 1
        Re = dtheta0*du/mu_ref
        


# Assumes the line has a comment, =, and D#
def get_param(line,verbose=False):
    original = line
    line = list(filter(None,line.split('!')))
    line = list(filter(None,line[0].split('=')))
    line = list(filter(None,line[1].split(' ')))
    try:
        line = list(filter(None,line[0].split('D')))
        mult = 10**(int(line[1]))
    except Exception as e:
        if verbose:
            print("Note: Missing D# multiplier for line:")
            print("\t{}".format(original))
            print("\tAssuming multiplier is 1.\n")
        mult = 1.0
    return float(line[0])*mult


# Gets Mc, rr, Re from namelist
def read_namelist(dirname,verbose):
    fname = dirname+"/input.dat"
    with open('{}'.format(fname),'r') as file:
        data_raw = file.read().splitlines()

    # What lines to get
    for line in data_raw:
        if 'Mc ' in line: Mc = get_param(line,verbose)
        if 'Re ' in line: Re = get_param(line,verbose)
        if 'rho_ratio ' in line: rr = get_param(line,verbose)

    if verbose:
        print("Params for this run:")
        print("\tMc = {}".format(Mc))
        print("\tRe = {}".format(Re))
        print("\trr = {}".format(rr))
    return Mc,Re,rr

# Gets Mc, rr, Re from namelist
def read_grid_params(dirname,verbose):
    fname = dirname+"/input.dat"
    with open('{}'.format(fname),'r') as file:
        data_raw = file.read().splitlines()

    # What lines to get
    for line in data_raw:
        if   'nx ' in line: Nx = get_param(line,verbose=False)
        elif 'ny ' in line: Ny = get_param(line,verbose=False)
        elif 'nz ' in line: Nz = get_param(line,verbose=False)
        elif 'Lx ' in line: Lx = get_param(line,verbose=False)
        elif 'Ly ' in line: Ly = get_param(line,verbose=False)
        elif 'Lz ' in line: Lz = get_param(line,verbose=False)

    if verbose:
        print("Params for this grid:")
        print("\tN = {}x{}x{}".format(Nx,Ny,Nz))
        print("\tL = {}x{}x{}".format(Lx,Ly,Lz))
    return Nx,Ny,Nz,Lx,Ly,Lz

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("Usage: {} dirname".format(sys.argv[0]))
        sys.exit()
        
    dirname = sys.argv[1]
    inp = inputs(dirname,verbose=True) 
    read_grid_params(dirname,verbose=True) 




