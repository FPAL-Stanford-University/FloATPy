from mpi4py import MPI
import numpy as np
import os
import sys

import floatpy.derivatives.compact.compact_derivative as cd
import floatpy.readers.padeops_reader as por
import floatpy.readers.parallel_reader as pdr
import floatpy.utilities.reduction as red


def reynolds_average(avg, f):
    return avg.average(f)

def favre_average(avg, rho, f, rho_bar=None):
    if rho_bar is None:
        rho_bar = avg.average(rho)

    return avg.average(rho*f) / rho_bar

def integrate_x(f, dx, grid_partition):
    comm = MPI.Comm.f2py(grid_partition.commx())

    myfint = np.array([ np.sum( f ) * dx ])
    fint = np.array([ myfint[0] ])
    comm.Allreduce(myfint, fint, op=MPI.SUM)

    return fint[0]

def integrate_xy(f, dx, dy, grid_partition):
    comm = MPI.Comm.f2py(grid_partition.commxy())

    myfint = np.array([ np.sum( f ) * dx * dy ])
    fint = np.array([ myfint[0] ])
    comm.Allreduce(myfint, fint, op=MPI.SUM)

    return fint[0]

def integrate_volume(f, dx, dy, dz, grid_partition):
    comm = MPI.Comm.f2py(grid_partition.comm3d())

    myfint = np.array([ np.sum( f ) * dx * dy * dz ])
    fint = np.array([ myfint[0] ])
    comm.Allreduce(myfint, fint, op=MPI.SUM)

    return fint[0]

def get_mixing_width_mixedness(Y_N2, Y_CO2, grid_partition, avg, dx, dy):
    # Y_N2, Y_CO2 = reader.readData( ('Massfraction_01', 'Massfraction_02') )

    # Get y-z average
    Y_CO2_bar = avg.average_z(Y_CO2)
    Y_CO2_bar = avg.average_y(Y_CO2_bar)

    mixwidth = integrate_x(4.*Y_CO2_bar*(1.-Y_CO2_bar), dx, grid_partition)

    numerator = avg.average(Y_N2 * Y_CO2)
    denominator = avg.average(Y_N2) * avg.average(Y_CO2)

    mixedness = integrate_xy(numerator, dx, dy, grid_partition) / integrate_xy(denominator, dx, dy, grid_partition)

    return mixwidth, mixedness

def TKE(rho, u, v, w, grid_partition, avg, dx, dy, dz, volume_integrated=True):
    # rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )

    # Get density average
    rho_bar = reynolds_average(avg, rho)

    # Get favre averaged velocities
    u_tilde = favre_average(avg, rho, u, rho_bar=rho_bar)
    v_tilde = favre_average(avg, rho, v, rho_bar=rho_bar)
    w_tilde = favre_average(avg, rho, w, rho_bar=rho_bar)
    
    # Get favre fluctuating velocities (array broadcasting)
    u_pprime = u - u_tilde
    v_pprime = v - v_tilde
    w_pprime = w - w_tilde

    # Get turbulent kinetic energy
    tke = 0.5 * rho * (u_pprime*u_pprime + v_pprime*v_pprime + w_pprime*w_pprime)

    if volume_integrated:
        tke_integrated = integrate_volume(tke, dx, dy, dz, grid_partition)
        return tke_integrated

    tke = reynolds_average(avg, tke)
    return tke

def circulation(vorticity_z, grid_partition, avg, dx, dy, dz):
    z_vorticity = reynolds_average(avg, vorticity_z)

    z_vorticity_pos = np.copy(z_vorticity)
    z_vorticity_pos[z_vorticity < 0.] = 0.

    z_vorticity_neg = np.copy(z_vorticity)
    z_vorticity_neg[z_vorticity > 0.] = 0.

    net_circulation = integrate_xy(z_vorticity,     dx, dy, grid_partition)
    pos_circulation = integrate_xy(z_vorticity_pos, dx, dy, grid_partition)
    neg_circulation = integrate_xy(z_vorticity_neg, dx, dy, grid_partition)

    return net_circulation, pos_circulation, neg_circulation

def get_enstrophy(vorticity, grid_partition, dx, dy, dz):
    enstrophy = vorticity[:,:,:,0]**2 + vorticity[:,:,:,1]**2 + vorticity[:,:,:,2]**2
    enstrophy_integrated = integrate_volume(enstrophy, dx, dy, dz, grid_partition)
    return enstrophy_integrated

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "Usage: "
        print "    python %s <filename_prefix> <outputfile> [start time index (optional)]" %sys.argv[0]
        sys.exit()

    filename_prefix = sys.argv[1]
    outputfile  = sys.argv[2]
    start_index = 0
    if len(sys.argv) >= 4:
        start_index = int(sys.argv[3])

    periodic_dimensions = (True,True,True)
    x_bc = (0,0)
    y_bc = (0,0)
    z_bc = (0,0)

    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    procs = comm.Get_size()

    # Set up the serial Miranda reader
    serial_reader = por.PadeopsReader(filename_prefix, periodic_dimensions=periodic_dimensions)

    # Set up the parallel reader
    reader = pdr.ParallelDataReader(comm, serial_reader)

    # Set up the reduction object
    avg = red.Reduction(reader.grid_partition, periodic_dimensions)

    # Available steps
    steps = sorted(reader.steps)

    x, y, z = reader.readCoordinates()
    dx = x[1,0,0] - x[0,0,0]
    dy = y[0,1,0] - y[0,0,0]
    dz = z[0,0,1] - z[0,0,0]

    # Set up the derivative object
    der = cd.CompactDerivative(reader.grid_partition, (dx, dy, dz), (10, 10, 10), periodic_dimensions)

    if rank == 0:
        line = ('#%26s' + ' %26s'*7) %('Time', 'Mixing width', 'Mixedness', 'TKE', \
                'Net circulation', 'Pos circulation', 'Neg circulation', 'Enstrophy')
        print line
        if start_index == 0:
            fout = open(outputfile, 'w+')
            fout.write( line + '\n' )
        else:
            fout = open(outputfile, 'a')

    for step in steps[start_index:]:
        reader.step = step

        rho, u, v, w = reader.readData( ('rho', 'u', 'v', 'w') )
        vorticity = der.curl(u, v, w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        tke = TKE(rho, u, v, w, reader.grid_partition, avg, dx, dy, dz)
        vorticity = der.curl(u, v, w, x_bc=x_bc, y_bc=y_bc, z_bc=z_bc)
        net_circulation, pos_circulation, neg_circulation = circulation(vorticity[:,:,:,2], reader.grid_partition, avg, dx, dy, dz)
        enstrophy = get_enstrophy(vorticity, reader.grid_partition, dx, dy, dz)

        if rank == 0:
            line = (' %26.14e'*6) %(reader.time, tke, net_circulation, pos_circulation, neg_circulation, enstrophy)
            print line
            fout.write( line + '\n' )
            fout.flush()

    if rank == 0:
        fout.close()
