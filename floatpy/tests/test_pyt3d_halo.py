from mpi4py import MPI
import numpy
import unittest

import floatpy.parallel.pyt3d.pyt3d as pyt3d

#------------------------------------
# Set the MPI communicator
#-------------------------------------
comm = MPI.COMM_WORLD
fcomm = MPI.COMM_WORLD.py2f()

nx = 8
ny = 16
nz = 32

periodic = numpy.array([True, True, True])
nghosts  = numpy.array([1,1,1], dtype=numpy.int32, order='F')
fail = numpy.array([True])
reorder = True

# gp = pyt3d.t3dmod.t3d(fcomm, nx, ny, nz, px, py, pz, periodic, reorder, fail, nghosts=nghosts, createcrosscommunicators=True)
gp = pyt3d.t3dmod.t3d(fcomm, nx, ny, nz, periodic, nghosts=nghosts)

# Get # of procs in x, y and z
px = gp.px()
py = gp.py()
pz = gp.pz()

# Get cartesian communicators
comm_x  = MPI.Comm.f2py( gp.commx () )
comm_y  = MPI.Comm.f2py( gp.commy () )
comm_z  = MPI.Comm.f2py( gp.commz () )
comm_xy = MPI.Comm.f2py( gp.commxy() )
comm_yz = MPI.Comm.f2py( gp.commyz() )
comm_xz = MPI.Comm.f2py( gp.commxz() )

# Get size of interior subdomain of this processor and start and end of the interior subdomain
sz3d = numpy.zeros(3, dtype=numpy.int32, order='F')
st3d = numpy.zeros(3, dtype=numpy.int32, order='F')
en3d = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_sz3d(sz3d)
gp.get_st3d(st3d)
gp.get_en3d(en3d)
st3d = st3d - 1  # Convert to 0 based indexing

# Get size of subdomain of this processor and start and end of the subdomain including ghost cells
sz3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
st3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
en3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_sz3dg(sz3dg)
gp.get_st3dg(st3dg)
gp.get_en3dg(en3dg)
st3dg = st3dg - 1  # Convert to 0 based indexing

# print "rank ", comm.rank, ": sz3D  = ", sz3d,  ", st3D  = ", st3d,  ", en3D  = ", en3d
# print "rank ", comm.rank, ": sz3Dg = ", sz3dg, ", st3Dg = ", st3dg, ", en3Dg = ", en3dg

array = numpy.zeros( (sz3dg[0], sz3dg[1], sz3dg[2]), dtype=numpy.float64, order='F' )
for k in range(nghosts[2],sz3dg[2]-nghosts[2]):
    for j in range(nghosts[1],sz3dg[1]-nghosts[1]):
        for i in range(nghosts[0],sz3dg[0]-nghosts[0]):
            array[i,j,k] = (i+st3dg[0]) + (j+st3dg[1])*nx + (k+st3dg[2])*nx*ny


# ----------- X direction halo communication -----------
gp.fill_halo_x( array )

mycorrect = numpy.array([True])
for k in range(nghosts[2],sz3dg[2]-nghosts[2]):
    for j in range(nghosts[1],sz3dg[1]-nghosts[1]):
        for i in range(sz3dg[0]):
            if numpy.absolute(array[i,j,k] - ( (i+st3dg[0]+nx)%nx + ((j+st3dg[1]+ny)%ny)*nx + ((k+st3dg[2]+nz)%nz)*nx*ny ) ) > 1.e-15:
                mycorrect[0] = False

correct = numpy.array([False])
comm.Reduce(mycorrect, correct, op=MPI.LAND, root=0)

if comm.rank == 0:
    if correct[0]:
        print "Halo communication in X works! :)"
    else:
        print "ERROR: Halo communication in X fails! :("
# ------------------------------------------------------

# ----------- Y direction halo communication -----------
gp.fill_halo_y( array )

mycorrect = numpy.array([True])
for k in range(nghosts[2],sz3dg[2]-nghosts[2]):
    for j in range(sz3dg[1]):
        for i in range(sz3dg[0]):
            if numpy.absolute(array[i,j,k] - ( (i+st3dg[0]+nx)%nx + ((j+st3dg[1]+ny)%ny)*nx + ((k+st3dg[2]+nz)%nz)*nx*ny ) ) > 1.e-15:
                mycorrect[0] = False

correct = numpy.array([False])
comm.Reduce(mycorrect, correct, op=MPI.LAND, root=0)

if comm.rank == 0:
    if correct[0]:
        print "Halo communication in Y works! :)"
    else:
        print "ERROR: Halo communication in Y fails! :("
# ------------------------------------------------------

# ----------- Z direction halo communication -----------
gp.fill_halo_z( array )

mycorrect = numpy.array([True])
for k in range(sz3dg[2]):
    for j in range(sz3dg[1]):
        for i in range(sz3dg[0]):
            if numpy.absolute(array[i,j,k] - ( (i+st3dg[0]+nx)%nx + ((j+st3dg[1]+ny)%ny)*nx + ((k+st3dg[2]+nz)%nz)*nx*ny ) ) > 1.e-15:
                mycorrect[0] = False

correct = numpy.array([False])
comm.Reduce(mycorrect, correct, op=MPI.LAND, root=0)

if comm.rank == 0:
    if correct[0]:
        print "Halo communication in Z works! :)"
    else:
        print "ERROR: Halo communication in Z fails! :("
# ------------------------------------------------------

