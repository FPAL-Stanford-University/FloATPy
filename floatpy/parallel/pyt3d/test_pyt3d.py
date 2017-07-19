from mpi4py import MPI
import numpy

import pyt3d

#------------------------------------
# Set the MPI communicator
#-------------------------------------
comm = MPI.COMM_WORLD
fcomm = MPI.COMM_WORLD.py2f()

nx = 8
ny = 16
nz = 32

px = 2
py = 2
pz = 2

if (px*py*pz != comm.size):
    print "This program needs to be run with %d processes. Rerun with the correct options." %(px*py*pz)

periodic = numpy.array([True, True, True])
nghosts  = numpy.array([1,1,1], dtype=numpy.int32, order='F')
fail = numpy.array([True])
reorder = True

gp = pyt3d.t3dmod.t3d(fcomm, nx, ny, nz, px, py, pz, periodic, reorder, fail, nghosts=nghosts, createcrosscommunicators=True)

sz3d = numpy.zeros(3, dtype=numpy.int32, order='F')
st3d = numpy.zeros(3, dtype=numpy.int32, order='F')
en3d = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_sz3d(sz3d)
gp.get_st3d(st3d)
gp.get_en3d(en3d)

print "rank ", comm.rank, ": sz3D = ", sz3d, ", st3D = ", st3d, ", en3D = ", en3d

sz3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
st3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
en3dg = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_sz3dg(sz3dg)
gp.get_st3dg(st3dg)
gp.get_en3dg(en3dg)

array = numpy.zeros( (sz3dg[0], sz3d[1], sz3d[2]), dtype=numpy.float64, order='F' )
for k in range(sz3d[2]):
    for j in range(sz3d[1]):  
        for i in range(nghosts[0],sz3dg[0]-nghosts[0]):
            array[i,j,k] = (i+st3dg[0]-1) + (j+st3d[1]-1)*nx + (k+st3d[2]-1)*nx*ny

gp.fill_halo_x( array )

mycorrect = numpy.array([True])
for k in range(sz3d[2]):
    for j in range(sz3d[1]):  
        for i in range(sz3dg[0]):
            if numpy.absolute(array[i,j,k] - ( (i+st3dg[0]-1+nx)%nx + ((j+st3d[1]-1+ny)%ny)*nx + ((k+st3d[2]-1+nz)%nz)*nx*ny ) ) > 1.e-15:
                mycorrect[0] = False

correct = numpy.array([False])
comm.Reduce(mycorrect, correct, op=MPI.LAND, root=0)

if comm.rank == 0:
    if correct[0]:
        print "Halo communication in X works! :)"
    else:
        print "Halo communication in X fails! :("

