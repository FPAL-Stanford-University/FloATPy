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
fail = numpy.array([True])
reorder = True

# Create the t3d object
gp = pyt3d.t3dmod.t3d(fcomm, nx, ny, nz, periodic)

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

# Get size of subdomain of this processor and start and end of the subdomain
sz3d = numpy.zeros(3, dtype=numpy.int32, order='F')
st3d = numpy.zeros(3, dtype=numpy.int32, order='F')
en3d = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_sz3d(sz3d)
gp.get_st3d(st3d)
gp.get_en3d(en3d)
st3d = st3d - 1  # Convert to 0 based indexing

# Get size of subdomain of this processor and start and end of the subdomain
szx = numpy.zeros(3, dtype=numpy.int32, order='F')
stx = numpy.zeros(3, dtype=numpy.int32, order='F')
enx = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_szx(szx)
gp.get_stx(stx)
gp.get_enx(enx)
stx = stx - 1  # Convert to 0 based indexing

# Get size of subdomain of this processor and start and end of the subdomain
szy = numpy.zeros(3, dtype=numpy.int32, order='F')
sty = numpy.zeros(3, dtype=numpy.int32, order='F')
eny = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_szy(szy)
gp.get_sty(sty)
gp.get_eny(eny)
sty = sty - 1  # Convert to 0 based indexing

# Get size of subdomain of this processor and start and end of the subdomain
szz = numpy.zeros(3, dtype=numpy.int32, order='F')
stz = numpy.zeros(3, dtype=numpy.int32, order='F')
enz = numpy.zeros(3, dtype=numpy.int32, order='F')
gp.get_szz(szz)
gp.get_stz(stz)
gp.get_enz(enz)
stz = stz - 1  # Convert to 0 based indexing

print "rank ", comm.rank, ": sz3D  = ", sz3d,  ", st3D  = ", st3d,  ", en3D  = ", en3d
print "rank ", comm.rank, ": szx   = ", szx ,  ", stx   = ", stx ,  ", enx   = ", enx 

array = numpy.zeros( (sz3d[0], sz3d[1], sz3d[2]), dtype=numpy.float64, order='F' )
for k in range(sz3d[2]):
    for j in range(sz3d[1]):
        for i in range(sz3d[0]):
            array[i,j,k] = (i+st3d[0]) + (j+st3d[1])*nx + (k+st3d[2])*nx*ny


print comm.rank, ": array = ", array[:,0,0]

# ----------- X direction halo communication -----------
array_x = numpy.zeros( (szx[0], szx[1], szx[2]), dtype=numpy.float64, order='F' )
gp.transpose_3d_to_x( array, array_x )

print comm.rank, ": array_x = ", array_x[:,0,0]

mycorrect = numpy.array([True])
for k in range(szx[2]):
    for j in range(szx[1]):
        for i in range(szx[0]):
            if numpy.absolute(array_x[i,j,k] - ( (i+stx[0]) + ((j+stx[1]))*nx + ((k+stx[2]))*nx*ny ) ) > 1.e-15:
                mycorrect[0] = False

correct = numpy.array([False])
comm.Reduce(mycorrect, correct, op=MPI.LAND, root=0)

if comm.rank == 0:
    if correct[0]:
        print "Halo communication in X works! :)"
    else:
        print "ERROR: Halo communication in X fails! :("
# ------------------------------------------------------
