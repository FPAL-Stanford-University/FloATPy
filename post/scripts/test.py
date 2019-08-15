from mpi4py import MPI
import numpy as np
import os
import sys
    
if __name__ == '__main__':
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    procs = comm.Get_size()

    sendbuf=[]
    root=0
    if comm.rank==0:
        m=np.array(range(comm.size*comm.size),dtype=float)
        m.shape=(comm.size,comm.size)
        print(m)
        sendbuf=m
    v=comm.scatter(sendbuf,root)
    print("I got this array:")
    print(v)
    v=v*v
    recvbuf=comm.gather(v,root)
    if comm.rank==0:
        print np.array(recvbuf)
