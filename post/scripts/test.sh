#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))
mpirun -f $COBALT_NODEFILE -n $PROCS python /home/kmatsuno/FloATPy/post/scripts/tke_vs_t.py /projects/ShockInducedMix/ShearLayerData/production/Mc24/rr1/512x512x256/shearlayer_
