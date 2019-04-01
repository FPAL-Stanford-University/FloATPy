#!/bin/python
import sys
import os
from datetime import datetime
from decimal import Decimal


# A class for accessing vars given a namelist
class inputs:
    def __init__(self,dirname,verbose=False):
        self.dirname = dirname

# Get the nodes requested from the cobaltlog 
def get_nnodes(line):
    original = line
    try:
        line = list(filter(None,line.split('-n')))
        line = list(filter(None,line[1].split(' ')))
        n = line[0]
    except Exception as e:
        print("Error getting N nodes for line")
        print("\t{}".format(original))
        n = None
    return float(n) 

# Get the time from the cobaltlog 
def get_time(line):
    original = line
    try:
        line = list(filter(None,line.split('+')))
        line = list(filter(None,line[0].split(' ')))
        time = line[3]
    except Exception as e:
        print("Error getting day, time for line")
        print("\t{}".format(original))
        time = "00:00:00"
    return time  


# Gets Mc, rr, Re from namelist
def get_cost(fname,verbose=False):
    with open('{}'.format(fname),'r') as file:
        data_raw = file.read().splitlines()

    # Get the number of nodes (1core/node) requested
    nodes = get_nnodes(data_raw[1])

    # Second to last line is the start time (after booting)
    # Last line is the finish (check if completed normally)
    t1 = get_time(data_raw[-2])
    t2 = get_time(data_raw[-1])
    FMT = '%H:%M:%S'
    runtime = datetime.strptime(t2,FMT) - datetime.strptime(t1,FMT)
    if runtime.days < 0:
        print("Runtime passed midnight,adjusting")
        run1 = datetime.strptime("23:59:59",FMT) - datetime.strptime(t1,FMT)
        run2 = datetime.strptime(t2,FMT) - datetime.strptime("00:00:00",FMT)
        runtime = run1+run2

    # Convert to decimal hours
    hrs = runtime.total_seconds()/360.

    cost = hrs*nodes
    if verbose:
        print("{}:".format(data_raw[0]))
        print("\t nodes = {}".format(nodes))
        print("\t t1 = {}".format(t1))
        print("\t t2 = {}".format(t2))
        print("\t cost = %.2E core-hours"%Decimal(cost))
    
    return cost


if __name__ == '__main__':
    if len(sys.argv)<2:
        print("Returns the cost in node-hrs (core-hrs) for all cobaltlogs")
        print("Usage: \n\t{} dirname".format(sys.argv[0]))
        sys.exit()
        
    dirname = sys.argv[1]

    cost = 0;
    for f in os.listdir(dirname):
        if f.endswith(".cobaltlog"):
           cost += get_cost(dirname+f,verbose=True)
    print("Total cost = %.2E core-hours"%Decimal(cost))
           


