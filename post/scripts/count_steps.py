#!/bin/python
import sys
import os
from datetime import datetime
from decimal import Decimal


# Get the nodes requested from the cobaltlog 
def get_steps(fname,verbose=True):

    # Get tje lasst chunk of data
    with open('{}'.format(fname),'r') as file:
        data_raw = file.read().splitlines()
    data_raw = data_raw[:-15]

    # Search for step number from bottom up 
    count = 0
    for line in data_raw[::-1]:
        if 'Step' in line: 
            line = list(filter(None,line.split('=')))
            count = max(count, int(line[-1])) 

    if verbose:
        print("File {}:".format(fname))
        print("\tsteps = {}".format(count))
    return count 


if __name__ == '__main__':
    if len(sys.argv)<2:
        print("Returns the steps from all output files")
        print("Usage: \n\t{} dirname [verbose (True)]".format(sys.argv[0]))
        sys.exit()
    elif len(sys.argv)==3:
        verbosity=False
    else:
        verbosity = True
        
    dirname = sys.argv[1]

    steps = 0;
    for f in os.listdir(dirname):
        if f.endswith(".output"):
            steps = max(steps, get_steps(dirname+f,verbose=verbosity))
    print("Total steps = {}".format(steps))
