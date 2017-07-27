#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`

case "$os" in
    Darwin)
		;;

    Linux)
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac

cd ${HOME}
git clone https://github.com/jameskermode/f90wrap
cd f90wrap
python setup.py install
