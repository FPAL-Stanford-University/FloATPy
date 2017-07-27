#!/bin/bash

if [[ -z "${MPIEXEC}" ]]; then
    MPIEXEC='mpiexec'
else
    MPIEXEC="${MPIEXEC}"
fi

# Get the directory of this file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for testfile in ${DIR}/test_*.py; do
    python ${testfile}
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

for testfile in ${DIR}/mpitest_*.py; do
    ${MPIEXEC} -n 4 python ${testfile}
done
