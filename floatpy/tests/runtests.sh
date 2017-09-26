#!/bin/bash

if [[ -z "${MPIEXEC}" ]]; then
    MPIEXEC='mpiexec'
else
    MPIEXEC="${MPIEXEC}"
fi
echo "Using '${MPIEXEC}' to run the parallel tests"
echo ""

# Get the directory of this file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for testfile in ${DIR}/test_*.py; do
    echo "python ${testfile}"
    python ${testfile}
    if [ $? -ne 0 ]; then
        exit 1
    fi
    echo ""
done

for testfile in ${DIR}/mpitest_*.py; do
    echo "${MPIEXEC} -n 4 python ${testfile}"
    ${MPIEXEC} -n 4 python ${testfile}
    if [ $? -ne 0 ]; then
        exit 1
    fi
    echo ""
done
