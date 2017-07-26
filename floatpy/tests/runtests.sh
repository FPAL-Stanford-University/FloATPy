#!/bin/bash

if [[ -z "${MPIEXEC}" ]]; then
    MPIEXEC='mpiexec'
else
    MPIEXEC="${MPIEXEC}"
fi

for testfile in test_*.py; do
    python ${testfile}
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

for testfile in mpitest_*.py; do
    ${MPIEXEC} -n 4 python ${testfile}
done
