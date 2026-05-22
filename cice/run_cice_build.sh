#!/bin/bash

MDIR=$(pwd)
echo ${MDIR}
CICE_VERSION=${MDIR}/CICE-CICE6.5.0
CASE_PATH=${MDIR}/cpl_lib_musoac

cd ${CICE_VERSION}
#./cice.setup -c $CASE_PATH -g 0p25 -m fx1000 -e fujitsu
./cice.setup -c $CASE_PATH -g 0p25 -m f1 -e nvhpc

cd ${CASE_PATH}
./cice.build makdep

./cice.build depends

./cice.build libcice.a

ln -sf compile include
