#!/bin/bash
#------------------------------------------------------------------------------
#
#  build CWBSUM (GFS and RSM) forecast model
# 
# 1. select targe machine(fx10, fx100, or pcc) 
# 2. run this building up script
#   ./build.sh [MACHINE]
#                                                                Deng-Shun Chen
#                                                                   2020-04-28
#-------------
# add feature:                                        
#   1. for RSM-IO                                                  CHEN,YING-JU
#      use "modulefile.tcogfs.${MACHINE}.rsm" files                  2021-07-30
#------------------------------------------------------------------------------
#
# MACHINE : fx10, fx100, pcc,fx1000
#
if [ $# == 1 ] ; then
  export MACHINE=${1}
else
  echo "usage: $0 [MACHINE]"
  exit
fi
machines='fx1000 fx100 fx10 pcc'
[[ $machines =~ (^|[[:space:]])$MACHINE($|[[:space:]]) ]] && known='True' || known='False'
if [ "${known}" == 'True' ] ; then
  echo "${HOSTNAME} : Build ${MACHINE} executable"
else
  echo "Fatal Error : $0: Unknown machine --> ${MACHINE}" ; exit
fi

if [ "${MACHINE}" == 'fx10' ] ; then 
  hostnames='login07 login08 login05 login06'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False' 
  if [ "${known}" == 'False' ] ; then
    echo "Fatal Error : Build ${MACHINE} executable, please move to login07/08 for inside HPC, login05/06 for outside HPC !" 
    exit
  fi
fi
if [ "${MACHINE}" == 'fx100' ] ; then
  hostnames='login11 login12 login15 login16'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False' 
  if [ "${known}" == 'False' ] ; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login11/12 for inside HPC, login15/16 for outside HPC !" 
   exit
  fi
fi
if [ "${MACHINE}" == 'fx1000' ] ; then
  hostnames='h6ln12 h6ln13 h6ln15 h6ln16 h6ln17 h6ln18 h6ln19 h6ln23'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False'
  if [ "${known}" == 'False' ] ; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login12/13/15/16/17/18/19 for inside HPC, login23 for outside HPC !"
   exit
  fi
fi

set -x

export MDIR=$(pwd)
#
# build RSM
#
cd rsm;
./build.sh ${MACHINE} 5km
cd $MDIR
#
# build TIMCOM regional
#
cd rsm/src/timcom.fd
./build.sh ${MACHINE}
cd $MDIR
#
# build TCo639L72
#
cd tco639l72;
./build.sh ${MACHINE}
cd $MDIR
#
# build TIMCOM glb
#
cd timcom_glb;
./build.sh ${MACHINE}
cd $MDIR

# build OMIP
#cd timcom_glb;
#./build_omip.sh ${MACHINE}
#cd $MDIR


cp /data/common/gfs/GEPSv3_lib/config/conf* ${MDIR}/jobs/

set +x
echo " "
echo " "
echo " "
echo "--------------------------------------------------"
echo "----- build GEPSv3 -------------------------------"
echo "----- check tco639l72/src/libtco.a ------"
echo "----- check rsm/run/compile_4cpl.log ------------------"
echo "----- check rsm/run/exe/librsm4cpl.a-----------"
echo "----- check rsm/src/timcom.fd/libtimcom_tai.a"
echo "----- check timcom_glb/src_lib/libtimcom.a--------------------"
echo " "
echo " "
echo " "
