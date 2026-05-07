#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS forecast model
#
# ‰Allowed machine / couple mode:
#   ./build.sh gpu 2cpl
#   ./build.sh fx1000 2cpl
#   ./build.sh fx1000 4cpl
#
#------------------------------------------------------------------------------

# Check for number of arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 [MACHINE] [CPL_MODE] [TCo???]"
  echo "Allowed combinations:"
  echo "  gpu 2cpl TCo383/TCo199"
  echo "  fx1000 2cpl TCo383/TCo199"
  echo "  fx1000 4cpl TCo383/TCo199"
  echo "  fx1000 2cpl_CICE TCo383/TCo199"
  exit 1
fi

MACHINE=$1
CPL_MODE=$2
RESN=$3
# Check for legal machine / couple modes
if [[ "$MACHINE" == "gpu" && "$CPL_MODE" == "2cpl" ]]; then
  :
elif [[ "$MACHINE" == "fx1000" && ( "$CPL_MODE" == "2cpl" || "$CPL_MODE" == "4cpl"  || "$CPL_MODE" == "2cpl_CICE" ) ]]; then
  :
else
  echo "Fatal Error: Unsupported combination: $MACHINE $CPL_MODE"
  echo "Allowed: ~a100 2cpl | fx1000 2cpl | fx1000 4cpl | fx1000 2cpl_CICE~ TCo199/TCo383"
  exit 1
fi

# FX1000 
# Check if building on fx1000 login node
# You may enable the check / modify the hostname 
#if [ "$MACHINE" == "fx1000" ]; then
#  if [[ ! "$HOSTNAME" =~ ^h6ln[0-9][0-9]$ ]]; then
#    echo "Fatal Error: Build fx1000 executable must be done on h6ln?? login node!"
#    exit 1
#  fi
#fi

set -x

# Set arguments
export FRAME=$CPL_MODE
export MACHINE=$MACHINE
export atmres=$RESN
if [ "$MACHINE" == "gpu" ]; then
  # GPU: NVIDIA HPC SDK 
  . $MODULESHOME/init/bash
  #. $LMOD_ROOT/lmod/init/bash
  module purge
  module use /package/x86_64/nvidia/hpc_sdk/modulefiles
  module load nvhpc-hpcx-cuda12/24.11
  module use ../modulefiles/modulefile.tcogfs.gpu
  echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
else
  # FX1000: Fujitsu
  . $MODULESHOME/init/bash
  #. $LMOD_ROOT/lmod/init/bash

  module purge
  module use ../modulefiles
  module load modulefile.tcogfs.${MACHINE}_${FRAME}
  module list
fi

# Building
if [ "$FRAME" == "2cpl_CICE" ]; then
  make -f Makefile.${MACHINE}.${atmres}ice
else
  make -f Makefile.${MACHINE}.${atmres}
fi

if [[ $? -ne 0 ]]; then
  echo "Build failed."
  exit 1
fi

echo "Build successful: $MACHINE $FRAME $atmres"

