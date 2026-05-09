#!/bin/bash
# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

set -e

# Error handling
trap 'echo "Error occurred at line $LINENO"; exit 1' ERR

# Ensure a MACHINE argument is provided
if [ $# -lt 1 ]; then
  echo "usage: $0 [MACHINE] [TARGET_MODE]"
  echo "  [MACHINE]   : fx10, fx100, pcc, fx1000, a100"
  echo "  [TARGET_MODE]: OMIP (optional, default is empty, meaning "LIB")"
  exit 1
fi

MACHINE=$1
TARGET_MODE=${2:-}  # Default to empty string if not provided

# Function to check if a machine is valid
check_machine_validity() {
  local machine=$1
  local known_machines=('fx1000' 'gpu')

  if [[ ! " ${known_machines[@]} " =~ " ${machine} " ]]; then
    echo "Fatal Error : Unknown machine --> ${machine}"
    exit 1
  fi

  # Check if building on fx1000 login node
# You may enable the check / modify the hostname
  if [ "$machine" == "fx1000" ] && [[ ! $HOSTNAME =~ h6ln[0-9][0-9] ]]; then
    echo "Fatal Error: Please move to a valid login node for HPC."
    exit 1
  fi
}

check_machine_validity "$MACHINE"

# Set USE_OMIP based on TARGET_MODE (default to "LIB" if empty)
if [ -z "$TARGET_MODE" ]; then
  USE_OMIP="OFF"  # Default is "LIB", so USE_OMIP="OFF"
elif [ "$TARGET_MODE" == "OMIP" ]; then
  USE_OMIP="ON"
else
  echo "Fatal Error: Invalid TARGET_MODE value. It should be either empty (default 'LIB') or 'OMIP'."
  exit 1
fi

# Set machine-specific configurations

if [ "$MACHINE" == 'gpu' ]; then
  echo "Loading NVIDIA HPCSDK for GPU..."
  . $MODULESHOME/init/bash
  #. $LMOD_ROOT/lmod/init/bash
  module use ../modulefiles
  module load modulefile.tcogfs.gpu
  module list
  USE_ARM="OFF"
  USE_GPU="ON"
  USE_CLM="OFF"  # =ON(with clm_r, restart)
  USE_ICE="OFF"  # =ON(with cpl_cice, GPU does NOT support this for now)
elif [ "$MACHINE" == 'fx1000' ]; then
  echo "Loading Fujitsu for fx1000..."
  . $MODULESHOME/init/bash
  #. $LMOD_ROOT/lmod/init/bash
  module use ../modulefiles
  module load modulefile.tcogfs.fx1000
  module list
  USE_ARM="ON"
  USE_GPU="OFF"
  USE_CLM="OFF"  # =ON(with clm_r, restart)
  USE_ICE="ON"  # =ON(with cpl_cice)
fi

# Run cmake and build
echo "Building with CMake..."

cmake -B build_${MACHINE} -S . -DUSE_OMIP=${USE_OMIP} -DUSE_ARM=${USE_ARM} -DUSE_GPU=${USE_GPU} -DUSE_CLM=${USE_CLM} -DUSE_ICE=${USE_ICE}
cmake --build build_${MACHINE} -j $(($(nproc) / 2)) -v
