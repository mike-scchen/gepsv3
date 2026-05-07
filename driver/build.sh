#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS forecast model
#
# 使用方式（僅支援以下三種）：
#   ./build.sh a100 2cpl
#   ./build.sh fx1000 2cpl
#   ./build.sh fx1000 4cpl
#
#------------------------------------------------------------------------------

# 檢查參數
if [ $# -ne 3 ]; then
  echo "Usage: $0 [MACHINE] [CPL_MODE] [TCo???]"
  echo "Allowed combinations:"
  echo "  a100 2cpl TCo383/TCo199"
  echo "  fx1000 2cpl TCo383/TCo199"
  echo "  fx1000 4cpl TCo383/TCo199"
  echo "  fx1000 2cpl_CICE TCo383/TCo199"
  exit 1
fi

MACHINE=$1
CPL_MODE=$2
RESN=$3
# 驗證合法組合
if [[ "$MACHINE" == "a100" && "$CPL_MODE" == "2cpl" ]]; then
  :
elif [[ "$MACHINE" == "fx1000" && ( "$CPL_MODE" == "2cpl" || "$CPL_MODE" == "4cpl"  || "$CPL_MODE" == "2cpl_CICE" ) ]]; then
  :
else
  echo "Fatal Error: Unsupported combination: $MACHINE $CPL_MODE"
  echo "Allowed: ~a100 2cpl | fx1000 2cpl | fx1000 4cpl | fx1000 2cpl_CICE~ TCo199/TCo383"
  exit 1
fi

# FX1000 特殊登入節點檢查
if [ "$MACHINE" == "fx1000" ]; then
  if [[ ! "$HOSTNAME" =~ ^h6ln[0-9][0-9]$ ]]; then
    echo "Fatal Error: Build fx1000 executable must be done on h6ln?? login node!"
    exit 1
  fi
fi

set -x

# 設定 FRAME
export FRAME=$CPL_MODE
export MACHINE=$MACHINE
export atmres=$RESN
if [ "$MACHINE" == "a100" ]; then
  # A100: 載入 NVIDIA HPC SDK module
  module use /package/x86_64/nvidia/hpc_sdk/modulefiles
  module load nvhpc-hpcx-cuda12/24.11
  echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
else
  # FX1000: 載入專案 module
  MDIR=$(cd ../tco639l72 && pwd)
  . /usr/share/Modules/init/bash
  module purge
  module use ${MDIR}/modulefiles

  module show modulefile.tcogfs.${MACHINE}_${FRAME}
  module load modulefile.tcogfs.${MACHINE}_${FRAME}
  module list
  module unuse ${MDIR}/modulefiles
fi

# 編譯
if [ "$FRAME" == "2cpl_CICE" ]; then
  make -f Makefile.${MACHINE}.${atmres}ice
else
  make -f Makefile.${MACHINE}.${atmres}
fi

if [[ $? -ne 0 ]]; then
  echo "Build failed."
  exit 1
fi

echo "✅ Build successful: $MACHINE $FRAME $atmres"

