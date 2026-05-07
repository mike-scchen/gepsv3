#!/bin/ksh
#PJM -L "node=1:noncont"
#PJM -L elapse=00:30:00
#PJM -L node-mem=unlimited
#PJM -L rscgrp=large-x
#PJM --no-stging
#PJM --mpi "proc=1"
#PJM -j
##
set -x


PWDIR=/IFS6/data2/datagfs/xb127/TCo_RSM_branch/tco639l72
############################## GFS START ###################################

SDATE=19090100
echo $SDATE

GFSWRK=/nwpr/gfs/xb127/data2/TCo_RSM_branch/RSM_arc/wrk_RSM_${SDATE}_NoahLSM

cd $GFSWRK

/usr/bin/time ${GFSWRK}/run_extrsfc.sh
echo "Finish RSM post processor - extract sfc field"

cp -r  ${PWDIR}/plot  $GFSWRK/.
/usr/bin/time ${GFSWRK}/plot/drawsfc.sh   $GFSWRK
/usr/bin/time ${GFSWRK}/plot/drawatm.sh   $GFSWRK
echo "Finish MSM 1km surface field plotting"

