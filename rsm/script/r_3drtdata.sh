#! /bin/sh
#
set -ex
#
PROG=r3drtdata
WORKDIR=$1
MACHINE=$2
WORK=$3
#. ${DISK}/run/configure
#to get the larger boundary of lat. lon., set the smaller distance of lat. and Earth radius.
#further more, extends 2 grid points to calculate the new boundary.
export slat=`echo "scale=2; $RCENLAT-(($RBTMGRD+2)*$RDELY/1000./109.)" | bc`
export nlat=`echo "scale=2; $RCENLAT+(($JGRD-$RBTMGRD+2)*$RDELY/1000./109.)" | bc`
export alat=`echo "scale=8; (6300.*2*a(1)*4*c(($slat+$nlat)/2.*a(1)*4/180))/360." | bc -l`
export wlon=`echo "scale=2; $RCENLON-((($RLFTGRD+2)*$RDELX/1000.)/$alat)" | bc`
export elon=`echo "scale=2; $RCENLON+((($IGRD-$RLFTGRD+2)*$RDELX/1000.)/$alat)" | bc`
#
#     namelist /naminfo/ TS_LAT, TE_LAT, TS_LON, TE_LON
#
echo " &naminfo" > info_3drt.parm
echo " TS_LAT=$slat," >> info_3drt.parm
echo " TE_LAT=$nlat," >> info_3drt.parm
echo " TS_LON=$wlon," >> info_3drt.parm
echo " TE_LON=$elon," >> info_3drt.parm
echo " condir='$WORKDIR'," >> info_3drt.parm
echo " FILE_NAME='nc3.nc'," >> info_3drt.parm
echo " &END" >> info_3drt.parm

# get the land sea mask of the RSM/MSM domain from "RMTN"
if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
ln -fs ${WORK}/rmtn_init/rmtnslm         fort.23
else
export FORT23=${WORK}/rmtn_init/rmtnslm
fi

if [ ${RUNRSMTERRAIN} = yes ]; then
    ln -fsv "${TEMP}/rsm_terrain_init/rsm-lsm.bin" fort.23
fi

#
ln -fs  ${FIXDIR}/Param_3DRT/* .
cat rsmlocationinit info_3drt.parm > r3drtdata.parm
#
ls -l
pwd
##
ln -fs $EXPEXE/${PROG}.x $PROG.x
./$PROG.x <r3drtdata.parm >stdout.r3drtdata  || exit
cat stdout.r3drtdata
rm -f fort.[0-9]* 2>/dev/null

################  end of file  ##################

#>>----3DRT data information----->>      
