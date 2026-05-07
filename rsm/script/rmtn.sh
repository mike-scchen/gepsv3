#! /bin/sh
#
set -ex
#
PROG=rmtn
mtnres=$1
type=$2
WORKDIR=$3
#yj# type=init or type=base
#
rm -f fort.[0-9]* top* 2>/dev/null
if [ $mtnres -ne 30 ] ; then   # not GTOPO30
  if [ $mtnres -ne 1 ] ; then  # not SRTM1ARC
   MTNDATA=top${mtnres}m
   MTNDIR=${FIXDIR}
   MTN_AVG=$MTNDIR/${MTNDATA}_avg.20i4.asc
   MTN_VAR=$MTNDIR/${MTNDATA}_var.20i4.asc
   MTN_MAX=$MTNDIR/${MTNDATA}_max.20i4.asc
   MTN_SLM=$MTNDIR/${MTNDATA}_slm.80i1.asc
   ##cp $MTNDIR/${MTNDATA}* . || exit
   ##uncompress ${MTNDATA}*
   cp rsmlocation${type} rmtn.parm
  fi
fi
#
#GTOPO30
if [ $mtnres -eq 30 ] ; then   
   MTN_AVG=dummy1
   MTN_VAR=dummy2
   MTN_MAX=dummy3
   MTN_SLM=dummy4
if [ $RMTN_STUB_ONLY != yes ]; then
   ln -fs ${FIXDIR}/GTOPO30/*.tar.gz .
   for i in `ls *.gz`
   do 
     gunzip -c $i | tar xvf -
   done
   rm *.tar.gz
fi
   MTNDIR=${RUNDIR}
##   MTNDIR=${RUNDIR}
   MTNDIR=${WORKDIR}
   echo " &namcondir" >condir.parm
   echo " condir='$MTNDIR'," >>condir.parm
   echo " &END" >>condir.parm
   cat rsmlocation${type} condir.parm >rmtn.parm
fi
#SRTM1ARC
if [ $mtnres -eq 1 ] ; then   
   MTN_AVG=dummy1
   MTN_VAR=dummy2
   MTN_MAX=dummy3
   MTN_SLM=dummy4
   ln -fs ${FIXDIR}/SRTM1ARC_TWN/*.zip .
   for i in `ls *.zip`
   do 
     unzip  -n ${FIXDIR}/SRTM1ARC_TWN/$i -d ${WORKDIR}
#    ssh login07 unzip  -n ${FIXDIR}/SRTM1ARC_TWN/$i -d ${WORKDIR}
   done
   rm *.zip
##   MTNDIR=${RUNDIR}
   MTNDIR=${WORKDIR}
   echo " &namcondir" >condir.parm
   echo " condir='$MTNDIR'," >>condir.parm
   echo " &END" >>condir.parm
   cat rsmlocation${type} condir.parm >rmtn.parm
fi
#
rm -f fort.[0-9]* 2>/dev/null
#
if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
 ln -fs $MTN_AVG       fort.11
 ln -fs $MTN_VAR       fort.12
 ln -fs $MTN_MAX       fort.13
 ln -fs $MTN_SLM       fort.14
else
 export FORT11=$MTN_AVG
 export FORT12=$MTN_VAR
 export FORT13=$MTN_MAX
 export FORT14=$MTN_SLM
fi
#
#
ls -l
#
if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
 ln -fs rmtnslm     fort.51
 ln -fs rmtnoro     fort.52
 ln -fs rmtnvar     fort.53
 ln -fs rmtnors     fort.54
 ln -fs rmtnoss     fort.55
else
 export FORT51=rmtnslm
 export FORT52=rmtnoro
 export FORT53=rmtnvar
 export FORT54=rmtnors
 export FORT55=rmtnoss
fi
#
rm -f $PROG.x
#ln -fs $SRCDIR/rsm_$PROG.fd/$PROG.x $PROG.x
ln -fs $EXPEXE/${PROG}_${type}.x $PROG.x
./$PROG.x <rmtn.parm >stdout.rmtn  || exit
cat stdout.rmtn
rm -f fort.[0-9]* 2>/dev/null

