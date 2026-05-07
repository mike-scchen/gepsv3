#!/bin/bash
set -x
##
pwd
PROG=rpgbnawips
INCMIN=5
hxs=0
FSTM=5
FEND=360
dhd=`expr $FEND / 60`
lmd=$[ $FEND % 60 ]
lmm=59
hd=` expr $hxs + 0 `
md=` expr $FSTM + 0 `

while [ $hd -le $dhd ]; do
if [ $hd -lt 10 ] ; then hh="0${hd}" ; else hh="${hd}" ; fi
if [ $hd -eq $dhd ] ; then lmm=$lmd ;fi
while [ $md -le $lmm ]; do
if [ $md -lt 10 ] ; then mm="0${md}" ; else mm="${md}" ; fi
#r_flxf00:09:00
h="${hh}:${mm}:00" 
echo $h
#r_pgbf=r_pgb.f$h
r_pgbf=r_pgb.f${h}
##############################
# post script
##############################
#UTLDIR=/home/rsmop/TCO_RSM_branch/tco639l72/RSM/utl
UTLDIR=/work1/yj8taiwania/model/tco639l72/RSM/utl

rm -f fort.[0-9]* 2>/dev/null
#ln -fs ../r_sigf$h     fort.11
ln -fs r_sigf$h     fort.11
#if [ $hd -eq 00 ];  then
##ln -fs ../r_flx.f$h     fort.21
#ln -fs r_flx.f$h     fort.21
#else
#ln -fs ../r_flxf$h     fort.21
ln -fs r_flxf$h     fort.21
#fi
#ln -fs r_sig.f$h     fort.11
#ln -fs r_flx.f$h     fort.21
ln -fs $r_pgbf       fort.51
ln -fs ctlprs        fort.61
ln -fs ctlslr        fort.62
ln -fs ctldlr        fort.63
#cp -a  ../ctlprs        fort.61
#cp -a  ../ctlslr        fort.62
#cp -a  ../ctldlr        fort.63
#if [ $NEWSLM = 1 ] ; then
#ln -fs $SLMFILE      fort.71
#fi
rm -f $PROG.x
ln -fs $UTLDIR/$PROG.x $PROG.x
./$PROG.x <rpgb.parm >stdout.rpgb$h 2>&1
#./$PROG.x <../rpgb.parm >stdout.rpgb$h 2>&1
if [ $? -ne 0 ] ; then
  exit 8
fi
#
#grep ',0' ctlprs | wc -l >wc.out
grep ',0' fort.61 | wc -l >wc.out
read n <wc.out && rm wc.out
echo "s|DATAFILE|^$r_pgbf|g" >inp
echo "s|MAPFILE|^$r_pgbf.map|g" >>inp
echo "s|TOTALNUM|$n|g" >>inp
sed -f inp ctlprs >$r_pgbf.ctlprs
#sed -f inp fort.61 >$r_pgbf.ctlprs
#rm -f fort.[0-9]* 2>/dev/null

md=`expr $md + $INCMIN`
if [ $md -eq 60 ] ; then md=100 ; fi
done
md=0
hd=`expr $hd + 1`
done


#cp -a r_pgb.f* ../.
#mv    r_pgb.f* ../.
#cd ../

exit

