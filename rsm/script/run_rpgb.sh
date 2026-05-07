#!/bin/bash
set -x
##
pwd
PROG=rpgbnawips
hxs=${1}

mkdir -p dir_rpgb${hxs}
cd dir_rpgb${hxs}

h=` expr $hxs + 0 `

if [ $h -lt 10 ] ; then h=0$h ; fi
r_pgbf=r_pgb.f$h
##############################
# post script
##############################

rm -f fort.[0-9]* 2>/dev/null
ln -fs ../r_sigf$h     fort.11
if [ $h -eq 00 ];  then
ln -fs ../r_flx.f$h     fort.21
else
ln -fs ../r_flxf$h     fort.21
fi
ln -fs $r_pgbf       fort.51
ln -fs ctlprs        fort.61
ln -fs ctlslr        fort.62
ln -fs ctldlr        fort.63
#if [ $NEWSLM = 1 ] ; then
#ln -fs $SLMFILE      fort.71
#fi
rm -f $PROG.x
#ln -fs $UTLDIR/$PROG.x $PROG.x
#ln -sf /nwpr/gfs/xb119/GEPS/exe/rpgbnawips.x $PROG.x 
ln -fs /nwpr/gfs/xb169/code/rsm/utl/rpgbnawips.x $PROG.x
./$PROG.x <../rpgb.parm >../stdout.rpgb$h 2>&1
if [ $? -ne 0 ] ; then
  exit 8
fi
#
sed -n '/vars/,/endvars/{/vars/b;/endvars/b;p}' ctlprs |wc -l > wc.out
read n <wc.out && rm wc.out
echo "s|DATAFILE|^$r_pgbf|g" >inp
echo "s|MAPFILE|^$r_pgbf.map|g" >>inp
echo "s|TOTALNUM|$n|g" >>inp
sed -f inp ctlprs >$r_pgbf.ctlprs
rm -f fort.[0-9]* 2>/dev/null


#cp -a r_pgb.f* ../.
mv    r_pgb.f* ../.
#mv    ctlslr ../.
#mv    ctldlr ../.
cd ../
rm -rf dir_rpgb${hxs}

exit

