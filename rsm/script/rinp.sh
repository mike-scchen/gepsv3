#/bin/sh

set -ex

#
PROG=rinp
dtype=$1
fh=$2
#yj##
#yj#if [ do$dtype = doG2R ] ; then
#yj#   sig2rg=.TRUE.
#yj#   sfc2rg=.TRUE.
#yj#   newmtn=.TRUE.
#yj#   pgb2rg=.FALSE.
#yj#fi
#yj#if [ do$dtype = doP2R ] ; then
#yj#   sig2rg=.FALSE.
#yj#   sfc2rg=.FALSE.
#yj#   newmtn=.TRUE.
#yj#   pgb2rg=.TRUE.
#yj#fi  
#yj#if [ do$dtype = doC2R -o do$dtype = doN2R ] ; then
#yj#   sig2rg=.TRUE.
#yj#   sfc2rg=.TRUE.
#yj#   newmtn=.TRUE.
#yj#   pgb2rg=.FALSE.
#yj#fi
#yj#
#yj#   
#yj#rm -f fort.[0-9]* 2>/dev/null
#yj##
#yj#    echo " Regional input starts:"
#yj#    echo " &NAMRIN                                                "  >rinpparm
#yj#    echo "    SIG2RG=$sig2rg,SFC2RG=$sfc2rg,PERCMTN=$PERCMTN,     " >>rinpparm
#yj#    echo "    NEWSIG=$NEWSIG,NEWMTN=$newmtn,NEWHOR=.FALSE.,       " >>rinpparm
#yj#    echo "    PGB2RG=$pgb2rg,NEWSST=$NEWSST,                      " >>rinpparm
#yj#    echo "    ivs=$IVS,iqvar=$IQVAR,igribversion=$IGRIBVERSION,   " >>rinpparm
#yj#    echo " &END                                                   " >>rinpparm
#yj##
#yj#    cat rsmlocation >> rinpparm
#yj#
#yj#    echo " &namsig                                         " >>rinpparm
#yj##    echo "    fcsttime=-9999,                              " >>rinpparm
#yj#    echo " &end                                            " >>rinpparm
#yj#    echo " &namsfc                                         " >>rinpparm
#yj##    echo "    fcsttime=-9999,                              " >>rinpparm
#yj#    echo " &end                                            " >>rinpparm
#yj#    echo " &namclim                                        " >>rinpparm
#yj#if [ $CLIM = 1 ]
#yj#then
#yj#    echo "    iclim=1,                                    " >>rinpparm
#yj#    echo "     fnmskh=\"$FNMSKH\",                             " >>rinpparm
#yj#    echo "     fnalbc=\"$FNALBC\",                             " >>rinpparm
#yj#    echo "     fnsotc=\"$FNSOTC\",                             " >>rinpparm
#yj#    echo "     fnvegc=\"$FNVEGC\",                             " >>rinpparm
#yj#    echo "     fnvetc=\"$FNVETC\",                             " >>rinpparm
#yj#    echo "     fnzorc=\"$FNZORC\",                             " >>rinpparm
#yj#    echo "    fntg3c=\"$FNTG3C\",                             " >>rinpparm
#yj#fi
#yj#    echo " &end                                            " >>rinpparm

if [ do$dtype = doG2R -o do$dtype = doC2R -o do$dtype = doN2R ] ; then
    ln -fs rb_sigf$fh   fort.11
    ln -fs rb_sfcf$fh   fort.12
fi
if [ do$dtype = doP2R ] ; then
    ln -fs rb_pgbf$fh   input.grib
fi
if [ do$NEWSST = do.TRUE. ] ; then
    ln -fs rb_sstf$fh   sst.grib
fi
     
#
    ln -fs rmtnslm             fort.13
    ln -fs rmtnoss             fort.14
    ln -fs ${FIXDIR}/siglevel.l${LEVR}.txt       fort.15
if [ do$dtype = doP2R ] ; then
    ln -fs ${FIXDIR}/pgblevel.l${LEVS}.txt       fort.17
fi

if [ $CLIM = 1 ]
then
    ln -fs ${FIXDIR}/$FNMSKH .
    ln -fs ${FIXDIR}/$FNALBC .
    ln -fs ${FIXDIR}/$FNSOTC .
    ln -fs ${FIXDIR}/$FNVEGC .
    ln -fs ${FIXDIR}/$FNVETC .
    ln -fs ${FIXDIR}/$FNZORC .
    ln -fs ${FIXDIR}/$FNTG3C .
fi
#
pwd
    ln -fs r_sigtmp            fort.51
    ln -fs r_sfctmp            fort.52
    ln -fs r_sigi              fort.61
    ln -fs r_sfci              fort.62
#
# (note: If NEWMTN=.TRUE. --> output=61, 62)
# (note: If NEWMTN=.FALSE. --> output=51, 52)
#
rm -f $PROG.x
#ln -fs $SRCDIR/rsm_$PROG.fd/$PROG.x $PROG.x
ln -fs $EXPEXE/$PROG.x $PROG.x
./$PROG.x <rinpparm > stdout.rinp$fh  || exit
#
NEWMTN=YES
if [ $NEWMTN = NO ]
then
    cp r_sigtmp r_sigi || exit 2
    cp r_sfctmp r_sfci || exit 3
fi
  rm -f r_sigtmp r_sfctmp
#
  cp r_sigi r_sigitdt
  cp r_sigi  r_sig.f${fh}
  cp r_sfci  r_sfc.f${fh}
#== yj added # ==start
  cp r_sigi rinp_r_sigitdt
  cp r_sigi rinp_r_sig.f${fh}
  cp r_sfci rinp_r_sfc.f${fh}
ls -latr r_sigi r_sfci
echo "=========check the output file==========rinp"
#== yj added # ==end
rm -f fort.[0-9]* 2>/dev/null

echo "========== exit the rinp.sh =========="
