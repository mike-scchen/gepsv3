#!/bin/sh
#
#-----------------------------------------------
# specify starting date and end date
#-----------------------------------------------
#define directories
set -x
syear=`echo $SDATE |cut -c1-4`          ## starting year
smonth=`echo $SDATE |cut -c5-6`         ## starting month
sday=`echo $SDATE |cut -c7-8`           ## starting day
CHOUR=`echo $SDATE |cut -c9-10`         ## starting hour
CDATE=$syear$smonth$sday
dayend=01                ## end date
hourend=$CHOUR           ## end hour
rmm=$smonth
ryyyy=$syear
if [ $sday -gt 15 ]; then
  rmm=`expr $rmm + 1`
  if [ $rmm -gt 12 ]; then
     rmm=`expr $rmm - 12`
     ryyyy=`expr $ryyyy + 1`
  fi
  if [ $rmm -lt 10 ]; then
     rmm=0$rmm
  fi
fi
#------------------------------------------------
# determine length of integration
#------------------------------------------------
yearend=$syear
if [ $sday -gt 15 ]; then
   monthend=`expr $smonth + $LENMON + 1`
else
   monthend=`expr $smonth + $LENMON`
fi
if [ $monthend -gt 12 ]; then
    monthend=`expr $monthend - 12`
    yearend=`expr $yearend + 1`
fi
if [ $monthend -lt 10 ]; then
     monthend=0$monthend
fi
eidate=$yearend$monthend$dayend$hourend
export END_HR=`$UTLDIR/nhour $eidate $CDATE$CHOUR`
export ENDHOUR=${ENDHOUR:-$END_HR}

#-----------------------------------------------
# set running space
#-----------------------------------------------
#yj#BASE_DIR=$TEMP/${SDATE}
#yj#export BASEDIR; BASEDIR=${BASEDIR:-$BASE_DIR}
mkdir -p $RUNDIR
cd $RUNDIR || exit
#
#-----------------------------------------------
# determine model run parameters
#-----------------------------------------------
# move some parameters to utl/rsm_default.option #yj#
#yj#RESVER=$LEVS
#yj#sfc_freq=24
#yj#RSWRHOUR=1
#yj#RLWRHOUR=1
#yj#RDFISEC=0.
#yj#RLXHSEC=${RLXHSEC:-1800.}
#yj#RLXMSEC=${RLXMSEC:-1800.}
#yj#DIFH=3
#yj#DIFM=2
#yj#ISEMIMP=1
#yj#IIMPRLX=1
#yj#IDMPJET=0
#yj#IOUTNHR=1
if [ $RSFC_MERGE = "yes" ]; then
   ISFC_MERGE=1
else
   ISFC_MERGE=0
fi
#
# NO NEED TO CHANGE BELOW THIS!
#
FCSTSEC=`expr $INCHOUR \* 3600`
PRNTSEC=`expr $PRTHOUR \* 3600`
RSWRSEC=`expr $RSWRHOUR \* 3600`
RLWRSEC=`expr $RLWRHOUR \* 3600`
BASESEC=`expr $INCBASE \* 3600`
DIFHSEC=`expr $DIFH \* $DELTAT_REG `
DIFMSEC=`expr $DIFM \* $DELTAT_REG `
if [ $NONHYD = "yes" ]; then
FILTA=0.80
else
FILTA=0.92
fi
IGRD1=`expr $IGRD + 1`
JGRD1=`expr $JGRD + 1`
RSFCSEC=`expr $sfc_freq \* 3600`;
#
#### define fcst parameters ######
 cat >rsmlocation <<EOF
 &NAMLOC
 RPROJ    = $RPROJ,
 RTRUTH   = $RTRUTH,
 RORIENT  = $RORIENT,
 RDELX    = $RDELX,
 RDELY    = $RDELY,
 RCENLAT  = $RCENLAT,
 RCENLON  = $RCENLON,
 RLFTGRD  = $RLFTGRD,
 RBTMGRD  = $RBTMGRD,
 CPROJ    = $CPROJ,
 CTRUTH   = $CTRUTH,
 CORIENT  = $CORIENT,
 CDELX    = $CDELX,
 CDELY    = $CDELY,
 CCENLAT  = $CCENLAT,
 CCENLON  = $CCENLON,
 CLFTGRD  = $CLFTGRD,
 CBTMGRD  = $CBTMGRD,
 CLAT2    = $CLAT2,
 CLAT1    = $CLAT1,
 CLON2    = $CLON2,
 CLON1    = $CLON1,
 &END

EOF
########
cat >rfcstparm <<EOF
 &NAMRSM
 DELTIME=$DELTAT_REG,
 FCSTSEC=$FCSTSEC,
 PRNTSEC=$PRNTSEC,
 RSWRSEC=$RSWRSEC,
 RLWRSEC=$RLWRSEC,
 BASESEC=$BASESEC,
 RDFISEC=$RDFISEC,
 FILTA=$FILTA,
 RLXHSEC=$RLXHSEC,
 RLXMSEC=$RLXMSEC,
 DIFHSEC=$DIFHSEC,
 DIFMSEC=$DIFMSEC,
 RSFCSEC=$RSFCSEC,
 ISEMIMP=$ISEMIMP,
 IIMPRLX=$IIMPRLX,
 IDMPJET=$IDMPJET,
 IMDLPHY=$IMDLPHY,
 IOUTNHR=$IOUTNHR,
 ISFCMRG=$ISFC_MERGE,
 
 &END
EOF
#
  cat rsmlocation >> rfcstparm
  cp rsmlocation $RUNDIR
  cp $STTPRM $RUNDIR/station.parm
#
##########################################################
#
  #
  #  Restart
  #
#yj#        echo 'Restart files existed!!!!!!!'
#yj#         rm fort.*
#yj#         ln -fs r_sigi fort.11
#yj#         $UTLDIR/fhour.x >dte.out
#yj#         read hour month day year FH <dte.out && rm dte.out
#yj#         if [ $FH -lt 10 ];then FH=0$FH;fi
#yj#         export FH
#yj#         export FEND=`expr $FH + $FHMAX`
#yj#         if [ $FEND -gt $ENDHOUR ]; then
#yj#            FEND=$ENDHOUR
#yj#         fi
        export FH=00
        export FEND=`expr $FH + $FHMAX`
        if [ $FEND -gt $ENDHOUR ]; then
           FEND=$ENDHOUR
        fi
#
# RSM INITIAL forecast
#
#     ln -fs $BASEDIR/sigf$CDATE$CHOUR rb_sigf00
#     ln -fs $BASEDIR/sfcf$CDATE$CHOUR rb_sfcf00
if [ -s r_sigi -a -s r_sigitdt -a -s r_sfci ] ; then
  rm r_sigi r_sigitdt r_sfci fort.* r_sig.f* r_sfc.f* 
  cp rinp_r_sig.f00 r_sigi || exit 8
  cp rinp_r_sigitdt r_sigitdt  || exit 8
  cp rinp_r_sfc.f00 r_sfci || exit 8
  echo 'Restart files existed!!!!!!!'
  echo '==========rmtn.x and rinp.x are not executed!=========='
else
     if [ do$G2R = doyes ] ; then
       ln -fs $BASEDIR/sigf00 rb_sigf00
       ln -fs $BASEDIR/sfcf00 rb_sfcf00
     fi
     if [ do$P2R = doyes ] ; then
       ln -fs $BASEDIR/pgbf00 rb_pgbf00
     else
       if [ do$C2R = doyes ] ; then
         ln -fs $BASEDIR/r_sig.f00 rb_sigf00
         ln -fs $BASEDIR/r_sfc.f00 rb_sfcf00
       fi
     fi
     if [ do$NEWSST = do.TRUE. ] ; then
       ln -fs $BASEDIR/sstf00 rb_sstf00
     fi

#
#
#   Regional mountain
#
if [ do$RUNRMTN = doyes ] ; then
    $SHSDIR/rmtn.sh $MTNRES || exit 8
fi
#  Initial field for rsm run
#
if [ do$RUNRINP = doyes ] ; then
    $SHSDIR/rinp.sh $NEST 00 || exit 8
fi
if [ do$POSTTYPE = dosync ]; then
    $SHSDIR/rpgb_post.sh 00 || exit 8
fi
fi
#
#######################################
# Forecast loop
########################################
h=$FH
#FEND=6
while [ $h -lt $FEND ]; do
  hx=`expr $h + $INCHOUR`
  if [ $hx -gt $FEND ]; then  hx=$FEND; fi
  hh=$hx
  if [ $hx -lt 10 ];then hx=0$hx;fi
  hhr=`expr $h + 0`
  while [ $hhr -le $hx ]; do
       if [ $hhr -lt 10 ]; then hhr=0$hhr; fi
         rfti=`$UTLDIR/ndate $hhr $CDATE$CHOUR`
#        ln -fs $BASEDIR/sigf$rfti rb_sigf$hhr
#        ln -fs $BASEDIR/sfcf$rfti rb_sfcf$hhr
       if [ do$G2R = doyes ] ; then
         ln -fs $BASEDIR/sigf$hhr rb_sigf$hhr
         ln -fs $BASEDIR/sfcf$hhr rb_sfcf$hhr
       fi
       if [ do$P2R = doyes ] ; then
         ln -fs $BASEDIR/pgbf$hhr rb_pgbf$hhr
       else
         if [ do$C2R = doyes ] ; then
           ln -fs $BASEDIR/r_sig.f$hhr rb_sigf$hhr
#yj           ln -fs $BASEDIR/r_sfc.f$hhr rb_sfcf$hhr
         fi
       fi
       if [ do$NEWSST = do.TRUE. ] ; then
         ln -fs $BASEDIR/sstf$hhr rb_sstf$hhr
       fi
       hhr=`expr $hhr + $INCBASE`
  done
#
# rinp for g2c and l2c
if [ do$RUNRINP2 = doyes ] ; then
  $SHSDIR/rinp.sh $NEST ${hx}|| exit 8
fi
#
# fcst
#yj#if [ do$LAMMPI = doyes ]; then
#yj#   cp -f $EXPDIR/lamhosts . || exit      # always copy from submit directory
#yj#   recon lamhosts
#yj#   lamboot -v lamhosts
#yj#fi

if [ do$RUNFCST = doyes ] ; then
  echo "Forecast starting from hour $h..." >>stdout
  $SHSDIR/rfcstsfc.sh $hx || exit 8

#yj#if [ do$LAMMPI = doyes ]; then
#yj#   lamclean
#yj#   lamhalt
#yj#fi

#
#  Interpolate global surface forecast to regional grid to merge
#  $sfc_freq should be equal to or larger than $INCHOUR
#
# prepare for next fcst step
  cp r_sig.f$hx r_sigi || exit 8
  mv r_sigftdt r_sigitdt  || exit 8
  cp r_sfc.f$hx r_sfci || exit 8
#
# output at every PRTHOUR:
hr=`expr $h + $PRTHOUR`
while [ $hr -lt $hx ];do
  if [ $hr -lt 10 ];then hr=0$hr;fi
  mv r_sigf$hr   r_sig.f$hr
  mv r_sfcf$hr   r_sfc.f$hr
  mv r_flxf$hr   r_flx.f$hr
  hr=`expr $hr + $PRTHOUR`
done
#
#  r_pgb
#
if [ do$POSTTYPE = dosync ]; then

# no need at all, by zyf 2010.01.19
##if [ $FH = 00 ] ; then
##    $USHDIR/rpgb_post.sh $FH || exit 8
##fi

hr=$h
hr=`expr $h`
while [ $hr -lt $hx ];do
  if [ $hr -lt 10 ];then hr=0$hr;fi
  $SHSDIR/rpgb_post.sh $hr || exit 8
  hr=`expr $hr + $PRTHOUR`
done
#
if [ $hx -eq $ENDHOUR ]; then
   if [ $hx -lt 10 ];then hx=0$hx;fi
  mv r_sigf$hx   r_sig.f$hx
  mv r_sfcf$hx   r_sfc.f$hx
  mv r_flxf$hx   r_flx.f$hx
  $SHSDIR/rpgb_post.sh $hx || exit 8
fi

fi # r_pgb dosync

fi # fcst 
#submit avgmean 
   h=$hx
done
#

