#!/bin/ksh
set -x

export FLIB_CNTL_BARRIER_ERR=FALSE
DISK=$1
SDATE=$2
RESTRHR=$3                              # restart hour (zero means NOT restart)
FCSTHR=$4                               # forecast hour
RUNDIR=$5

# configure
. $DISK/run/configure_4cpl

##yj##---- user setting----
#diffusion on div/vor(DIFM) and T,P(DIFH)
DIFH=6     # org:6,when it's smaller, the diffusion is stronger. 
DIFM=2
DIFFT=10   # times of diffusion (if > 1, stronger diffusion)
DIFHSEC=`expr $DIFH \* $DELTAT_REG \/ $DTIMES `
DIFMSEC=`expr $DIFM \* $DELTAT_REG \/ $DTIMES \/ $DIFFT`
##yj##---- user setting----

####$SHSDIR/rsm_fcst.sh  
  mkdir -p $RUNDIR
  cd $RUNDIR || exit

if [ $CWB_MPMD = "no" ]&&[ $CWB_CFS = "no" ]; then
# link base field files to $RUNDIR  
# ln -fs  ${BASEDIR}/rsm_data* $RUNDIR
# ln -fs  ${BASEDIR}/rsm_idate* $RUNDIR
# rename "_${SDATE}00" "" rsm_data_*.f*
# rename "_${SDATE}" "" rsm_idate_*
# ln -fs  ${BASEDIR}/rsm_xlat_${SDATE}00.txt $RUNDIR/g2r_xlat.txt
# ln -fs  ${BASEDIR}/rsm_xlon_${SDATE}00.txt $RUNDIR/g2r_xlon.txt
  echo "link data code in 04_RSMTIM.ksh"
fi

#-----------------------------------------------
# determine model run parameters
#--------------n---------------------------------
if [ $RSFC_MERGE = "yes" ]; then
  ISFC_MERGE=1
else
  ISFC_MERGE=0
fi
FCSTSEC=`expr $FCSTHR \* 3600`
PRNTSEC=`expr $PRTHOUR \* 3600`
RSWRSEC=`expr $RSWRHOUR \* 3600`
RLWRSEC=`expr $RLWRHOUR \* 3600`
BASESEC=`expr $INCBASE \* 3600`
cplf_timcom=`expr $cplf_hr \* 3600`

if [ $NONHYD = "yes" ]; then
  FILTA=0.80
else
  FILTA=0.92
fi
IGRD1=`expr $IGRD + 1`
JGRD1=`expr $JGRD + 1`
RSFCSEC=`expr $sfc_freq \* 3600`;
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
  RESTRHR=$RESTRHR,
  READCLIMHR=$READCLIMHR,
  ISEMIMP=$ISEMIMP,
  IIMPRLX=$IIMPRLX,
  IDMPJET=$IDMPJET,
  IMDLPHY=$IMDLPHY,
  IOUTNHR=$IOUTNHR,
  ISFCMRG=$ISFC_MERGE,
  IRCLIM=$IRCLIM,
  cplf_timcom=$cplf_timcom,

  &END
EOF
#
  cat rsmlocation >> rfcstparm

  echo "---in run.sh of rsm---"
  cp $STTPRM $RUNDIR/station.parm

########## RUNRINP ###########
#if [ do$P2R = doyes ] ; then
   sig2rg=.FALSE.
   sfc2rg=.FALSE.
   newmtn=.TRUE.
   pgb2rg=.TRUE.
#fi  
#
    echo " Regional input starts:"
    echo " &NAMRIN                                                "  >rinpparm
    echo "    SIG2RG=$sig2rg,SFC2RG=$sfc2rg,PERCMTN=$PERCMTN,     " >>rinpparm
    echo "    NEWSIG=$NEWSIG,NEWMTN=$newmtn,       " >>rinpparm
    echo "    PGB2RG=$pgb2rg,NEWSST=$NEWSST,                      " >>rinpparm
    echo " &END                                                   " >>rinpparm
#
#### define rinp parameters ######
cat >rsmlocationbase <<EOF
 &NAMLOCBASE
  RPROJ    = $BRPROJ,
  RTRUTH   = $BRTRUTH,
  RORIENT  = $BRORIENT,
  RDELX    = $BRDELX,
  RDELY    = $BRDELY,
  RCENLAT  = $BRCENLAT,
  RCENLON  = $BRCENLON,
  RLFTGRD  = $BRLFTGRD,
  RBTMGRD  = $BRBTMGRD,
  CPROJ    = $BCPROJ,
  CTRUTH   = 0,
  CORIENT  = 0,
  CDELX    = 0,
  CDELY    = 0,
  CCENLAT  = 0,
  CCENLON  = 0,
  CLFTGRD  = 0,
  CBTMGRD  = 0,
  CLAT2    = $BCLAT2,
  CLAT1    = $BCLAT1,
  CLON2    = $BCLON2,
  CLON1    = $BCLON1,
 &END

EOF
    cat rsmlocationbase >> rinpparm
#    
cat >rsmlocationinit <<EOF
 &NAMLOCINIT
  RPROJ    = $RPROJ,
  RTRUTH   = $RTRUTH,
  RORIENT  = $RORIENT,
  RDELX    = $RDELX,
  RDELY    = $RDELY,
  RCENLAT  = $RCENLAT,
  RCENLON  = $RCENLON,
  RLFTGRD  = $RLFTGRD,
  RBTMGRD  = $RBTMGRD,
  CPROJ    = $BCPROJ,
  CTRUTH   = 0,
  CORIENT  = 0,
  CDELX    = 0,
  CDELY    = 0,
  CCENLAT  = 0,
  CCENLON  = 0,
  CLFTGRD  = 0,
  CBTMGRD  = 0,
  CLAT2    = $BCLAT2,
  CLAT1    = $BCLAT1,
  CLON2    = $BCLON2,
  CLON1    = $BCLON1,
 &END

EOF
    cat rsmlocationinit >> rinpparm

    echo " &namsig                                         " >>rinpparm
    echo " &end                                            " >>rinpparm
    echo " &namsfc                                         " >>rinpparm
    echo " &end                                            " >>rinpparm
    echo " &namclim                                        " >>rinpparm
if [ $CLIM = 1 ]
then
    echo "    iclim=1,                                    " >>rinpparm
    echo "     fnmskh=\"$FNMSKH\",                             " >>rinpparm
    echo "     fnalbc=\"$FNALBC\",                             " >>rinpparm
    echo "     fnsotc=\"$FNSOTC\",                             " >>rinpparm
    echo "     fnvegc=\"$FNVEGC\",                             " >>rinpparm
    echo "     fnvetc=\"$FNVETC\",                             " >>rinpparm
    echo "     fnzorc=\"$FNZORC\",                             " >>rinpparm
    echo "     fntg3c=\"$FNTG3C\",                             " >>rinpparm
    echo "     fnslpc=\"$FNSLPC\",                             " >>rinpparm
    echo "     fnabsc=\"$FNABSC\",                             " >>rinpparm
    echo "     fnvmnc=\"$FNVMNC\",                             " >>rinpparm
    echo "     fnvmxc=\"$FNVMXC\",                             " >>rinpparm
fi
    echo " &end                                            " >>rinpparm

#####
#
##########################################################
#
# RSM INITIAL forecast
#
   if [ do$P2R = doyes ] ; then
   if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
      ln -fs ${FIXDIR}/terr_grib2.bin fort.16 
      ln -fs ${FIXDIR}/slmk_grib2.bin fort.17
      ln -fs ${FIXDIR}/terr_grib2.bin fort.26 
      ln -fs ${FIXDIR}/slmk_grib2.bin fort.27
   else
     export FORT16=${FIXDIR}/terr_grib2.bin
     export FORT17=${FIXDIR}/slmk_grib2.bin
     export FORT26=${FIXDIR}/terr_grib2.bin
     export FORT27=${FIXDIR}/slmk_grib2.bin
   fi
   fi
####
#   Regional mountain
########## RUNRMTN ###########
if [ do$RUNRMTN = doyes ] ; then
    if [ -s ${TEMP}/rmtn_init/rmtnvar ]&&[ -s ${TEMP}/rmtn_base/rmtnvar ];then
      echo "RMTN data exist, skip!"
    else
      mtntype=init
      mkdir -p ${TEMP}/rmtn_${mtntype}
      cd ${TEMP}/rmtn_${mtntype} || exit 8
      cp $RUNDIR/rsmlocation${mtntype} .
      sed -i '1c &NAMLOC' rsmlocation${mtntype}
      ${SHSDIR}/rmtn.sh $MTNRES $mtntype ${TEMP}/rmtn_${mtntype}|| exit 8

      mtntype=base
      if [ $SRTM1ARC = yes ]; then
        # use GTOPO30 at base field
        export MTNRES=30
      fi
      mkdir -p ${TEMP}/rmtn_${mtntype}
      cd ${TEMP}/rmtn_${mtntype} || exit 8
      cp $RUNDIR/rsmlocation${mtntype} .
      sed -i '1c &NAMLOC' rsmlocation${mtntype}
      ${SHSDIR}/rmtn.sh $MTNRES $mtntype ${TEMP}/rmtn_${mtntype}|| exit 8
    fi
    cd $RUNDIR || exit
fi
###
echo "after rmtn"

if [ ${RUNRSMTERRAIN} = yes ]; then
    cd "${TEMP}" || exit

    if [ -f "rsm_terrain_init/rsm-oro.bin" ] && [ -f "rsm_terrain_base/rsm-oro.bin" ]; then
        echo "rsm terrain (mtnlm7_slm30g) will not run because output files already exist"
    else
        OUTPUTTYPE=init

        echo "Running rsm terrain (mtnlm7_slm30g) for ${OUTPUTTYPE}"
        mkdir -pv "rsm_terrain_${OUTPUTTYPE}" && cd "rsm_terrain_${OUTPUTTYPE}" || exit
        ln -fsv "${FIXDIR}/thirty.second.antarctic.new.bin" fort.15  || exit
        if [ ${GMTED2010} = yes ]; then
            ln -fsv "${FIXDIR}/gmted2010.30sec.fine"        fort.235 || exit
        elif [ ${GTOPO30} = yes ]; then
            ln -fsv "${FIXDIR}/gtopo30_gg.fine"             fort.235 || exit
        else
            echo "When RUNRSMTERRAIN is yes, one of GMTED2010 and GTOPO30 must also be yes."
            exit
        fi
        cat "${EXPDIR}/rsm-domain.nml" >rsm-domain.nml || exit
        ln -fsv "${EXPEXE}/mtnlm7_slm30g.exe" . || exit
        echo "${IGRD} ${JGRD}" | ./mtnlm7_slm30g.exe -Wl,-T15,-T235
        cd ..

        mkdir -pv "rmtn_${OUTPUTTYPE}" && cd "rmtn_${OUTPUTTYPE}" || exit
        cp -fv "${RUNDIR}/rsmlocation${OUTPUTTYPE}" .
        sed -i '1c &NAMLOC' "rsmlocation${OUTPUTTYPE}"
        ln -fsv "../rsm_terrain_${OUTPUTTYPE}/rsm-oro.bin" .
        ${SHSDIR}/rmtn.sh ${MTNRES} ${OUTPUTTYPE} ${TEMP}/rmtn_${OUTPUTTYPE} || exit
        cp -fv rmtnoss "../rsm_terrain_${OUTPUTTYPE}/rsm-oro-smoothed.bin"
        cd ..

        OUTPUTTYPE=base

        echo "Running rsm terrain (mtnlm7_slm30g) for ${OUTPUTTYPE}"
        mkdir -pv "rsm_terrain_${OUTPUTTYPE}" && cd "rsm_terrain_${OUTPUTTYPE}" || exit
        ln -fsv "${FIXDIR}/thirty.second.antarctic.new.bin" fort.15  || exit
        ln -fsv "${FIXDIR}/gmted2010.30sec.fine"            fort.235 || exit
        cat >rsm-domain.nml <<RSMDOMAIN
&rsmdomain
    rproj   = ${BRPROJ},
    rtruth  = ${BRTRUTH},
    rorient = ${BRORIENT},
    rdelx   = ${BRDELX},
    rdely   = ${BRDELY},
    rcenlat = ${BRCENLAT},
    rcenlon = ${BRCENLON},
    rlftgrd = ${BRLFTGRD},
    rbtmgrd = ${BRBTMGRD}
/
RSMDOMAIN
        ln -fsv "${EXPEXE}/mtnlm7_slm30g.exe" . || exit
        echo "${BIGRD} ${BJGRD}" | ./mtnlm7_slm30g.exe -Wl,-T15,-T235
        cd ..

        mkdir -pv "rmtn_${OUTPUTTYPE}" && cd "rmtn_${OUTPUTTYPE}" || exit
        cp -fv "${RUNDIR}/rsmlocation${OUTPUTTYPE}" .
        sed -i '1c &NAMLOC' "rsmlocation${OUTPUTTYPE}"
        ln -fsv "../rsm_terrain_${OUTPUTTYPE}/rsm-oro.bin" .
        ${SHSDIR}/rmtn.sh ${MTNRES} ${OUTPUTTYPE} ${TEMP}/rmtn_${OUTPUTTYPE} || exit
        cp -fv rmtnoss "../rsm_terrain_${OUTPUTTYPE}/rsm-oro-smoothed.bin"
        cd ..
    fi

    cd "${RUNDIR}" || exit
fi

########## RUNLSMDATA ###########
#   Regional LSM data
########## RUNLSMDATA ###########
if [ do$RUNLSMDATA = doyes ] ; then
    if [ -s ${TEMP}/rlsmdata/out_landuse.bin ];then
      echo "LSM data exist, skip!"
    else
      mkdir -p ${TEMP}/rlsmdata
      cd ${TEMP}/rlsmdata || exit 8
      cp $RUNDIR/rsmlocationinit .
      sed -i '1c &NAMLOC' rsmlocationinit
      dotype=landuse         # landuse , soiltyp , vegfrac
      ${SHSDIR}/rlsmdata.sh ${TEMP}/rlsmdata $dotype $MACHINE $TEMP || exit 8
      dotype=soiltyp         # landuse , soiltyp , vegfrac
      ${SHSDIR}/rlsmdata.sh ${TEMP}/rlsmdata $dotype $MACHINE $TEMP || exit 8
      dotype=vegfrac         # landuse , soiltyp , vegfrac
      ${SHSDIR}/rlsmdata.sh ${TEMP}/rlsmdata $dotype $MACHINE $TEMP || exit 8
      ln -fs  ${TEMP}/rlsmdata/out_landuse.bin ${FIXDIR}/.
      ln -fs  ${TEMP}/rlsmdata/out_soiltype.bin ${FIXDIR}/.
      ln -fs  ${TEMP}/rlsmdata/out_vegfrac.bin ${FIXDIR}/.
    fi
  cd $RUNDIR || exit
fi
#####-end Regional LSM data
#
########## RUN3DRTDATA ###########
#   Regional 3DRT data
########## RUN3DRTDATA ###########
if [ do$RUN3DRTDATA = doyes ] ; then
    if [ -s ${TEMP}/r_3drtdata/rsm_3drt.bin ];then
      echo "3DRT data exist, skip!"
    else
      mkdir -p ${TEMP}/r_3drtdata
      cd ${TEMP}/r_3drtdata || exit 8
      cp $RUNDIR/rsmlocationinit .
      sed -i '1c &NAMLOC' rsmlocationinit
      ${SHSDIR}/r_3drtdata.sh ${TEMP}/r_3drtdata $MACHINE $TEMP || exit 8
      ln -fs  ${TEMP}/r_3drtdata/rsm_3drt.bin ${FIXDIR}/rsm_3drt.bin
    fi
  cd $RUNDIR || exit
fi
if [ do$do3DRT = doyes ] ; then
    ln -fs  ${FIXDIR}/rsm_3drt.bin .
fi
#####-end 3DRT data
   if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
     ln -fsv "${TEMP}/rsm_terrain_base/rsm-lsm.bin" fort.13
     ln -fsv "${TEMP}/rsm_terrain_base/rsm-oro-smoothed.bin" fort.14
     ln -fs ${FIXDIR}/siglevel.l${LEVR}.txt       fort.15
     ln -fsv "${TEMP}/rsm_terrain_init/rsm-lsm.bin" fort.23
     ln -fsv "${TEMP}/rsm_terrain_init/rsm-oro-smoothed.bin" fort.24
     ln -fs ${FIXDIR}/siglevel.l${LEVR}.txt       fort.25
     if [ ${RUNRMTN} = yes ]; then
       ln -fs ${TEMP}/rmtn_base/rmtnslm         fort.13
       ln -fs ${TEMP}/rmtn_base/rmtnoss         fort.14
       ln -fs ${TEMP}/rmtn_init/rmtnslm         fort.23
       ln -fs ${TEMP}/rmtn_init/rmtnoss         fort.24
     fi
   else
     export FORT13="${TEMP}/rsm_terrain_base/rsm-lsm.bin"
     export FORT14="${TEMP}/rsm_terrain_base/rsm-oro-smoothed.bin"
     export FORT15=${FIXDIR}/siglevel.l${LEVR}.txt
     export FORT23="${TEMP}/rsm_terrain_init/rsm-lsm.bin"
     export FORT24="${TEMP}/rsm_terrain_init/rsm-oro-smoothed.bin"
     export FORT25=${FIXDIR}/siglevel.l${LEVR}.txt
     if [ ${RUNRMTN} = yes ]; then
       export FORT13=${TEMP}/rmtn_base/rmtnslm
       export FORT14=${TEMP}/rmtn_base/rmtnoss
       export FORT23=${TEMP}/rmtn_init/rmtnslm
       export FORT24=${TEMP}/rmtn_init/rmtnoss
     fi
   fi

  if [ $CLIM = 1 ] ;
  then
    ln -fs ${FIXDIR}/$FNMSKH .
    ln -fs ${FIXDIR}/$FNALBC .
    ln -fs ${FIXDIR}/$FNSOTC .
    ln -fs ${FIXDIR}/$FNVEGC .
    ln -fs ${FIXDIR}/$FNVETC .
  if [ do$NOAHLSM = dono  ]; then
    ln -fs ${FIXDIR}/$FNZORC .
  fi
    ln -fs ${FIXDIR}/$FNTG3C .
    ln -fs ${FIXDIR}/$FNSLPC .
    ln -fs ${FIXDIR}/$FNABSC .
    ln -fs ${FIXDIR}/$FNVMNC .
    ln -fs ${FIXDIR}/$FNVMXC .
  fi

# wind_fitch
  if [ do$wind_fitch = doyes ] ; then
    ln -fs ${FIXDIR}/wind_fitch/wind* .
  fi

#
#######################################
# Forecast loop
########################################
#yj#h = 0
#yj#hx = h + $INCHOUR
#yjif [ do$RUNFCST = doyes ] ; then
#yj#  $SHSDIR/rfcstsfc.sh $hx || exit 8
###yj###   ${WORK}/rfcstsfc.sh || exit 8
 if [ $MACHINE = Fujitsu_fx1000 ] ; then
 TUNE1=$FIXDIR/global_cldtune.f77_lt
 O3CLIM=$FIXDIR/global_o3clim.txt
 O3PROD=$FIXDIR/global_o3prod.f77_lt
 O3LOSS=$FIXDIR/global_o3loss.f77_lt
 AERDIR=$FIXDIR
 else
 TUNE1=$FIXDIR/global_cldtune.f77
 O3CLIM=$FIXDIR/global_o3clim.txt
 O3PROD=$FIXDIR/global_o3prod.f77
 O3LOSS=$FIXDIR/global_o3loss.f77
 AERDIR=$FIXDIR
 fi
   if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
     ln -fs $TUNE1     fort.43
     ln -fs $O3CLIM    fort.45  #yj -ori=48
     ln -fsv "${TEMP}/rsm_terrain_init/rsm-hprime.bin" fort.42
     ln -fs $O3PROD    fort.35  #yj -ori=28
     ln -fs $O3LOSS    fort.36  #yj -ori=29
     ln -fs r_sigf     fort.30
     ln -fs r_sigftdt  fort.31
     ln -fs r_sfcf     fort.32
     if [ ${RUNRMTN} = yes ]; then
       ln -fs ${TEMP}/rmtn_init/rmtnvar    fort.42
     fi
   else
     export FORT43=$TUNE1
     export FORT45=$O3CLIM
     export FORT42="${TEMP}/rsm_terrain_init/rsm-hprime.bin"
     export FORT35=$O3PROD
     export FORT36=$O3LOSS
     export FORT30=r_sigf
     export FORT31=r_sigftdt
     export FORT32=r_sfcf
     if [ ${RUNRMTN} = yes ]; then
     export FORT42=${TEMP}/rmtn_init/rmtnvar
     fi
   fi

for m in 01 02 03 04 05 06 07 08 09 10 11 12
do
  ln -fs $AERDIR/global_aeropac3a.m$m.txt aeropac3a.m$m
done
   if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
    ln -fs r_sigbase      fort.34
    ln -fs r_sigf         fort.80
    ln -fs r_sigftdt      fort.81
    ln -fs r_sfcf         fort.82
    ln -fs r_flx.f00      fort.74
    ln -fs r_init         fort.78
   else
    export FORT34=r_sigbase
    export FORT80=r_sigf
    export FORT81=r_sigftdt
    export FORT82=r_sfcf
    export FORT74=r_flx.f00
    export FORT78=r_init
   fi
cat ${RUNDIR}/rinpparm ${RUNDIR}/rfcstparm > rfcstparm.all
if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
  ln -fs rfcstparm.all fort.5
else
  export FORT5=rfcstparm.all
fi

echo ${WORK} ${RUNDIR}
cp -p ${WORK}/configure_4cpl ${RUNDIR}/exp_configure_4cpl
cp -p ${WORK}/run_4cpl.sh ${RUNDIR}/exp_run.sh

# post process
###################################
# Defining namelist
###################################

IGRD1=`expr $IGRD + 1`
JGRD1=`expr $JGRD + 1`
INTERP=${INTERP:-no}

#...........Grid pacific...........................................
if [ $INTERP = yes ] ; then
echo " &NAMPGB                                          "   > rpgb.parm
echo "         NCPUS=4, KO=$KO,                         "  >> rpgb.parm
echo "         IO=$IO, JO=$JO,                          "  >> rpgb.parm
echo "         PROJO=$PROJO,                            "  >> rpgb.parm
echo "         RLAT1O=$RLAT1O, RLON1O=$RLON1O,          "  >> rpgb.parm
echo "         RLAT2O=$RLAT2O, RLON2O=$RLON2O,          "  >> rpgb.parm
echo "         NTRAC=1, NCLD=$NCLD,NEWSLM=$NEWSLM,      "  >> rpgb.parm
echo " &END                                             "  >> rpgb.parm
else
echo " &NAMPGB                                 "   > rpgb.parm
echo "         NCPUS=4, KO=$KO,                 "  >> rpgb.parm
echo "         IO=$IGRD1, JO=$JGRD1,            "  >> rpgb.parm
echo "         NTRAC=1, NCLD=$NCLD,             "  >> rpgb.parm
echo "         PO= 1000.,975.,950.,925.,900.,   " >> rpgb.parm
echo "        850.,800.,750.,700.,600.,         " >> rpgb.parm
echo "        500.,400.,300.,250.,200.,         " >> rpgb.parm
echo "        150.,100., 70., 50., 30.,         " >> rpgb.parm
echo "        20., 10.,  0.,                    " >> rpgb.parm
echo " &END                                 "     >> rpgb.parm
fi  

if [ $? -ne 0 ] ; then
   exit 8
fi
