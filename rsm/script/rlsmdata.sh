#! /bin/sh
#
set -ex
#
PROG=rlsmdata
WORKDIR=$1
dotype=$2
MACHINE=$3
WORK=$4
#
iocheck=0
#  namelist /naminfo/ dotype,dlon,dlat,pwest,psouth,c_min,c_max,projt,  &
#                     wordsize,tile_x,tile_y,tile_z,units
if [ $dotype = landuse ] ; then
  iocheck=`expr $iocheck + 1`
  pathn="modis_landuse_20class_15s"
  pathn2=" "
  echo " &naminfo" > typeinfo.parm
  echo " dotype='$dotype'," >> typeinfo.parm
  echo " dlon=0.00416667," >> typeinfo.parm
  echo " dlat=0.00416667," >> typeinfo.parm
  echo " pwest=-179.9979167," >> typeinfo.parm
  echo " psouth=-89.9979167," >> typeinfo.parm
  echo " c_min=1," >> typeinfo.parm
  echo " c_max=20," >> typeinfo.parm
  echo " projt=regular_ll," >> typeinfo.parm
  echo " wordsize=1," >> typeinfo.parm
  echo " tile_x=2400," >> typeinfo.parm
  echo " tile_y=2400," >> typeinfo.parm
  echo " tile_z=1," >> typeinfo.parm
  echo " units=category," >> typeinfo.parm
  echo " scale_fractor=1.0," >> typeinfo.parm
  echo " missing_value=-999.," >> typeinfo.parm
  echo " &END" >> typeinfo.parm
fi
if [ $dotype = soiltyp ] ; then
  iocheck=`expr $iocheck + 1`
  pathn="bnu_soiltype_top+NCU"
  pathn2=" "
  echo " &naminfo" > typeinfo.parm
  echo " dotype='$dotype'," >> typeinfo.parm
  echo " dlon=0.00833333," >> typeinfo.parm
  echo " dlat=0.00833333," >> typeinfo.parm
  echo " pwest=-179.99583," >> typeinfo.parm
  echo " psouth=-89.99583," >> typeinfo.parm
  echo " c_min=1," >> typeinfo.parm
  echo " c_max=16," >> typeinfo.parm
  echo " projt=regular_ll," >> typeinfo.parm
  echo " wordsize=1," >> typeinfo.parm
  echo " tile_x=1200," >> typeinfo.parm
  echo " tile_y=1200," >> typeinfo.parm
  echo " tile_z=1," >> typeinfo.parm
  echo " units=category," >> typeinfo.parm
  echo " scale_fractor=1.0," >> typeinfo.parm
  echo " missing_value=-999.," >> typeinfo.parm
  echo " &END" >> typeinfo.parm
  echo
fi
if [ $dotype = vegfrac ] ; then
  iocheck=`expr $iocheck + 1`
  pathn="veg_frac_lsa"
  pathn2=" "
  echo " &naminfo" > typeinfo.parm
  echo " dotype='$dotype'," >> typeinfo.parm
  echo " dlon=0.01," >> typeinfo.parm
  echo " dlat=0.01," >> typeinfo.parm
  echo " pwest=-179.995," >> typeinfo.parm
  echo " psouth=-89.995," >> typeinfo.parm
  echo " c_min=0.," >> typeinfo.parm
  echo " c_max=100.," >> typeinfo.parm
  echo " projt=regular_ll," >> typeinfo.parm
  echo " wordsize=1," >> typeinfo.parm
  echo " tile_x=1000," >> typeinfo.parm
  echo " tile_y=1000," >> typeinfo.parm
  echo " tile_z=12," >> typeinfo.parm
  echo " units=fraction," >> typeinfo.parm
  echo " scale_fractor=0.01," >> typeinfo.parm
  echo " missing_value=200.," >> typeinfo.parm
  echo " &END" >> typeinfo.parm
fi
     if [ $MACHINE = Fujitsu_fx10 ]||[ $MACHINE = Fujitsu_fx100 ]||[ $MACHINE = Fujitsu_fx1000 ] ; then
       ln -fs ${WORK}/rmtn_init/rmtnslm         fort.23
     else
       export FORT23=${WORK}/rmtn_init/rmtnslm
     fi

if [ ${RUNRSMTERRAIN} = yes ]; then
    ln -fsv "${TEMP}/rsm_terrain_init/rsm-lsm.bin" fort.23
fi


if [ $iocheck -ne 1 ] ; then
  echo "Error finding surface data type to process"
fi


#
ln -fs  ${FIXDIR}/GEOG_V391/* .
echo " &namcondir" >condir.parm
echo " condir='$WORKDIR'," >>condir.parm
echo " pathn ='$pathn'," >>condir.parm
echo " pathn2 ='$pathn2'," >>condir.parm
echo " &END" >>condir.parm
cat rsmlocationinit condir.parm typeinfo.parm > rlsmdata.parm
#
ls -l
pwd
##
ln -fs $EXPEXE/${PROG}.x $PROG.x
./$PROG.x <rlsmdata.parm >stdout.rlsmdata_$dotype  || exit
cat stdout.rlsmdata_$dotype
rm -f fort.[0-9]* 2>/dev/null

################  end of file  ##################

#>>----modis_landuse_20class_15s info----->>      
#type=categorical
#category_min=1  (c_min)
#category_max=20 (c_max)
#projection=regular_ll (projt)
#dx=0.00416667 (dlon)
#dy=0.00416667 (dlat)
#known_x=1.0
#known_y=1.0
#known_lat=-89.9979167  (psouth)
#known_lon=-179.9979167 (pwest)
#wordsize=1   (wordsize)
#tile_x=2400  
#tile_y=2400
#tile_z=1
#units="category"
#description="MODIS modified-IGBP landuse - 500 meter"
#mminlu="MODIFIED_IGBP_MODIS_NOAH"
#iswater=17
#isice=15
#isurban=13
#<<----modis_landuse_20class_15s info-----<<      
#>>----    bnu_soiltype_top+NCU      ----->>      
#type=categorical
#category_min=1
#category_max=16
#projection=regular_ll
#dx=0.00833333
#dy=0.00833333
#known_x=1.0
#known_y=1.0
#known_lat=-89.99583
#known_lon=-179.99583
#wordsize=1
#tile_x=1200
#tile_y=1200
#tile_z=1
#units="category"
#description="16-category top-layer soil type"
#>>----    bnu_soiltype_top+NCU      ----->>      
#<<----    greenfrac_fpar_modis      ----->>      
#type=continuous
#projection=regular_ll
#dx=0.00833333
#dy=0.00833333
#known_x=1.0
#known_y=1.0
#known_lat=-89.995833
#known_lon=-179.995833
#wordsize=1
#missing_value=200.
#tile_x=1200
#tile_y=1200
#tile_z=12
#scale_factor=0.01
#units="fraction"
#description="MODIS FPAR"
#>>----    greenfrac_fpar_modis      ----->>      

#>>---          veg_frac_lsa          ---->>
#type=continuous
#projection=regular_ll
#dx=0.01
#dy=0.01
#known_x=1.0
#known_y=1.0
#known_lat=-89.995
#known_lon=-179.995
#wordsize=1
#missing_value=200.
#tile_x=1000
#tile_y=1000
#tile_z=12
#scale_factor=0.01
#units="fraction"
#description="MODIS FPAR"
#<<---          veg_frac_lsa          ----<<

