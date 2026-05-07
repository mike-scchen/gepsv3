#!/bin/ksh
# build library and some executable files

#
# at_MACHINE : fx1000, pcc
#
if [ $# == 2 ] ; then
  export at_MACHINE=${1}
  export do_case=${2}
else
  echo "usage: $0 [at_MACHINE] [do_case]"
  echo "  at_MACHINE can be fx1000, or pcc"
  echo "  do_case can be test, 5km, 12km, user, no"
  echo "  !!!Notice!!! When do_case equals user, you need to setup your own case and compile!"
  echo "  FOR EXAMPLE: ./build.sh fx1000 test"
  echo "               |-> This will build a test case for you"
  echo "  !!!Notice!!! When do_case equals no, only the libraries, utilies and fix data are linked"
  echo "  FOR EXAMPLE: ./build.sh fx1000 no"
  echo "               |-> No experiment is setup" 
  exit
fi

export PWDPATH=` pwd `

# compile options for pcc:linux_intel, fx1000:Fujitsu_fx1000
if [ ${at_MACHINE} = pcc ] ; then
  export MACHINE=linux_intel
elif [ ${at_MACHINE} = fx1000 ] ; then
  export MACHINE=Fujitsu_fx1000
else
  echo "can't support your machine"
  exit 88
fi

#---------------------lib--------------------------
export LIBPATH=`cd ./lib ; pwd `

if [ ${MACHINE} = linux_intel ] ; then
export LINKPATH=/home/xb127/pkg
elif [ ${MACHINE} = Fujitsu_fx1000 ] ; then
export LINKPATH=/nwpr/gfs/xb127/pkg/fx1000
else
echo "can't support your machine"
echo "please define your path for NCEPLIB and netcdf library"
echo " export LINKPATH=..."
exit 88
fi

if [ ! -d "${LINKPATH}" ];then echo "package path doesn't exist, please build the pkg";exit 88;fi

#clean all old lib and unlink every file in 
rm -rf ${LIBPATH}/*.a ${LIBPATH}/*.la ${LIBPATH}/incmod
find -maxdepth 1 -type l -delete

#build lib comes with RSM
#cd ${LIBPATH}/src || exit 81
#. compile.sh  || exit 82

#link library to LIBPATH
cd ${LIBPATH} || exit 81
[ -d "${LINKPATH}/jasper-1.900.1" ] && ln -fs ${LINKPATH}/jasper-1.900.1 . || echo "${LINKPATH}/jasper-1.900.1 not exist"
if [ ${MACHINE} = linux_intel ] ; then
  [ -d "${LINKPATH}/libpng-1.2.56" ] && ln -fs ${LINKPATH}/libpng-1.2.56 . || echo "${LINKPATH}/libpng-1.2.56 not exist"
elif [ ${MACHINE} = Fujitsu_fx1000 ] ; then
  [ -d "${LINKPATH}/libpng-1.6.37" ] && ln -fs ${LINKPATH}/libpng-1.6.37 . || echo "${LINKPATH}/libpng-1.6.37 not exist"
else
  [ -d "${LINKPATH}/libpng-1.2.50" ] && ln -fs ${LINKPATH}/libpng-1.2.50 . || echo "${LINKPATH}/libpng-1.2.50 not exist"
fi

if [ ${MACHINE} = Fujitsu_fx1000 ] ; then
[ -d "${LINKPATH}/zlib-1.2.3" ] && ln -fs ${LINKPATH}/zlib-1.2.3 . || echo "${LINKPATH}/zlib-1.2.3 not exist"
else
[ -d "${LINKPATH}/zlib-1.2.8" ] && ln -fs ${LINKPATH}/zlib-1.2.8 . || echo "${LINKPATH}/zlib-1.2.8 not exist"
fi


#link new ncep library to LIBPATH
cd ${LIBPATH} || exit 81
[ -d "${LINKPATH}/NCEPLIB" ] && ln -fs ${LINKPATH}/NCEPLIB . || echo "${LINKPATH}/NCEPLIB not exist"

#link netcdf library to LIBPATH
cd ${LIBPATH} || exit 81
[ -d "${LINKPATH}/netcdf-3.6.3" ] && ln -fs ${LINKPATH}/netcdf-3.6.3 . || echo "${LINKPATH}/netcdf-3.6.3 not exist"


#
#-------------------cwb mpmd-----------------------
#

cd ${PWDPATH}/src/cwb_mpmd.fd;
mkdir -p $LIBPATH/cwb_mpmd;
if [ ${at_MACHINE} = pcc ] ; then
  make -f Makefile_pcc
elif [ ${at_MACHINE} = fx1000 ] ; then
  make -f Makefile_fx1000
fi
#
#-------------------post processor-----------------
#
PGRB_VERSION=2
if [[ ${PGRB_VERSION} = 2 ]];then
  cd ${PWDPATH}/src/rsm_pgrb2.fd
else
  cd ${PWDPATH}/src/rsm_pgrb.fd
fi
if [ ${at_MACHINE} = pcc ] ; then
  make -f makefile_awips_ifort
elif [ ${at_MACHINE} = fx1000 ] ; then
  make -f makefile_awips_fx1000
fi

#
#--------------required data base------------------
#
# landuse, soiltype, vegetation fraction from WRF data base
mkdir -p ${PWDPATH}/fix/GEOG_V391
if [ ${at_MACHINE} = pcc ] ; then
  ln -fs /home/xb127/DATA/DATA_GEOG_V391/* ${PWDPATH}/fix/GEOG_V391/
elif [ ${at_MACHINE} = fx1000 ] ; then
#  ln -fs /nwpr/gfs/xb127/DATA/GEOG_V391/* ${PWDPATH}/fix/GEOG_V391/
  ln -sf /data/common/gfs/GEPSv3_lib/data/GEOG_V391/* ${PWDPATH}/fix/GEOG_V391/
  ln -fs /data/common/gfs/GEPSv3_lib/data/veg_frac_lsa ${PWDPATH}/fix/GEOG_V391/
else
  ln -fs /nwpr/gfs/xb127/data2/DATA/GEOG_V391/* ${PWDPATH}/fix/GEOG_V391/
fi
#
# 3drt data from Dr. Lee Wei-Liang
mkdir -p ${PWDPATH}/fix/Param_3DRT
if [ ${at_MACHINE} = pcc ] ; then
  ln -fs /home/xb127/DATA/PARA_3DRT/Param_3DRT/* ${PWDPATH}/fix/Param_3DRT
elif [ ${at_MACHINE} = fx1000 ] ; then
#  ln -fs /nwpr/gfs/xb127/DATA/PARA_3DRT/Param_3DRT/* ${PWDPATH}/fix/Param_3DRT
  ln -fs /data/common/gfs/GEPSv3_lib/data/PARA_3DRT/Param_3DRT/* ${PWDPATH}/fix/Param_3DRT
else
  ln -fs /nwpr/gfs/xb127/data2/DATA/PARA_3DRT/Param_3DRT/* ${PWDPATH}/fix/Param_3DRT
fi

#
#--------------base field data from GFS------------
#
mkdir -p ${PWDPATH}/wrk/DATA_GFS;
if [ ${at_MACHINE} = pcc ] ; then
  echo "============= Error ============"
  echo "GFS data not prepared for pcc, please run GFS to output data for RSM, or copy from fx1000!!!"
  echo "GFS data not prepared for pcc, please run GFS to output data for RSM, or copy from fx1000!!!"
  echo "GFS data not prepared for pcc, please run GFS to output data for RSM, or copy from fx1000!!!"
  echo "--- after you get the GFS and , please ues this build.sh to build RSM again.---"
  exit 88
elif [ ${at_MACHINE} = fx1000 ] ; then
  echo "GFS data at /nwpr/gfs/xb127/DATA/DATA_GFS/synOP_ZC_21060200_initT511"
#  ln -fs /nwpr/gfs/xb127/DATA/DATA_GFS/gfs_21100100_xnew ${PWDPATH}/wrk/DATA_GFS/
  ln -fs /data/common/gfs/GEPSv3_lib/data/DATA_GFS/gfs_21100100_xnew ${PWDPATH}/wrk/DATA_GFS/
else
  echo "GFS data at /nwpr/gfs/xb127/data2/DATA/DATA_GFS/synOP_ZC_21060200_initT511"
  ln -fs /nwpr/gfs/xb127/DATA/DATA_GFS/gfs_21100100_xnew ${PWDPATH}/wrk/DATA_GFS/
fi

#
#----------convert RSM fix file to little endian---
#
if [ ${at_MACHINE} = fx1000 ] ; then
  cd ${PWDPATH}/fix/tran_b2l.dir || exit 88
  ksh "do.sh" || exit 99
fi


#
#-------------------model components---------------
#
cat > ${PWDPATH}/run/set_machine << EOF
export MACHINE=${MACHINE}
export DEBUG=no
EOF

if [ ${do_case} = no ]; then
  echo " "
  echo "only libraries and some fix data are built and linkded"
  echo "to setup experiments, please specify [do case] when you run build.sh"
  echo "or create your own run/exp/configure_exp and link it to run/configure"
  echo " "
else
if [ ${do_case} = user ]; then
  cd ${PWDPATH}/run || exit 88
  export new_userfile=configure_${do_case}
  (set -C &&  cat < ${PWDPATH}/run/exp/configure_test > ${PWDPATH}/run/exp/${new_userfile})
  if [ $? -ne 0 ] ; then
    export new_userfile=configure_${do_case}_`date +%Y%m%d%H%S`
    echo " "
    echo "---There is already a file called configure_user in run/exp/---"
    echo "---to avoid overwritten, a new user configure file is named as ${new_userfile}---"
    (set -C &&  cat < ${PWDPATH}/run/exp/configure_test > ${PWDPATH}/run/exp/${new_userfile})
  fi
  ln -fs ${PWDPATH}/run/exp/${new_userfile} ${PWDPATH}/run/configure;
  echo " "
  echo "No default case is set. Go to run/configure to setup your own case"
  echo "After setup configure, please run ./compile "
  echo "The platform/machine is set to be ${at_MACHINE}"
  echo " "
else
  if [ ${do_case} = test ]||[ ${do_case} = 5km ]||[ ${do_case} = 12km ]||[ ${do_case} = SUMcV2.0 ]||[ ${do_case} = SUMcV2.0_vote ]; then
    cd ${PWDPATH}/run || exit 88
    ln -fs "${PWDPATH}/run/exp/configure_${do_case}"      "${PWDPATH}/run/configure"
    ln -fs "${PWDPATH}/run/exp/rsm-domain-${do_case}.nml" "${PWDPATH}/run/rsm-domain.nml"
    ./compile_4cpl > compile_4cpl.log 2>&1 &
    echo "The model - RSM is building ... at ${at_MACHINE}"
    echo "The model - RSM is building ... at ${at_MACHINE}"
    echo "The model - RSM is building ... at ${at_MACHINE}"
    echo "check run/compile.log to see if compilation is finished or not"
    echo " "

    echo "         / ====|                                            "
    echo "        /      |                                            "
    echo "       /       |                                            "
    echo "      /        |                                            "
    echo "     /         |                                            "
    echo "    /          |                                            "
    echo "   |           |                                            "
    echo "   |  hammer   |                                            "
    echo "   |           |                                            "
    echo "   |   as      |==================================          "
    echo "   |           |  Check run/compile_4cpl.log !!!     ||          "
    echo "   |   a       |  or RSM may crush like smashed ||          "
    echo "   \           /==================================          "
    echo "    | reminder|                                             "
    echo "   /           \                                            "
    echo "   |           |                                            "
    echo "   =============                                            "
    echo "     #$%^!!%$^*                                             "

  else
    echo " "
    echo "===Fail to build==="
    echo "invalid [do_case] is found"
    echo "[do_case] should be test, 5km, 12km, or user"
    echo "===Fail to build==="
    echo " "
  fi
fi
fi

