# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

set(CMAKE_Fortran_COMPILER mpifrtpx)

set(NETCDF_PATH "/data/common/gfs/GEPSv3_lib/fx1000/netcdf-4.7.4")
set(HDF5_PATH "/data/common/gfs/GEPSv3_lib/fx1000/hdf5-1.10.7")
if(${USE_OMIP})
    set(ZLIB_PATH "/data/common/gfs/GEPSv3_lib/fx1000/zlib-1.2.11_fastest")
else()
    set(ZLIB_PATH "/data/common/gfs/GEPSv3_lib/fx1000/zlib-1.2.8")
endif()
set(PNETCDF_PATH "/data/common/gfs/GEPSv3_lib/fx1000/pnetcdf-1.12.3")
set(MCT_PATH "/data/common/gfs/GEPSv3_lib/fx1000/MCT")

set(CURRENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
#file(REAL_PATH "${CURRENT_DIR}/../" opath)
#set(CPL_PATH "${opath}/coupler/fx1000")
set(CPL_PATH "/data/common/gfs/GEPSv3_lib/coupler/fx1000")
set(CICE_PATH "/data/common/gfs/GEPSv3_lib/cice/cpl_lib_musoac")

message("CPL_PATH=${CPL_PATH}")

set(TIMCOM_INC_PATH
    ${CICE_PATH}/compile
    ${NETCDF_PATH}/include
    ${HDF5_INC_PATH}/include
    ${ZLIB_PATH}/include
    ${PNETCDF_PATH}/include
    ${MCT_PATH}/include
    ${CPL_PATH}
)
set(TIMCOM_LIB_PATH
    ${CICE_PATH}
    ${NETCDF_PATH}/lib
    ${HDF5_PATH}/lib
    ${ZLIB_PATH}/lib
    ${PNETCDF_PATH}/lib
    ${MCT_PATH}/lib
    ${CPL_PATH}
)
set(TIMCOM_OPTION
    -fw
    -Kfast,parallel,ocl,autoobjstack
    -CcdRR8
    -SSL2BLAMP
    -Ncheck_intrfunc
)
set(TIMCOM_LINK_LIB
    -lcpl
    -lpnetcdf
    -lnetcdff
    -lnetcdf
    -lmct
    -lmpeu
    -lhdf5_hl
    -lhdf5
    -lz
)

if(${USE_ICE})
    set(TIMCOM_LINK_LIB -lcice ${TIMCOM_LINK_LIB})
endif()

