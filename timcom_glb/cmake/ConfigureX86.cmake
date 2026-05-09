# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

set(CMAKE_Fortran_COMPILER mpif90)
set(GEPS_LIBS $ENV{GEPS_LIBS})

set(NETCDF_PATH "${GEPS_LIBS}/x86_64/nvidia/netcdf-4.9.0")
set(PNETCDF_PATH "${GEPS_LIBS}/x86_64/nvidia/pnetcdf-1.12.3")
set(MCT_PATH "${GEPS_LIBS}/x86_64/nvidia/mct-2.11.0")
#file(REAL_PATH "${CURRENT_DIR}/../" opath)
#set(CPL_PATH "${opath}/coupler/a100")
set(CPL_PATH "${GEPS_LIBS}/x86_64/nvidia/coupler")
#GPU port can't use CICE for now
#set(CICE_PATH "${GEPS_LIBS}/x86_64/nvidia/cpl_lib_musoac_a100")

message("CPL_PATH=${CPL_PATH}")


set(TIMCOM_INC_PATH
    ${CICE_PATH}/compile
    ${NETCDF_PATH}/include
    ${PNETCDF_PATH}/include
    ${MCT_PATH}/include
    ${CPL_PATH}
)
set(TIMCOM_LIB_PATH
    ${CICE_PATH}
    ${NETCDF_PATH}/lib
    ${PNETCDF_PATH}/lib
    ${MCT_PATH}/lib
    ${CPL_PATH}
)
set(TIMCOM_OPTION
    -Ofast
    -Mr8
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
    -lblas
)

if(${USE_NVTX})
    list(APPEND EX_FUNC
        density_t_s_tt_ts
        potential_temperature_t
        temperature
        p_inprod
        mpi_exch_r8type2d
        p_sqprod
        p_sipsqelim
        nccl_check_helper
        ncclGroupStart
        ncclGroupEnd
        ncclSend
        ncclRecv
        ncclAllReduce
        make_local_coefficients
        eos_mkcoef
        ocn_cuda_check_helper
    )
    string(REPLACE ";" "," EX_FUNC "${EX_FUNC}")
    set(TIMCOM_OPTION ${TIMCOM_OPTION} -Minstrument -Minstrument-exclude-func-list=${EX_FUNC} -traceback -lnvhpcwrapnvtx)
endif()

if(${USE_GPU})
    set(TIMCOM_OPTION ${TIMCOM_OPTION} -DUSE_GPU -acc=gpu -gpu=cc${GPU_ARCH},cuda${CUDA_RUNTIME_VERSION} -cuda -cudalib=nccl -Minfo=accel)
endif()

if(${USE_ICE})
    set(TIMCOM_LINK_LIB -lcice ${TIMCOM_LINK_LIB})
endif()
