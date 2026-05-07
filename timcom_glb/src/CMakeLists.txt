# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

file(GLOB TIMCOM_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.F90")
list(REMOVE_ITEM TIMCOM_SRC "timcom_cpl.F90" "timcom_output.F90" "timcom_omip.F90")
if((NOT ${USE_ARM}) AND ${USE_GPU})
    file(GLOB TIMCOM_GPU_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "gpu/*.F90")
    list(APPEND TIMCOM_SRC ${TIMCOM_GPU_SRC})
endif()
add_library(timcom STATIC ${TIMCOM_SRC})

target_compile_options(timcom PUBLIC ${TIMCOM_OPTION})
target_include_directories(timcom PUBLIC ${TIMCOM_INC_PATH})
target_link_options(timcom PUBLIC ${TIMCOM_OPTION})
target_link_directories(timcom PUBLIC ${TIMCOM_LIB_PATH})
target_link_libraries(timcom ${TIMCOM_LINK_LIB})

if(${USE_OMIP})
    add_executable(timcom_omip timcom_omip.F90)
    target_link_libraries(timcom_omip PUBLIC timcom)
endif()

if(NOT ${USE_OMIP} AND ${USE_CLM})
    # 條件：USE_OMIP 為 OFF，USE_CLM 為 ON
    message(STATUS "Condition: USE_OMIP=${USE_OMIP} AND USE_CLM=${USE_CLM}")
    target_compile_options(timcom PUBLIC -Dowdate -Dclm_r -Dcpl -Dcpl_cice)
    message("USE FLAG : -Dowdate -Dclm_r -Dcpl -Dcpl_cice")
elseif(NOT ${USE_OMIP} AND NOT ${USE_CLM})
    # 條件：USE_OMIP 為 OFF，USE_CLM 為 OFF
    message(STATUS "Condition: USE_OMIP=${USE_OMIP} AND USE_CLM=${USE_CLM}")
    target_compile_options(timcom PUBLIC -Dcpl)
endif()

target_include_directories(timcom INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
