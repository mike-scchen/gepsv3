# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

cmake_minimum_required(VERSION 3.20)
project(TIMCOM LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_BUILD_TYPE "Release")

# Override the automatic optimization flag by CMake
set(CMAKE_Fortran_FLAGS_RELEASE "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "")

option(USE_ARM "enable to build library or executable for ARM platform (FX1000)" OFF)
option(USE_OMIP "enable to generate the executable that only runs the ocean model." OFF)
option(USE_NVTX "enable to insert NVTX range (valid only when USE_ARAM is OFF)" OFF)
option(USE_GPU "enable to use gpu (valid only when USE_ARAM is OFF)" ON)
option(USE_CLM "The climate mode is dedicated, including restart" OFF)
option(USE_ICE "enable to couple with CICE model" OFF)
set(GPU_ARCH "80" CACHE STRING "GPU architecture (V100=70; A100=80; H100=90)")
set(CUDA_RUNTIME_VERSION "12.6" CACHE STRING "CUDA runtime version")
set(TEST_MPI_NP 8 CACHE STRING "Number of MPI processes for tests")


message(STATUS "Build Summary")
message(STATUS "USE_ARM: ${USE_ARM}")
message(STATUS "USE_OMIP: ${USE_OMIP}")
message(STATUS "USE_GPU: ${USE_GPU}")
message(STATUS "USE_CLM: ${USE_CLM}")
message(STATUS "USE_ICE: ${USE_ICE}")

if(NOT ${USE_ARM})
    message(STATUS "USE_NVTX: ${USE_NVTX}")
    message(STATUS "USE_GPU: ${USE_GPU}")
endif()

if(${USE_ARM})
    include(ConfigureARM)
else()
    include(ConfigureX86)
endif()

add_subdirectory(src)

if(${USE_OMIP} AND (NOT ${USE_ARM}) AND ${USE_GPU})
    include(CTest)
    add_subdirectory(test)
endif()
