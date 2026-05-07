# ##############################################################################
# Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

file(GLOB TIMCOM_TEST_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "test_*.F90")
string(REPLACE ".F90" ";" TEST_LIST ${TIMCOM_TEST_SRC})

foreach(TEST ${TEST_LIST})
    add_executable(${TEST} ${TEST}.F90 helper.F90)
    add_test(
        NAME ${TEST}
        COMMAND mpiexec -n ${TEST_MPI_NP} ${TEST} 000)
    target_link_libraries(${TEST} PRIVATE timcom)
    add_dependencies(${TEST} timcom)
endforeach()
