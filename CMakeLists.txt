cmake_minimum_required(VERSION 3.12)
project(MInDes VERSION 0.1.0 LANGUAGES C CXX)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(CTest)
enable_testing()

file(GLOB_RECURSE FFTW_LIB lib/*.a)
link_libraries(${FFTW_LIB})

file(GLOB_RECURSE MODULES_SRC modules/*.cpp modules/*.c)
file(GLOB_RECURSE MODULES_HEAD modules/*.h modules/*.hh)
SET(MODULE_FILE ${MODULES_SRC} ${MODULES_HEAD})

file(GLOB_RECURSE SOLVERS_SRC solvers/*.cpp solvers/*.c)
file(GLOB_RECURSE SOLVERS_HEAD solvers/*.h solvers/*.hh)
SET(SOLVER_FILE ${SOLVERS_SRC} ${SOLVERS_HEAD})

SET(TARGET_ALL  ${MODULE_FILE} ${SOLVER_FILE} MInDes_EXE.cpp)

add_executable(${CMAKE_PROJECT_NAME}  ${TARGET_ALL})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
