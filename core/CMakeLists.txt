cmake_minimum_required (VERSION 2.8)
project (pydmrg)

set(BLOCKDIR $ENV{HOME}/Software/src/QC/Block)
message("Block src dir is  ${BLOCKDIR}\n")

set(LAPACK_LIBRARIES "-L${MKL_LIB_PATH} -I${MKLINCLUDE} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_avx -lm")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LIBRARIES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fPIC -O0")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(
  ${BLOCKDIR}
  ${BLOCKDIR}/include
  ${BLOCKDIR}/newmat10
  ${BLOCKDIR}/modules/twopdm
  ${BLOCKDIR}/modules/onepdm
  ${BLOCKDIR}/modules/generate_blocks)

find_package(OpenMP)
if (OPENMP_FOUND)
  set(_OPENMP 1)
endif()

# without MPI
set(SERIAL 1)

set(FAST_MTP 1)

find_package(PythonLibs REQUIRED)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
  OUTPUT_VARIABLE NUMPY_INCLUDE)
include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE})

set(Boost_USE_STATIC_LIBS off)
set(Boost_USE_MULTITHREADED off)
#target_link_libraries(${Boost_LIBRARY_DIRS})
set(Boost_COMPONENTS serialization filesystem system)
find_package(Boost 1.47.0 COMPONENTS ${Boost_COMPONENTS})
include_directories(${Boost_INCLUDE_DIRS})

configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h")

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

link_directories(
  ${BLOCKDIR}/newmat10
  ${BLOCKDIR})

add_library(_dmrg SHARED 
  _dmrg.cpp quanta.cc rotationmat.cc wavefn.cc stateinfo.cc spinblock.cc)
set_target_properties(_dmrg PROPERTIES PREFIX ""
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../)

target_link_libraries(_dmrg
  qcdmrg newmat
  ${Boost_LIBRARIES} ${LAPACK_LIBRARIES})

