project(landau_hydro CXX)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_landauhydro COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/landauhydro/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/landauhydro/software/include ${GITHOME_MSU}/commonutils/software/include )

link_directories(${GITHOME_MSU}/landauhydro/software/lib ${GITHOME_MSU}/commonutils/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/landauhydro/software/lib landau_hydro)

#################################################

set(SOURCE_FILES ./hydromain.cc)
SET_SOURCE_FILES_PROPERTIES(hydromain.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a
)
SET_SOURCE_FILES_PROPERTIES(hydroman.cc
	PROPERTIES OBJECT_DEPENDS  ${GITHOME_MSU}/landauhydro/software/lib/liblandau_hydro.a
)
project(hydro CXX)
add_executable(hydro ${SOURCE_FILES})
add_dependencies(hydro extern_landauhydro extern_commonutils)
target_link_libraries(hydro ${GITHOME_MSU}/landauhydro/software/lib/liblandau_hydro.a ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################

