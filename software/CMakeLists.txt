cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(landau_hydro CXX C)
set(GITHOME_MSU $ENV{GITHOME_MSU})
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} include ${GITHOME_MSU}/commonutils/software/include)

file(GLOB landau_hydro_src ${CMAKE_CURRENT_SOURCE_DIR}/src/**.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(landau_hydro STATIC ${landau_hydro_src})
target_link_libraries(landau_hydro ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
message("Libraries: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}")



