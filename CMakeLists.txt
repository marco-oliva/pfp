cmake_minimum_required(VERSION 3.5.0)

project(vcf-bwt)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

################################################################################
# Version number
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hpp")
configure_file(
        ${CMAKE_SOURCE_DIR}/include/version.hpp.in
        ${CMAKE_BINARY_DIR}/generated/version.hpp)

include_directories(${CMAKE_BINARY_DIR}/generated)

################################################################################
# Compiler and linker flags

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")

################################################################################
# Check dependencies

# HTS
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(HTSlib REQUIRED)
if(HTSlib_FOUND)
    include_directories(${HTSlib_INCLUDE_DIRS})
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    link_libraries(${OpenMP_CXX_LIBRARIES})
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # assumes clang build
        # we can't reliably detect when we're using clang, so for the time being we assume
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

################################################################################
# Lib

include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/generated)

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
file(GLOB MALLOC_COUNT_SOURCES ${CMAKE_SOURCE_DIR}/thirdparty/malloc_count/*.c)
file(GLOB MURMUR_SOURCES ${CMAKE_SOURCE_DIR}/thirdparty/murmur/*.cpp)

# 32 and 64 bits parse
add_library(pfp ${SOURCES} ${MALLOC_COUNT_SOURCES} ${MURMUR_SOURCES})
add_library(pfp64 ${SOURCES} ${MALLOC_COUNT_SOURCES} ${MURMUR_SOURCES})
target_compile_options(pfp64 PRIVATE "-DPFP_LONG_TYPE=ON")

set(VCF_LIB_DEPS "z" "${HTSlib_LIBRARIES}" "dl")

################################################################################
# Tests

add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests pfp ${VCF_LIB_DEPS})

add_executable(unit_tests64 tests/unit_tests.cpp)
target_compile_options(unit_tests64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(unit_tests64 pfp64 ${VCF_LIB_DEPS})

################################################################################
# Tools

# PFP++

# 32 and 64 bits parse
add_executable(pfp++ pfp++.cpp)
target_link_libraries(pfp++ pfp ${VCF_LIB_DEPS})
add_executable(pfp++64 pfp++.cpp)
target_compile_options(pfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(pfp++64 pfp64 ${VCF_LIB_DEPS})

# AuPair
add_executable(AuPair au_pair.cpp)
target_link_libraries(AuPair pfp ${VCF_LIB_DEPS})
add_executable(AuPair64 au_pair.cpp)
target_link_libraries(AuPair64 pfp64 ${VCF_LIB_DEPS})

# PFP++ merge
add_executable(mpfp++ merge.cpp)
target_link_libraries(mpfp++ pfp ${VCF_LIB_DEPS})
add_executable(mpfp++64 merge.cpp)
target_compile_options(mpfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(mpfp++64 pfp64 ${VCF_LIB_DEPS})

# VCT to Fasta
add_executable(vcf_to_fa vcf_to_fa.cpp)
target_link_libraries(vcf_to_fa pfp ${VCF_LIB_DEPS})

# Unparse
add_executable(unparse unparse.cpp)
target_link_libraries(unparse pfp ${VCF_LIB_DEPS})
add_executable(unparse64 unparse.cpp)
target_compile_options(unparse64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(unparse64 pfp64 ${VCF_LIB_DEPS})

# Integrity Check
add_executable(check check_integrity.cpp)
target_link_libraries(check pfp ${VCF_LIB_DEPS})
add_executable(check64 check_integrity.cpp)
target_compile_options(check64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(check64 pfp64 ${VCF_LIB_DEPS})







