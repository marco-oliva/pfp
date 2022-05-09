cmake_minimum_required(VERSION 3.15)

project(PFP++)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "3")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
set(PROJECT_VERSION ${VERSION})

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

# HTS
if(HTSlib_HIPERGATOR)
    set(HTSlib_INCLUDE_DIRS "/apps/htslib/1.15/include/htslib")
    set(HTSlib_LIBRARIES "/apps/htslib/1.15/lib/libhts.a")
    set(LIBLZMA_LIBRARIES "/apps/cmake/3.21.3/lib/liblzma.so")
    set(BZIP2_LIBRARIES "/apps/cmake/3.21.3/lib/libbz2.so")
    set(ZLIB_LIBRARIES "/apps/cmake/3.21.3/lib/libz.so")
    set(CURL_LIBRARIES "/apps/cmake/3.21.3/lib/libcurl.so")
    set(CRYPTO_LIBRARIES "/apps/cmake/3.21.3/lib/libcrypto.so")
    set(HPG_INCLUDES "/apps/cmake/3.21.3/include")
    set(HTSlib_INCLUDE_DIRS ${HTSlib_INCLUDE_DIRS} ${HPG_INCLUDES})
    set(HTSlib_LIBRARIES ${HTSlib_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${CURL_LIBRARIES} ${CRYPTO_LIBRARIES})
    include_directories(${HTSlib_INCLUDE_DIRS})
    message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
elseif(HTSlib_CUSTOM)
    set(HTSlib_INCLUDE_DIRS "")
    set(HTSlib_LIBRARIES "")
    set(LIBLZMA_LIBRARIES "")
    set(BZIP2_LIBRARIES "")
    set(ZLIB_LIBRARIES "")
    set(CURL_LIBRARIES "")
    set(HTSlib_LIBRARIES ${HTSlib_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${CURL_LIBRARIES})
    include_directories(${HTSlib_INCLUDE_DIRS})
    message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
else()
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
    find_package(HTSlib REQUIRED)
    if(HTSlib_FOUND)
        include_directories(${HTSlib_INCLUDE_DIRS})
        message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
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

set(VCF_LIB_DEPS ${HTSlib_LIBRARIES} dl pthread ssl crypto m)

################################################################################
# Tests

add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests pfp ${VCF_LIB_DEPS})

add_executable(unit_tests64 tests/unit_tests.cpp)
target_compile_options(unit_tests64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(unit_tests64 pfp64 ${VCF_LIB_DEPS})

enable_testing()
add_test(NAME unit_tests COMMAND unit_tests)
add_test(NAME unit_tests64 COMMAND unit_tests64)
################################################################################
# Tools

# PFP++

# 32 and 64 bits parse
add_executable(pfp++ pfp++.cpp)
target_link_libraries(pfp++ pfp ${VCF_LIB_DEPS})
add_executable(pfp++64 pfp++.cpp)
target_compile_options(pfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(pfp++64 pfp64 ${VCF_LIB_DEPS})

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

################################################################################
# Install

install(TARGETS pfp++ pfp++64 mpfp++ mpfp++64 vcf_to_fa unparse unparse64 check check64)

################################################################################
# CPack

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${VERSION}")

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "STGZ;TGZ;DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Florida")
set(CPACK_PACKAGE_CONTACT "marco.oliva@ufl.edu")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "PFP++")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-sources")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Marco Oliva")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE)
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEB_COMPONENT_INSTALL YES)
include(CPack)





