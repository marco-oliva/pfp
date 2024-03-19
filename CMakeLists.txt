cmake_minimum_required(VERSION 3.15)

project(PFP++)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "3")
SET(VERSION_PATCH "9")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
set(PROJECT_VERSION ${VERSION})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
# Options

option(ENABLE_MALLOC_COUNT "Keep track of allocations with malloc_count."   OFF)
option(ENABLE_MIMALLOC "Use mimalloc if installed."                         OFF)
option(BUILD_TESTS "Build unit test."                                       OFF)
option(HTSlib_HIPERGATOR "Use hipergator default location for htslib."      OFF)
option(HTSlib_CUSTOM "Use custom location for htslib."                      OFF)

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
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.cpp")
configure_file(
        ${PROJECT_SOURCE_DIR}/src/version.cpp.in
        ${PROJECT_BINARY_DIR}/generated/version.cpp)

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

Include(FetchContent)

# Catch
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.0.1
)
# Making available only if building tests
FetchContent_GetProperties(Catch2)
if (NOT Catch2_POPULATED AND (BUILD_TESTS OR ${CMAKE_BUILD_TYPE} MATCHES "Debug"))
    FetchContent_MakeAvailable(Catch2)
endif()

# Hayai Benchmarks
FetchContent_Declare(
        hayai
        GIT_REPOSITORY https://github.com/nickbruun/hayai.git
        GIT_TAG        v1.0.2
)
FetchContent_GetProperties(hayai)
if(NOT hayai_POPULATED)
    FetchContent_Populate(hayai)
endif()

# mio mmap
FetchContent_Declare(
        mio
        GIT_REPOSITORY https://github.com/mandreyel/mio.git
)
FetchContent_GetProperties(mio)
if(NOT mio_POPULATED)
    FetchContent_Populate(mio)
endif()

# Murmur Hash
FetchContent_Declare(
        smhasher
        GIT_REPOSITORY https://github.com/aappleby/smhasher.git
)
FetchContent_GetProperties(smhasher)
if(NOT smhasher_POPULATED)
    FetchContent_Populate(smhasher)
endif()
set(MURMUR_SOURCES ${smhasher_SOURCE_DIR}/src/MurmurHash3.cpp)

# CLI11
FetchContent_Declare(
        cli11
        GIT_REPOSITORY "https://github.com/CLIUtils/CLI11.git"
        GIT_TAG "v1.9.0"
)
FetchContent_GetProperties(cli11)
if(NOT cli11_POPULATED)
    FetchContent_Populate(cli11)
endif()

# SPDLOG
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.10.0
)
FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
endif()

# Zstr
FetchContent_Declare(
        zstr
        GIT_REPOSITORY https://github.com/mateidavid/zstr.git
        GIT_TAG        v1.0.6
)
FetchContent_GetProperties(zstr)
if(NOT zstr_POPULATED)
    FetchContent_Populate(zstr)
endif()

# malloc count
FetchContent_Declare(
        malloc_count
        GIT_REPOSITORY https://github.com/bingmann/malloc_count.git
        GIT_TAG        v0.7.1
)
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
    FetchContent_Populate(malloc_count)
    add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
    target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")
    target_link_libraries(malloc_count dl)
endif()

# kseq
FetchContent_Declare(
        kseq
        GIT_REPOSITORY https://github.com/lh3/seqtk.git
        GIT_TAG        v1.3
)
FetchContent_GetProperties(kseq)
if(NOT kseq_POPULATED)
    FetchContent_Populate(kseq)
endif()

# Openmp
find_package(OpenMP)
if(NOT OpenMP_FOUND AND NOT OPENMP_FOUND AND APPLE)
    # Before we can hunt around for OpenMP on Apple we need to know where
    # Homebrew is, if installed.
    set(HOMEBREW_PREFIX "$ENV{HOMEBREW_PREFIX}")
    if("${HOMEBREW_PREFIX}" STREQUAL "")
        # Ask Homebrew if it exists and where it is.
        # Returns an empty string if this command can't run.
        execute_process(COMMAND brew --prefix OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE HOMEBREW_PREFIX)
    endif()
    if("${HOMEBREW_PREFIX}" STREQUAL "")
        # Use defaults
        if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm")
            set(HOMEBREW_PREFIX "/opt/homebrew")
        else()
            set(HOMEBREW_PREFIX "/usr/local")
        endif()
    endif()
endif()
if(OpenMP_FOUND OR OPENMP_FOUND)
    message(STATUS "HAVE_OPENMP: OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(OMP_LIBRARIES "${OpenMP_CXX_LIBRARIES}")
elseif(APPLE AND EXISTS "/opt/local/lib/libomp")
    # official Apple compiler with macports' libomp
    message(STATUS "HAVE_OPENMP: OpenMP[MacPorts]")
    set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp -I/opt/local/include/libomp")
    set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/local/include/libomp")
    set(OMP_LIBRARIES "-L/opt/local/lib/libomp -lomp")
elseif(APPLE AND EXISTS "${HOMEBREW_PREFIX}/opt/libomp")
    # official Apple compiler with homebrew's kegg libomp
    message(STATUS "HAVE_OPENMP: OpenMP[Homebrew kegg]")
    set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
    set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
    set(OMP_LIBRARIES "-L${HOMEBREW_PREFIX}/opt/libomp/lib -lomp")
elseif(APPLE AND EXISTS "/usr/local/lib/libomp.dylib")
    # official Apple compiler with homebrew's libomp
    message(STATUS "HAVE_OPENMP: OpenMP[Homebrew]")
    set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
    set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set(OMP_LIBRARIES "-L/usr/local/lib -lomp")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND EXISTS "/usr/local/lib64/libgomp.a")
    # custom static openmp for official linux builds
    message(STATUS "HAVE_OPENMP: OpenMP[static]")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set(OMP_LIBRARIES "/usr/local/lib64/libgomp.a")
else()
    message(STATUS "HAVE_OPENMP: NOT FOUND")
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
    message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
elseif(HTSlib_CUSTOM)
    set(HTSlib_INCLUDE_DIRS "")
    set(HTSlib_LIBRARIES "")
    set(LIBLZMA_LIBRARIES "")
    set(BZIP2_LIBRARIES "")
    set(ZLIB_LIBRARIES "")
    set(CURL_LIBRARIES "")
    set(HTSlib_LIBRARIES ${HTSlib_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${CURL_LIBRARIES})
    message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
else()
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
    find_package(HTSlib REQUIRED)
    if(HTSlib_FOUND)
        message(STATUS "HTSlib: ${HTSlib_LIBRARIES}")
    endif()
endif()

# MiMalloc
if (ENABLE_MIMALLOC)
    set(MI_BUILD_STATIC OFF)
    set(MI_BUILD_OBJECT OFF)
    set(MI_BUILD_TESTS  OFF)
    FetchContent_Declare(
            mimalloc
            GIT_REPOSITORY https://github.com/microsoft/mimalloc.git
            GIT_TAG        v2.0.9
    )
    FetchContent_GetProperties(mimalloc)
    if (NOT mimalloc_POPULATED)
        FetchContent_MakeAvailable(mimalloc)
        message(STATUS "Building mimalloc. It will be dynamically linked.")
    endif()
    set(MIMALLOC_D_LIB mimalloc)
endif()

################################################################################
# Lib

set(PFP_INCLUDES
        ${mio_SOURCE_DIR}/single_include
        ${smhasher_SOURCE_DIR}/src
        ${HTSlib_INCLUDE_DIRS}
        ${cli11_SOURCE_DIR}/include
        ${spdlog_SOURCE_DIR}/include
        ${zstr_SOURCE_DIR}/src
        ${kseq_SOURCE_DIR}
        ${hayai_SOURCE_DIR}/src)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PFP_INCLUDES})

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB SOURCES_GENERATED ${PROJECT_BINARY_DIR}/generated/*.cpp)

# 32 and 64 bits parse
add_library(pfp ${SOURCES} ${SOURCES_GENERATED} ${MURMUR_SOURCES})
add_library(pfp64 ${SOURCES} ${SOURCES_GENERATED} ${MURMUR_SOURCES})
target_compile_options(pfp64 PRIVATE "-DPFP_LONG_TYPE=ON")

set(PFP_PUBLIC_HEADERS
        ${PROJECT_SOURCE_DIR}/include/internals.hpp
        ${PROJECT_SOURCE_DIR}/include/pfp_algo.hpp
        ${PROJECT_SOURCE_DIR}/include/utils.hpp
        ${PROJECT_SOURCE_DIR}/include/vcf.hpp
        ${PROJECT_SOURCE_DIR}/include/version.hpp)

set_target_properties(pfp PROPERTIES PUBLIC_HEADER "${PFP_PUBLIC_HEADERS}")
set_target_properties(pfp64 PROPERTIES PUBLIC_HEADER "${PFP_PUBLIC_HEADERS}")

# Set libraries
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(VCF_LIB_DEPS ${HTSlib_LIBRARIES} dl pthread ssl crypto)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(VCF_LIB_DEPS ${HTSlib_LIBRARIES} dl pthread)
endif()

if (ENABLE_MIMALLOC AND ENABLE_MALLOC_COUNT)
    message(STATUS "mimalloc enabled, disabling malloc count")
    set(ENABLE_MALLOC_COUNT FALSE)
endif()

if (ENABLE_MALLOC_COUNT)
    message(STATUS "Using malloc_count")
    set(MALLOC_COUNT_LIBRARIES malloc_count)
else()
    set(MALLOC_COUNT_LIBRARIES "")
endif()

set(PFP_LIBS ${VCF_LIB_DEPS} ${OMP_LIBRARIES} ${MALLOC_COUNT_LIBRARIES} ${MIMALLOC_D_LIB})

################################################################################
# Tests

if (BUILD_TESTS OR ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(STATUS "Building unit tests")
    add_executable(pfp_unit_tests tests/unit_tests.cpp)
    target_link_libraries(pfp_unit_tests pfp ${PFP_LIBS} Catch2::Catch2)

    add_executable(pfp_unit_tests64 tests/unit_tests.cpp)
    target_compile_options(pfp_unit_tests64 PRIVATE "-DPFP_LONG_TYPE=ON")
    target_link_libraries(pfp_unit_tests64 pfp64 ${PFP_LIBS} Catch2::Catch2)

    enable_testing()
    add_test(NAME pfp_unit_tests COMMAND pfp_unit_tests)
    add_test(NAME pfp_unit_tests64 COMMAND pfp_unit_tests64)

    message(STATUS "Building benchmarks")
    add_executable(pfp_benchmarks benchmarks/kr_hash.cpp)
    target_link_libraries(pfp_benchmarks pfp ${PFP_LIBS})
endif()

################################################################################
# Tools

# PFP++

# 32 and 64 bits parse
add_executable(pfp++ pfp++.cpp)
target_link_libraries(pfp++ pfp ${PFP_LIBS})
add_executable(pfp++64 pfp++.cpp)
target_compile_options(pfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(pfp++64 pfp64 ${PFP_LIBS})

# PFP++ merge
add_executable(mpfp++ merge.cpp)
target_link_libraries(mpfp++ pfp ${PFP_LIBS})
add_executable(mpfp++64 merge.cpp)
target_compile_options(mpfp++64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(mpfp++64 pfp64 ${PFP_LIBS})

# VCT to Fasta
add_executable(vcf_to_fa vcf_to_fa.cpp)
target_link_libraries(vcf_to_fa pfp ${PFP_LIBS})

# Integrity Check
add_executable(check check_integrity.cpp)
target_link_libraries(check pfp ${PFP_LIBS})
add_executable(check64 check_integrity.cpp)
target_compile_options(check64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(check64 pfp64 ${PFP_LIBS})

# Extract properties from precomputed pfp
add_executable(exprop exprop.cpp)
target_link_libraries(exprop pfp ${PFP_LIBS})
add_executable(exprop64 exprop.cpp)
target_compile_options(exprop64 PRIVATE "-DPFP_LONG_TYPE=ON")
target_link_libraries(exprop64 pfp64 ${PFP_LIBS})

################################################################################
# Install

install(TARGETS pfp pfp64 pfp++ pfp++64 mpfp++ mpfp++64 vcf_to_fa check check64 exprop exprop64
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        PUBLIC_HEADER DESTINATION lib/static/include)
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





