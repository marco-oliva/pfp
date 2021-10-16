# - Try to find HTSlib
# Once done, this will define
#
#  HTSlib_FOUND - system has HTSlib
#  HTSlib_INCLUDE_DIRS - the HTSlib include directories
#  HTSlib_LIBRARIES - link these to use HTSlib

include(LibFindMacros)

# HTSlib dependencies
libfind_package(HTSLib ZLIB REQUIRED)
libfind_package(HTSLib LibLZMA REQUIRED)
libfind_package(HTSLib BZip2 REQUIRED)

include_directories(${LIBLZMA_INCLUDE_DIRS})
include_directories(${BZIP2_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS})

#target_link_libraries(DukasCompiler ${LIBLZMA_LIBRARIES})

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(HTSlib_PKGCONF HTSlib)

# Locate include dir
find_path(HTSlib_INCLUDE_DIR
        NAMES hts.h sam.h
        PATHS ${HTSlib_PKGCONF_INCLUDE_DIRS}
        PATH_SUFFIXES include include/htslib
        )

# Find the library itself
find_library(HTSlib_LIBRARY
        NAMES hts
        PATHS ${HTSlib_PKGCONF_LIBRARY_DIRS}
        )

libfind_process(HTSlib)

set(HTSlib_LIBRARIES ${HTSlib_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES})