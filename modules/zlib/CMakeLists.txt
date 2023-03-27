cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME
    zlib_Library
)
set(LIBRARY_NAME
    zlib
)
set(LIBRARY_HEADERS_DIR
    .
)
set(LIBRARY_HEADERS
    ${LIBRARY_HEADERS_DIR}/crc32.h
    ${LIBRARY_HEADERS_DIR}/deflate.h
    ${LIBRARY_HEADERS_DIR}/gzguts.h
    ${LIBRARY_HEADERS_DIR}/inffast.h
    ${LIBRARY_HEADERS_DIR}/inffixed.h
    ${LIBRARY_HEADERS_DIR}/inflate.h
    ${LIBRARY_HEADERS_DIR}/inftrees.h
    ${LIBRARY_HEADERS_DIR}/trees.h
    ${LIBRARY_HEADERS_DIR}/zconf.h
    ${LIBRARY_HEADERS_DIR}/zlib.h
    ${LIBRARY_HEADERS_DIR}/zutil.h
)
set(LIBRARY_SOURCE_DIR
    .
)
set(LIBRARY_SOURCE
    ${LIBRARY_SOURCE_DIR}/adler32.c
    ${LIBRARY_SOURCE_DIR}/compress.c
    ${LIBRARY_SOURCE_DIR}/crc32.c
    ${LIBRARY_SOURCE_DIR}/deflate.c
    ${LIBRARY_SOURCE_DIR}/gzclose.c
    ${LIBRARY_SOURCE_DIR}/gzlib.c
    ${LIBRARY_SOURCE_DIR}/gzread.c
    ${LIBRARY_SOURCE_DIR}/gzwrite.c
    ${LIBRARY_SOURCE_DIR}/infback.c
    ${LIBRARY_SOURCE_DIR}/inffast.c
    ${LIBRARY_SOURCE_DIR}/inflate.c
    ${LIBRARY_SOURCE_DIR}/inftrees.c
    ${LIBRARY_SOURCE_DIR}/trees.c
    ${LIBRARY_SOURCE_DIR}/uncompr.c
    ${LIBRARY_SOURCE_DIR}/zutil.c
)

project(${PROJECT_NAME})

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_HEADERS}
    ${LIBRARY_SOURCE}
)

target_include_directories(${LIBRARY_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${LIBRARY_HEADERS_DIR}>
    $<INSTALL_INTERFACE:${LIBRARY_HEADERS_DIR}>
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:..>
)
