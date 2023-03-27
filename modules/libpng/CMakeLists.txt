cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME
    libpng_Library
)
set(LIBRARY_NAME
    libpng
)
set(LIBRARY_HEADERS_DIR
    .
)
set(LIBRARY_HEADERS
    ${LIBRARY_HEADERS_DIR}/png.h
    ${LIBRARY_HEADERS_DIR}/pngconf.h
    ${LIBRARY_HEADERS_DIR}/pngdebug.h
    ${LIBRARY_HEADERS_DIR}/pnginfo.h
    ${LIBRARY_HEADERS_DIR}/pnglibconf.h
    ${LIBRARY_HEADERS_DIR}/pngpriv.h
    ${LIBRARY_HEADERS_DIR}/pngstruct.h
)
set(LIBRARY_SOURCE_DIR
    .
)
set(LIBRARY_SOURCE
    ${LIBRARY_SOURCE_DIR}/png.c
    ${LIBRARY_SOURCE_DIR}/pngerror.c
    ${LIBRARY_SOURCE_DIR}/pngget.c
    ${LIBRARY_SOURCE_DIR}/pngmem.c
    ${LIBRARY_SOURCE_DIR}/pngpread.c
    ${LIBRARY_SOURCE_DIR}/pngread.c
    ${LIBRARY_SOURCE_DIR}/pngrio.c
    ${LIBRARY_SOURCE_DIR}/pngrtran.c
    ${LIBRARY_SOURCE_DIR}/pngrutil.c
    ${LIBRARY_SOURCE_DIR}/pngset.c
    ${LIBRARY_SOURCE_DIR}/pngtrans.c
    ${LIBRARY_SOURCE_DIR}/pngwio.c
    ${LIBRARY_SOURCE_DIR}/pngwrite.c
    ${LIBRARY_SOURCE_DIR}/pngwtran.c
    ${LIBRARY_SOURCE_DIR}/pngwutil.c
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
