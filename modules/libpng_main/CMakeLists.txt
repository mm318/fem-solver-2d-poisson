cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME
    libpng_main_Library
)
set(LIBRARY_NAME
    libpng_main
)
set(LIBRARY_HEADERS_DIR
    .
)
set(LIBRARY_HEADERS
    ${LIBRARY_HEADERS_DIR}/png_writer.h
)
set(LIBRARY_SOURCE_DIR
    .
)
set(LIBRARY_SOURCE
    ${LIBRARY_SOURCE_DIR}/png_writer.c
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
