cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME libpng_main_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME libpng_main)

set(LIBRARY_HEADERS
    png_writer.h
)

set(LIBRARY_SOURCES
    png_writer.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
)

add_dependencies(${LIBRARY_NAME} libpng)