cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME plplot_main_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME plplot_main)

set(LIBRARY_HEADERS
    interpolated_plot.h
)

set(LIBRARY_SOURCES
    interpolated_plot.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
)

target_include_directories(${LIBRARY_NAME} BEFORE
    PRIVATE "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
)

add_dependencies(${LIBRARY_NAME} plplot)

add_subdirectory(tests)
