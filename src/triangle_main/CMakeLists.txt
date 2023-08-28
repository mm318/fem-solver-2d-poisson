cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME triangle_main_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME triangle_main)

set(LIBRARY_HEADERS
    get_mesh.h
    mesh_vis.h
    node_input.h
)

set(LIBRARY_SOURCES
    get_mesh.c
    mesh_vis.c
    node_input.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
)

target_include_directories(${LIBRARY_NAME} BEFORE
    PRIVATE "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

add_dependencies(${LIBRARY_NAME} triangle)

add_subdirectory(tests)
