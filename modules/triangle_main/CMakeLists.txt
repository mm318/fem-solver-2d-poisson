cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME
    triangle_main_Library
)
set(LIBRARY_NAME
    triangle_main
)
set(LIBRARY_HEADERS_DIR
    .
)
set(LIBRARY_HEADERS
    ${LIBRARY_HEADERS_DIR}/get_mesh.h
    ${LIBRARY_HEADERS_DIR}/mesh_vis.h
    ${LIBRARY_HEADERS_DIR}/node_input.h
)
set(LIBRARY_SOURCE_DIR
    .
)
set(LIBRARY_SOURCE
    ${LIBRARY_HEADERS_DIR}/get_mesh.c
    ${LIBRARY_HEADERS_DIR}/mesh_vis.c
    ${LIBRARY_HEADERS_DIR}/node_input.c
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

add_subdirectory(tests)
