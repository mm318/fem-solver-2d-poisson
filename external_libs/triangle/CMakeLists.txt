cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME triangle_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME triangle)

set(LIBRARY_HEADERS
    triangle.h
)

set(LIBRARY_SOURCES
    triangle.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
)

set_target_properties(${LIBRARY_NAME} PROPERTIES
    PUBLIC_HEADER "${LIBRARY_HEADERS}"
)

install(TARGETS ${LIBRARY_NAME}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${LIBRARY_NAME}"
)
