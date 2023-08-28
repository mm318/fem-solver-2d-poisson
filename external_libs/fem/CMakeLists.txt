cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME fem_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME fem)

set(LIBRARY_HEADERS
    fem2d_poisson.h
    prob_def.h
    fem_solver_debug.h
)

set(LIBRARY_SOURCE
    dgb.c
    fem2d_poisson.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCE}
)

target_include_directories(${LIBRARY_NAME} BEFORE
    PRIVATE "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}"
)

set_target_properties(${LIBRARY_NAME} PROPERTIES
    PUBLIC_HEADER "${LIBRARY_HEADERS}"
)

install(TARGETS ${LIBRARY_NAME}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${LIBRARY_NAME}"
)
