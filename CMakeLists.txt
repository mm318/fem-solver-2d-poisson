cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(FEM_SOLVER_PROJECT_NAME FEM_Solver)
project(${FEM_SOLVER_PROJECT_NAME})


#########################
# External Dependencies #
#########################

include(ExternalProject)

ExternalProject_Add(zlib
    URL "${CMAKE_SOURCE_DIR}/external_libs/zlib.tar.bz2"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

ExternalProject_Add(triangle
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/external_libs/triangle"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

ExternalProject_Add(libpng
    DEPENDS zlib
    URL "${CMAKE_SOURCE_DIR}/external_libs/libpng.tar.bz2"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

ExternalProject_Add(libgd
    DEPENDS libpng
    URL "${CMAKE_SOURCE_DIR}/external_libs/libgd.tar.bz2"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

ExternalProject_Add(plplot
    DEPENDS libgd
    URL "${CMAKE_SOURCE_DIR}/external_libs/plplot.tar.bz2"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

ExternalProject_Add(fem
    DEPENDS triangle
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/external_libs/fem"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
)

link_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")


###########
# Project #
###########

set(FEM_SOLVER_LIBS 
    triangle_main
    plplot_main
    libpng_main
)

set(FEM_SOLVER_SOURCES
    src/main/fem_solver.c
    src/main/prob_def.c
)

foreach(LIBRARY ${FEM_SOLVER_LIBS})
    add_subdirectory("src/${LIBRARY}")
endforeach()

add_executable(${FEM_SOLVER_PROJECT_NAME}
    ${FEM_SOLVER_SOURCES}
)

target_include_directories(${FEM_SOLVER_PROJECT_NAME} BEFORE
    PRIVATE "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
            "src"
)

target_link_libraries(${FEM_SOLVER_PROJECT_NAME}
    ${FEM_SOLVER_LIBS}
    libfem.a
    libplplot.a
    liblibgd.a
    liblibpng.a
    libtriangle.a
    libzlib.a
    m
)

set_target_properties(${FEM_SOLVER_PROJECT_NAME} ${FEM_SOLVER_LIBS}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

file(COPY bin/ DESTINATION bin/)
