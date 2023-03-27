cmake_minimum_required(VERSION 3.0.0)

set(FEM_SOLVER_PROJECT_NAME
    FEM_Solver
)
set(FEM_SOLVER_SRC_DIR
    modules/main
)
set(FEM_SOLVER_SRC
    ${FEM_SOLVER_SRC_DIR}/fem_solver.c
)
set(FEM_SOLVER_LIBS_DIR
    modules
)
set(FEM_SOLVER_LIBS 
    fem
    libgd
    libpng
    libpng_main
    plplot
    plplot_main
    triangle
    triangle_main

    zlib
)

project(${FEM_SOLVER_PROJECT_NAME})

###########
# Project #
###########
add_executable(${FEM_SOLVER_PROJECT_NAME} ${FEM_SOLVER_SRC})

foreach(LIBRARY ${FEM_SOLVER_LIBS})
    add_subdirectory("${FEM_SOLVER_LIBS_DIR}/${LIBRARY}")
endforeach(LIBRARY)

set_target_properties(${FEM_SOLVER_PROJECT_NAME} ${FEM_SOLVER_LIBS}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

target_link_libraries(${FEM_SOLVER_PROJECT_NAME} ${FEM_SOLVER_LIBS})
