cmake_minimum_required(VERSION 3.0.0)

set(TEST_PROJECT_NAME
    plplot_main_Tests
)
set(LIBRARY_TESTS_SOURCE
    test_plplot.c
)

project(${TEST_PROJECT_NAME})

enable_testing()

add_executable(${TEST_PROJECT_NAME} ${LIBRARY_TESTS_SOURCE})

set_target_properties(${TEST_PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/tests"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/tests"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
)

target_link_libraries(${TEST_PROJECT_NAME}
    plplot
    libgd
    libpng
    zlib
    m
    ${LIBRARY_NAME} # NOTE: This is defined from project above
)

add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})
