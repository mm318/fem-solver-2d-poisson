cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(TEST_PROJECT_NAME triangle_main_Tests)
project(${TEST_PROJECT_NAME})

set(LIBRARY_TESTS_SOURCES
    test_triangle.c
)

add_executable(${TEST_PROJECT_NAME}
    ${LIBRARY_TESTS_SOURCES}
)

target_include_directories(${TEST_PROJECT_NAME} BEFORE
    PRIVATE "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
)

target_link_libraries(${TEST_PROJECT_NAME}
    ${LIBRARY_NAME} # NOTE: This is defined from project above
    libtriangle.a
    m
)

set_target_properties(${TEST_PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/tests"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/tests"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
)

enable_testing()

add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})
