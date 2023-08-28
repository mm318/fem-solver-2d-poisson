cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME libpng_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME libpng)

set(LIBRARY_HEADERS
    png.h
    pngconf.h
    pngdebug.h
    pnginfo.h
    pnglibconf.h
    pngpriv.h
    pngstruct.h
)

set(LIBRARY_SOURCES
    png.c
    pngerror.c
    pngget.c
    pngmem.c
    pngpread.c
    pngread.c
    pngrio.c
    pngrtran.c
    pngrutil.c
    pngset.c
    pngtrans.c
    pngwio.c
    pngwrite.c
    pngwtran.c
    pngwutil.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
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
