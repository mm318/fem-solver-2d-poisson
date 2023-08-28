cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME libgd_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME libgd)

set(LIBRARY_HEADERS
    bmp.h
    config.h
    entities.h
    gd.h
    gd_color.h
    gd_color_map.h
    gd_io.h
    gd_io_stream.h
    gd_nnquant.h
    gd_tga.h
    gdfontg.h
    gdfontl.h
    gdfontmb.h
    gdfonts.h
    gdfontt.h
    gdfx.h
    gdhelpers.h
    gdpp.h
    jisx0208.h
    wbmp.h
)

set(LIBRARY_SOURCES
    gd.c
    gd_bmp.c
    gd_color.c
    gd_color_map.c
    gd_crop.c
    gd_filter.c
    gd_gd.c
    gd_gd2.c
    gd_gif_in.c
    gd_gif_out.c
    gd_io.c
    gd_io_dp.c
    gd_io_file.c
    gd_io_ss.c
    gd_jpeg.c
    gd_nnquant.c
    gd_pixelate.c
    gd_png.c
    gd_security.c
    gd_ss.c
    gd_tga.c
    gd_tiff.c
    gd_topal.c
    gd_transform.c
    gd_wbmp.c
    gdfontg.c
    gdfontl.c
    gdfontmb.c
    gdfonts.c
    gdfontt.c
    gdft.c
    gdfx.c
    gdhelpers.c
    gdkanji.c
    gdtables.c
    gdxpm.c
    wbmp.c
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
