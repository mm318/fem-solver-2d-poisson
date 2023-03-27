cmake_minimum_required(VERSION 3.0.0)

set(PROJECT_NAME
    libgd_Library
)
set(LIBRARY_NAME
    libgd
)
set(LIBRARY_HEADERS_DIR
    .
)
set(LIBRARY_HEADERS
    ${LIBRARY_HEADERS_DIR}/bmp.h
    ${LIBRARY_HEADERS_DIR}/config.h
    ${LIBRARY_HEADERS_DIR}/entities.h
    ${LIBRARY_HEADERS_DIR}/gd.h
    ${LIBRARY_HEADERS_DIR}/gd_color.h
    ${LIBRARY_HEADERS_DIR}/gd_color_map.h
    ${LIBRARY_HEADERS_DIR}/gd_io.h
    ${LIBRARY_HEADERS_DIR}/gd_io_stream.h
    ${LIBRARY_HEADERS_DIR}/gd_nnquant.h
    ${LIBRARY_HEADERS_DIR}/gd_tga.h
    ${LIBRARY_HEADERS_DIR}/gdfontg.h
    ${LIBRARY_HEADERS_DIR}/gdfontl.h
    ${LIBRARY_HEADERS_DIR}/gdfontmb.h
    ${LIBRARY_HEADERS_DIR}/gdfonts.h
    ${LIBRARY_HEADERS_DIR}/gdfontt.h
    ${LIBRARY_HEADERS_DIR}/gdfx.h
    ${LIBRARY_HEADERS_DIR}/gdhelpers.h
    ${LIBRARY_HEADERS_DIR}/gdpp.h
    ${LIBRARY_HEADERS_DIR}/jisx0208.h
    ${LIBRARY_HEADERS_DIR}/wbmp.h
)
set(LIBRARY_SOURCE_DIR
    .
)
set(LIBRARY_SOURCE
    ${LIBRARY_SOURCE_DIR}/gd.c
    ${LIBRARY_SOURCE_DIR}/gd_bmp.c
    ${LIBRARY_SOURCE_DIR}/gd_color.c
    ${LIBRARY_SOURCE_DIR}/gd_color_map.c
    ${LIBRARY_SOURCE_DIR}/gd_crop.c
    ${LIBRARY_SOURCE_DIR}/gd_filter.c
    ${LIBRARY_SOURCE_DIR}/gd_gd.c
    ${LIBRARY_SOURCE_DIR}/gd_gd2.c
    ${LIBRARY_SOURCE_DIR}/gd_gif_in.c
    ${LIBRARY_SOURCE_DIR}/gd_gif_out.c
    ${LIBRARY_SOURCE_DIR}/gd_io.c
    ${LIBRARY_SOURCE_DIR}/gd_io_dp.c
    ${LIBRARY_SOURCE_DIR}/gd_io_file.c
    ${LIBRARY_SOURCE_DIR}/gd_io_ss.c
    ${LIBRARY_SOURCE_DIR}/gd_jpeg.c
    ${LIBRARY_SOURCE_DIR}/gd_nnquant.c
    ${LIBRARY_SOURCE_DIR}/gd_pixelate.c
    ${LIBRARY_SOURCE_DIR}/gd_png.c
    ${LIBRARY_SOURCE_DIR}/gd_security.c
    ${LIBRARY_SOURCE_DIR}/gd_ss.c
    ${LIBRARY_SOURCE_DIR}/gd_tga.c
    ${LIBRARY_SOURCE_DIR}/gd_tiff.c
    ${LIBRARY_SOURCE_DIR}/gd_topal.c
    ${LIBRARY_SOURCE_DIR}/gd_transform.c
    ${LIBRARY_SOURCE_DIR}/gd_wbmp.c
    ${LIBRARY_SOURCE_DIR}/gdfontg.c
    ${LIBRARY_SOURCE_DIR}/gdfontl.c
    ${LIBRARY_SOURCE_DIR}/gdfontmb.c
    ${LIBRARY_SOURCE_DIR}/gdfonts.c
    ${LIBRARY_SOURCE_DIR}/gdfontt.c
    ${LIBRARY_SOURCE_DIR}/gdft.c
    ${LIBRARY_SOURCE_DIR}/gdfx.c
    ${LIBRARY_SOURCE_DIR}/gdhelpers.c
    ${LIBRARY_SOURCE_DIR}/gdkanji.c
    ${LIBRARY_SOURCE_DIR}/gdtables.c
    ${LIBRARY_SOURCE_DIR}/gdxpm.c
    ${LIBRARY_SOURCE_DIR}/wbmp.c
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
