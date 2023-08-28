cmake_minimum_required(VERSION 3.0.0)

include(GNUInstallDirs)

set(PROJECT_NAME plplot_Library)
project(${PROJECT_NAME})

set(LIBRARY_NAME plplot)

set(LIBRARY_HEADERS
    csa/nan.h
    fonts/plhershey-unicode.h
    include/config.h
    include/dirent_msvc.h
    include/disptab.h
    include/drivers.h
    include/gcw.h
    include/ltdl_win32.h
    include/metadefs.h
    include/pdf.h
    include/plConfig.h
    include/plcore.h
    include/pldebug.h
    include/plDevs.h
    include/pldll.h
    include/plevent.h
    include/plfci-truetype.h
    include/plfci-type1.h
    include/plfreetype.h
    include/plplot.h
    include/plplotcanvas-hacktext.h
    include/plplotcanvas.h
    include/plplotP.h
    include/plstrm.h
    include/plunicode-type1.h
    include/plxwd.h
    include/ps.h
    include/qt.h
    qsastime/qsastime.h
    qsastime/qsastimedll.h
    qsastime/qsastimeP.h
    qsastime/tai-utc.h
    src/mt19937ar.h
)

set(LIBRARY_SOURCES
    drivers/gd.c
    qsastime/qsastime.c
    src/mt19937ar.c
    src/pdfutils.c
    src/plaffine.c
    src/plarc.c
    src/plargs.c
    src/plbox.c
    src/plbuf.c
    src/plcont.c
    src/plcore.c
    src/plctrl.c
    src/plcvt.c
    src/pldeprecated.c
    src/pldtik.c
    src/plf2ops.c
    src/plfill.c
    src/plfreetype.c
    src/plgradient.c
    src/plgridd.c
    src/plhist.c
    src/plimage.c
    src/pllegend.c
    src/plline.c
    src/plmap.c
    src/plot3d.c
    src/plpage.c
    src/plsdef.c
    src/plshade.c
    src/plstdio.c
    src/plstripc.c
    src/plsym.c
    src/pltick.c
    src/pltime.c
    src/plvect.c
    src/plvpor.c
    src/plwind.c
)

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_SOURCES}
)

target_include_directories(${LIBRARY_NAME} BEFORE
    PRIVATE "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}"
)

set_target_properties(${LIBRARY_NAME} PROPERTIES
    PUBLIC_HEADER "${LIBRARY_HEADERS}"
)

install(TARGETS ${LIBRARY_NAME}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${LIBRARY_NAME}"
)
