// Configured (by CMake) macros for PLplot that are required for the
// core build but _not_ required for the build of the installed
// examples (and presumably any user applications).  Therefore, the
// configured config.h should not be installed.  In contrast,
// include/plConfig.h.cmake (note, plConfig.h #includes config.h for
// the core build because HAVE_CONFIG_H is #defined in that case)
// contains configured macros that are required for the core build,
// installed examples build, and build of user applications.
// Therefore, in contrast to config.h, plConfig.h should be installed.
//
// Maintenance issue: in makes no sense to configure duplicate macros
// for both config.h and plConfig.h.  Therefore, when adding a macro
// decide which file to put it in depending on whether the result is
// needed for the installed examples build or not.  Furthermore, move
// configured macros from one file to the other as needed depending on
// that criterion, but do not copy them.
//
//

// Location of executables
#define BIN_DIR                  "./bin/"

// Location of Build tree
#define BUILD_DIR                "./build/"

// Location of package data files
#define DATA_DIR                 "./data/"

// Location of dynamically loaded drivers
#define DRV_DIR                  "./drv/"

// Name of the default cmap0 palette to use
#define PL_DEFAULT_CMAP0_FILE    "cmap0_white_bg.pal"

// Name of the default cmap1 palette to use
#define PL_DEFAULT_CMAP1_FILE    "cmap1_blue_red.pal"

// Define if support for deprecated plplot functions should be compiled
#define PL_DEPRECATED

// Define if there is support for dynamically loaded drivers
// #cmakedefine ENABLE_DYNDRIVERS

// Define to 1 if you have the <cmath> header file.
#define HAVE_CMATH 1

// Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
#undef HAVE_DIRENT_H

// Define to 1 if you have the <dlfcn.h> header file.
// #cmakedefine HAVE_DLFCN_H 1

// Define if finite is available
// #cmakedefine PL_HAVE_FINITE

// Define if _finite is available
// #cmakedefine PL__HAVE_FINITE

// Define if [freetype] is available
// #cmakedefine HAVE_FREETYPE

// Define if [agg] is available
// #cmakedefine HAVE_AGG

// Define to 1 if you have the <glib.h> header file.
// #cmakedefine HAVE_GLIB_H 1

// Define to 1 if you have the <glib-object.h> header file.
// #cmakedefine HAVE_GLIB_OBJECT_H 1

// Define to 1 if you have the <gtk/gtk.h> header file.
// #cmakedefine HAVE_GTK_GTK_H 1

// Define to 1 if you have the <inttypes.h> header file.
#define HAVE_INTTYPES_H 1

// Define if [incr], [Tcl] is available
// #cmakedefine HAVE_ITCL

// Define to 1 if you have the <itclDecls.h> header file.
// #cmakedefine HAVE_ITCLDECLS_H 1

// Define if Tk is available
// #cmakedefine ENABLE_tk

// Define if [incr], [Tk] is available
// #cmakedefine HAVE_ITK

// Define to 1 if you have the <jni.h> header file.
// #cmakedefine HAVE_JNI_H 1

// Define to 1 if you have the <libart_lgpl/libart.h> header file.
// #cmakedefine HAVE_LIBART_LGPL_LIBART_H 1

// Define to 1 if you have the <libgnomecanvas/libgnomecanvas.h> header file.
// #cmakedefine HAVE_LIBGNOMECANVAS_LIBGNOMECANVAS_H 1

// Define to 1 if you have the <libgnomeprint/gnome-print.h> header file.
// #cmakedefine HAVE_LIBGNOMEPRINT_GNOME_PRINT_H 1

// Define if libunicode is available
// #cmakedefine HAVE_LIBUNICODE

// Define to 1 if you have the <math.h> header file.
#define HAVE_MATH_H 1

// Define to 1 if you have the <memory.h> header file.
#define HAVE_MEMORY_H 1

// Define to 1 if the function mkstemp is available.
// #cmakedefine PL_HAVE_MKSTEMP 1

// Define to 1 if you have the <ndir.h> header file, and it defines `DIR'.
// #cmakedefine HAVE_NDIR_H 1

// Define if python numpy is available
// #cmakedefine HAVE_NUMPY

// Define if libpango is available
// #cmakedefine HAVE_PANGO

// Define if popen is available
// #cmakedefine HAVE_POPEN

// Define if _NSGetArgc is available
// #cmakedefine HAVE_NSGETARGC

// Define if pthreads is available
#define PL_HAVE_PTHREAD

// Define if Qhull is available
// #define PL_HAVE_QHULL

// Define to 1 if you have the <stdlib.h> header file.
#define HAVE_STDLIB_H 1

// Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
// #cmakedefine HAVE_SYS_DIR_H 1

// Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
// #cmakedefine HAVE_SYS_NDIR_H 1

// Define to 1 if you have the <sys/stat.h> header file.
#define HAVE_SYS_STAT_H 1

// Define to 1 if you have the <sys/types.h> header file.
#define HAVE_SYS_TYPES_H 1

// Define to 1 if you have <sys/wait.h> that is POSIX.1 compatible.
// #cmakedefine HAVE_SYS_WAIT_H 1

// Define to 1 if you have the <termios.h> header file.
// #cmakedefine HAVE_TERMIOS_H 1

// Define to 1 if you have the <crt_externs.h> header file.
// #cmakedefine HAVE_CRT_EXTERNS_H 1

// Define to 1 if the function unlink is available
// #cmakedefine PL_HAVE_UNLINK 1

// Define to 1 if you have the `vfork' function.
// #cmakedefine HAVE_VFORK 1

// Define to 1 if you have the <vfork.h> header file.
// #cmakedefine HAVE_VFORK_H 1

// Include sys/type.h if needed
// #cmakedefine NEED_SYS_TYPE_H

// Name of package
#define PACKAGE    "plplot"

// Define if the win32 ltdl implementation should be used
// #cmakedefine LTDL_WIN32

// Portable definition for PTHREAD_MUTEX_RECURSIVE
#define PLPLOT_MUTEX_RECURSIVE             PTHREAD_MUTEX_RECURSIVE

// Directory containing fonts that are accessible from freetype
#define PL_FREETYPE_FONT_DIR               "Times New Roman"

// MONO font accessible from freetype
#define PL_FREETYPE_MONO                   "Courier New"

// MONO_BOLD font accessible from freetype
#define PL_FREETYPE_MONO_BOLD              "Times New Roman"

// MONO_BOLD_ITALIC font accessible from freetype
#define PL_FREETYPE_MONO_BOLD_ITALIC       "Times New Roman"

// MONO_BOLD_OBLIQUE font accessible from freetype
#define PL_FREETYPE_MONO_BOLD_OBLIQUE      "Times New Roman"

// MONO_ITALIC font accessible from freetype
#define PL_FREETYPE_MONO_ITALIC            "Times New Roman"

// MONO_OBLIQUE font accessible from freetype
#define PL_FREETYPE_MONO_OBLIQUE           "Times New Roman"

// SANS font accessible from freetype
#define PL_FREETYPE_SANS                   "Times New Roman"

// SANS_BOLD font accessible from freetype
#define PL_FREETYPE_SANS_BOLD              "Times New Roman"

// SANS_BOLD_ITALIC font accessible from freetype
#define PL_FREETYPE_SANS_BOLD_ITALIC       "Times New Roman"

// SANS_BOLD_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SANS_BOLD_OBLIQUE      "Times New Roman"

// SANS_ITALIC font accessible from freetype
#define PL_FREETYPE_SANS_ITALIC            "Times New Roman"

// SANS_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SANS_OBLIQUE           "Times New Roman"

// SCRIPT font accessible from freetype
#define PL_FREETYPE_SCRIPT                 "Times New Roman"

// SCRIPT_BOLD font accessible from freetype
#define PL_FREETYPE_SCRIPT_BOLD            "Times New Roman"

// SCRIPT_BOLD_ITALIC font accessible from freetype
#define PL_FREETYPE_SCRIPT_BOLD_ITALIC     "Times New Roman"

// SCRIPT_BOLD_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SCRIPT_BOLD_OBLIQUE    "Times New Roman"

// SCRIPT_ITALIC font accessible from freetype
#define PL_FREETYPE_SCRIPT_ITALIC          "Times New Roman"

// SCRIPT_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SCRIPT_OBLIQUE         "Times New Roman"

// SERIF font accessible from freetype
#define PL_FREETYPE_SERIF                  "Times New Roman"

// SERIF_BOLD font accessible from freetype
#define PL_FREETYPE_SERIF_BOLD             "Times New Roman"

// SERIF_BOLD_ITALIC font accessible from freetype
#define PL_FREETYPE_SERIF_BOLD_ITALIC      "Times New Roman"

// SERIF_BOLD_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SERIF_BOLD_OBLIQUE     "Times New Roman"

// SERIF_ITALIC font accessible from freetype
#define PL_FREETYPE_SERIF_ITALIC           "Times New Roman"

// SERIF_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SERIF_OBLIQUE          "Times New Roman"

// Symbol font accessible from freetype
#define PL_FREETYPE_SYMBOL                 "Times New Roman"

// SYMBOL_BOLD font accessible from freetype
#define PL_FREETYPE_SYMBOL_BOLD            "Times New Roman"

// SYMBOL_BOLD_ITALIC font accessible from freetype
#define PL_FREETYPE_SYMBOL_BOLD_ITALIC     "Times New Roman"

// SYMBOL_BOLD_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SYMBOL_BOLD_OBLIQUE    "Times New Roman"

// SYMBOL_ITALIC font accessible from freetype
#define PL_FREETYPE_SYMBOL_ITALIC          "Times New Roman"

// SYMBOL_OBLIQUE font accessible from freetype
#define PL_FREETYPE_SYMBOL_OBLIQUE         "Times New Roman"

// Define as the return type of signal handlers (`int' or `void').
#define RETSIGTYPE                         int

// Location of Source tree
#define SOURCE_DIR                         "./src/"

// Define to 1 if you have the ANSI C header files.
#define STDC_HEADERS 1

// Location of Tcl stuff
#define TCL_DIR    ""

// Version number of package
#define VERSION    "5.9.9"

// Define if csa is desired
// #cmakedefine WITH_CSA

// Define if want to use general fill_intersection_polygon approach
// rather than the traditional code to fill the intersection of a polygon with
// the clipping limits.
#define USE_FILL_INTERSECTION_POLYGON

// Define to `char *' if <sys/types.h> does not define.
// #define caddr_t     char*

// Define to `int' if <sys/types.h> does not define.
// #define pid_t       int

// Define as `fork' if `vfork' does not work.
// #define vfork       fork

