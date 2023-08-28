// -*-C-*-
// $Id: plDevs.h.cmake 11269 2010-10-22 00:25:36Z airwin $
//
//  Maurice LeBrun
//  IFS, University of Texas at Austin
//  18-Jul-1994
//
//  Contains macro definitions that determine what device drivers are
//  compiled into the PLplot library.  On a Unix system, typically the
//  configure script builds plDevs.h from plDevs.h.in.  Elsewhere, it's
//  best to hand-configure a plDevs.h file and keep it with the
//  system-specific files.
//
//  Copyright (C) 2004  Andrew Roach
//  Copyright (C) 2005  Thomas J. Duck
//  Copyright (C) 2006  Andrew Ross
//  Copyright (C) 2006  Alan W. Irwin
//
//  This file is part of PLplot.
//
//  PLplot is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Library General Public License as published
//  by the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  PLplot is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Library General Public License for more details.
//
//  You should have received a copy of the GNU Library General Public License
//  along with PLplot; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
//

#ifndef __PLDEVS_H__
#define __PLDEVS_H__

// #cmakedefine PDL_aqt
// #cmakedefine PLD_plmeta
// #cmakedefine PLD_null
// #cmakedefine PLD_xterm
// #cmakedefine PLD_tek4010
// #cmakedefine PLD_tek4010f
// #cmakedefine PLD_tek4107
// #cmakedefine PLD_tek4107f
// #cmakedefine PLD_mskermit
// #cmakedefine PLD_vlt
// #cmakedefine PLD_versaterm
// #cmakedefine PLD_conex
// #cmakedefine PLD_linuxvga
// #cmakedefine PLD_dg300
#define PLD_png
// #define PLD_jpeg
// #define PLD_gif
// #cmakedefine PLD_cgm
// #cmakedefine PLD_ps
// #cmakedefine PLD_xfig
// #cmakedefine PLD_ljiip
// #cmakedefine PLD_ljii
// #cmakedefine PLD_lj_hpgl
// #cmakedefine PLD_hp7470
// #cmakedefine PLD_hp7580
// #cmakedefine PLD_imp
// #cmakedefine PLD_xwin
// #cmakedefine PLD_tk
// #cmakedefine PLD_pbm
// #cmakedefine PLD_gcw
// #cmakedefine PLD_gnome
// #cmakedefine PLD_pstex
// #cmakedefine PLD_psttf
// #cmakedefine PLD_ntk
// #cmakedefine PLD_tkwin
// #cmakedefine PLD_mem
// #define PLD_wingcc
// #cmakedefine PLD_wxwidgets
// #cmakedefine PLD_wxpng
// #cmakedefine PLD_svg
// #cmakedefine PLD_pdf
// #cmakedefine PLD_xcairo
// #cmakedefine PLD_pdfcairo
// #cmakedefine PLD_pscairo
// #cmakedefine PLD_svgcairo
// #cmakedefine PLD_pngcairo
// #cmakedefine PLD_memcairo
// #cmakedefine PLD_extcairo
// #cmakedefine PLD_wincairo
// #cmakedefine PLD_bmpqt
// #cmakedefine PLD_jpgqt
// #cmakedefine PLD_pngqt
// #cmakedefine PLD_ppmqt
// #cmakedefine PLD_tiffqt
// #cmakedefine PLD_svgqt
// #cmakedefine PLD_epsqt
// #cmakedefine PLD_pdfqt
// #cmakedefine PLD_qtwidget
// #cmakedefine PLD_extqt
// #cmakedefine PLD_memqt

#endif  // __PLDEVS_H__
