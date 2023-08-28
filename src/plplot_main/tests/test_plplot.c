// $Id: plcdemos.h 11289 2010-10-29 20:44:17Z airwin $
//
//      Everything needed by the C demo programs.
//      Created to avoid junking up plplot.h with this stuff.
//

#include <math.h>
#include <string.h>
#include <ctype.h>

// #include "lib/plplot/include/plConfig.h"
#include "plplot/plplot.h"

// define PI if not defined by math.h

// Actually M_PI seems to be more widely used so we deprecate PI.
#ifndef PI
#define PI    3.1415926535897932384
#endif

#ifndef M_PI
#define M_PI    3.1415926535897932384
#endif

// various utility macros

#ifndef MAX
#define MAX( a, b )    ( ( ( a ) > ( b ) ) ? ( a ) : ( b ) )
#endif

#ifndef MIN
#define MIN( a, b )    ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#endif

#ifndef ROUND
#define ROUND( a )    (PLINT) ( ( a ) < 0. ? ( ( a ) - .5 ) : ( ( a ) + .5 ) )
#endif

// Declarations for save string functions

#ifdef PL_HAVE_SNPRINTF
// In case only _snprintf is declared (as for Visual C++ and
// Borland compiler toolset) we redefine the function names
#ifdef _PL_HAVE_SNPRINTF
#define snprintf    _snprintf
#define snscanf     _snscanf
#endif // _PL_HAVE_SNPRINTF
#else // !PL_HAVE_SNPRINTF
// declare dummy functions which just call the unsafe
// functions ignoring the size of the string
int plsnprintf( char *buffer, int n, const char *format, ... );
int plsnscanf( const char *buffer, int n, const char *format, ... );
#define snprintf    plsnprintf
#define snscanf     plsnscanf
#endif // PL_HAVE_SNPRINTF

// Add in missing isnan definition if required
#if defined ( PL__HAVE_ISNAN )
#  define isnan    _isnan
#  if defined ( _MSC_VER )
#    include <float.h>
#  endif
#endif

#if !defined ( PL_HAVE_ISNAN )
#  define isnan( x )    ( ( x ) != ( x ) )
#endif


// #include "plcdemos.h"

// $Id: x11c.c 11680 2011-03-27 17:57:51Z airwin $
//
//      Mesh plot demo.
//
// Copyright (C) 2004  Rafael Laboissiere
//
// This file is part of PLplot.
//
// PLplot is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PLplot is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with PLplot; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//


#define XPTS      35            // Data points in x
#define YPTS      46            // Data points in y
#define LEVELS    10

static int   opt[] = { DRAW_LINEXY, DRAW_LINEXY };

static PLFLT alt[] = { 33.0, 17.0 };
static PLFLT az[] = { 24.0, 115.0 };

static char  *title[4] = {
    "#frPLplot Example 11 - Alt=33, Az=24, Opt=3",
    "#frPLplot Example 11 - Alt=17, Az=115, Opt=3",
};

static void cmap1_init()
{
    PLFLT i[2], h[2], l[2], s[2];

    i[0] = 0.0;         // left boundary
    i[1] = 1.0;         // right boundary

    h[0] = 240;         // blue -> green -> yellow ->
    h[1] = 0;           // -> red

    l[0] = 0.6;
    l[1] = 0.6;

    s[0] = 0.8;
    s[1] = 0.8;

    plscmap1n( 256 );
    c_plscmap1l( 0, 2, i, h, l, s, NULL );
}

//--------------------------------------------------------------------------
// main
//
// Does a series of mesh plots for a given data set, with different
// viewing options in each plot.
//--------------------------------------------------------------------------

int main()
{
    int   i, j, k;
    PLFLT *x, *y, **z;
    PLFLT xx, yy;
    int   nlevel = LEVELS;
    PLFLT clevel[LEVELS];
    PLFLT zmin, zmax, step;

    // Parse and process command line arguments
    int argc = 4;
    const char *argv[] = {"-o", "plplot_test.png", "-geometry", "1200x900"};
    (void) plparseopts(&argc, argv, PL_PARSE_NOPROGRAM|PL_PARSE_FULL);

    // Initialize plplot
    plinit();

    x = (PLFLT *) calloc( XPTS, sizeof ( PLFLT ) );
    y = (PLFLT *) calloc( YPTS, sizeof ( PLFLT ) );

    plAlloc2dGrid( &z, XPTS, YPTS );
    for ( i = 0; i < XPTS; i++ ) {
        x[i] = 3. * (double) ( i - ( XPTS / 2 ) ) / (double) ( XPTS / 2 );
    }

    for ( i = 0; i < YPTS; i++ )
        y[i] = 3. * (double) ( i - ( YPTS / 2 ) ) / (double) ( YPTS / 2 );

    for ( i = 0; i < XPTS; i++ ) {
        xx = x[i];
        for ( j = 0; j < YPTS; j++ ) {
            yy      = y[j];
            z[i][j] = 3. * ( 1. - xx ) * ( 1. - xx ) * exp( -( xx * xx ) - ( yy + 1. ) * ( yy + 1. ) ) -
                      10. * ( xx / 5. - pow( xx, 3. ) - pow( yy, 5. ) ) * exp( -xx * xx - yy * yy ) -
                      1. / 3. * exp( -( xx + 1 ) * ( xx + 1 ) - ( yy * yy ) );

            if ( 0 ) { // Jungfraujoch/Interlaken
                if ( z[i][j] < -1. )
                    z[i][j] = -1.;
            }
        }
    }

    plMinMax2dGrid( (const PLFLT **) z, XPTS, YPTS, &zmax, &zmin );
    step = ( zmax - zmin ) / ( nlevel + 1 );
    for ( i = 0; i < nlevel; i++ )
        clevel[i] = zmin + step + step * i;

    cmap1_init();
    for ( k = 0; k < 2; k++ ) {
        for ( i = 0; i < 4; i++ ) {
            pladv( 0 );
            plcol0( 1 );
            plvpor( 0.0, 1.0, 0.0, 0.9 );
            plwind( -1.0, 1.0, -1.0, 1.5 );
            plw3d( 1.0, 1.0, 1.2, -3.0, 3.0, -3.0, 3.0, zmin, zmax, alt[k], az[k] );
            plbox3( "bnstu", "x axis", 0.0, 0,
                    "bnstu", "y axis", 0.0, 0,
                    "bcdmnstuv", "z axis", 0.0, 4 );

            plcol0( 2 );

#if 0
            // wireframe plot
            if ( i == 0 )
                plmesh( x, y, (const PLFLT **) z, XPTS, YPTS, opt[k] );

            // magnitude colored wireframe plot
            else if ( i == 1 )
                plmesh( x, y, (const PLFLT **) z, XPTS, YPTS, opt[k] | MAG_COLOR );

            // magnitude colored wireframe plot with sides
            else if ( i == 2 )
                plot3d( x, y, (const PLFLT **) z, XPTS, YPTS, opt[k] | MAG_COLOR, 1 );

            // magnitude colored wireframe plot with base contour
            else if ( i == 3 )
                plmeshc( x, y, (const PLFLT **) z, XPTS, YPTS, opt[k] | MAG_COLOR | BASE_CONT,
                         clevel, nlevel );
#endif

            plmeshc( x, y, (const PLFLT **) z, XPTS, YPTS, opt[k] | MAG_COLOR | BASE_CONT, clevel, nlevel );

            plcol0( 3 );
            plmtex( "t", 1.0, 0.5, 0.5, title[k] );
        }
    }

    // Clean up
    free( (void *) x );
    free( (void *) y );
    plFree2dGrid( z, XPTS, YPTS );
    plend();

    return 0;
}

