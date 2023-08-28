// $Id: x21c.c 11680 2011-03-27 17:57:51Z airwin $
//      Grid data demo
//
// Copyright (C) 2004  Joao Cardoso
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
//

#include <math.h>

#include "plplot/plplot.h"


#ifndef M_PI
#define M_PI    3.1415926535897932384
#endif

#ifndef MAX
#define MAX( a, b )    ( ( ( a ) > ( b ) ) ? ( a ) : ( b ) )
#endif

#ifndef MIN
#define MIN( a, b )    ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#endif


// Options data structure definition.
static PLINT    pts;
static PLINT    xp        = 100;
static PLINT    yp        = 100;
static PLINT    nl        = 32;
static int      knn_order = 20;
// static PLFLT    threshold = 1.001;
// static PLFLT    wmin      = -1e3;
// static int      randn     = 0;
// static int      rosen     = 0;

static PLFLT xm, xM, ym, yM;


void create_grid(PLFLT **xi, int px, PLFLT **yi, int py)
{
    PLFLT *x, *y;
    int   i;

    x = *xi = (PLFLT*) malloc( px * sizeof ( PLFLT ) );
    y = *yi = (PLFLT*) malloc( py * sizeof ( PLFLT ) );

    for ( i = 0; i < px; i++ )
        *x++ = xm + ( xM - xm ) * i / ( px - 1. );

    for ( i = 0; i < py; i++ )
        *y++ = ym + ( yM - ym ) * i / ( py - 1. );
}

void free_grid( PLFLT *xi, PLFLT *yi )
{
    free((void *) xi);
    free((void *) yi);
}

void get_data(double * const points_list, double * const u, const int num_points,
              PLFLT **xi, PLFLT **yi, PLFLT **zi, PLFLT *x_min, PLFLT *x_max,
              PLFLT *y_min, PLFLT *y_max, PLFLT *z_min, PLFLT *z_max)
{
    int i;
    PLFLT *x, *y, *z;
    PLFLT x_lmin = 999999, x_lmax = -999999;
    PLFLT y_lmin = 999999, y_lmax = -999999;
    PLFLT z_lmin = 999999, z_lmax = -999999;

    x = (PLFLT*) malloc(num_points*sizeof(PLFLT));
    y = (PLFLT*) malloc(num_points*sizeof(PLFLT));
    z = (PLFLT*) malloc(num_points*sizeof(PLFLT));

    pts = 0;
    for(i=0; i<num_points; i++) {
        if(!isnan(u[i])) {
            x[i] = (PLFLT) points_list[2*i];
            y[i] = (PLFLT) points_list[2*i+1];
            z[i] = (PLFLT) u[i];

            if(x[i] < x_lmin) {
                x_lmin = x[i];
            } else if(x[i] > x_lmax) {
                x_lmax = x[i];
            }

            if(y[i] < y_lmin) {
                y_lmin = y[i];
            } else if(y[i] > y_lmax) {
                y_lmax = y[i];
            }

            if(z[i] < z_lmin) {
                z_lmin = z[i];
            } else if(z[i] > z_lmax) {
                z_lmax = z[i];
            }

            pts++;
        }
    }

    *xi = x;
    *yi = y;
    *zi = z;
    *x_min = x_lmin;
    *x_max = x_lmax;
    *y_min = y_lmin;
    *y_max = y_lmax;
    *z_min = z_lmin;
    *z_max = z_lmax;
}

void free_data( PLFLT *x, PLFLT *y, PLFLT *z )
{
    free((void *) x);
    free((void *) y);
    free((void *) z);
}


int density_plot(double * const points_list, double * const u, const int num_points,
                 int * const line_list, const int num_lines, char * const name,
                 char * const filename)
{
    int i, j;
    int ii, jj;
    int argc = 4;
    const char *argv[] = {"-o", filename, "-geometry", "1200x900"};

    PLFLT *x, *y, *z, *clev;
    PLFLT *xg, *yg, **zg;
    PLFLT zmin, zmax, lzm, lzM;
    PLFLT dist, d;

    int p1, p2;
    PLFLT x_temp[2], y_temp[2];

    // plMergeOpts(options, "interpolate options", NULL);
    plparseopts(&argc, argv, PL_PARSE_NOPROGRAM|PL_PARSE_FULL);

    // Initialize plplot
    plinit();

    get_data(points_list, u, num_points, &x, &y, &z,
             &xm, &xM, &ym, &yM, &zmin, &zmax); // the sampled data
    create_grid(&xg, xp, &yg, yp);  // grid the data at
    plAlloc2dGrid(&zg, xp, yp);     // the output gridded data

    plgriddata(x, y, z, pts, xg, xp, yg, yp, zg, GRID_NNIDW, knn_order);

    for ( i = 0; i < xp; i++ ) {
        for ( j = 0; j < yp; j++ ) {
            if (isnan(zg[i][j])) {
                // average (IDW) over the 8 neighbors
                zg[i][j] = 0.;
                dist = 0.;

                for ( ii = i - 1; ii <= i + 1 && ii < xp; ii++ ) {
                    for ( jj = j - 1; jj <= j + 1 && jj < yp; jj++ ) {
                        if ( ii >= 0 && jj >= 0 && !isnan( zg[ii][jj] ) ) {
                            d         = ( abs( ii - i ) + abs( jj - j ) ) == 1 ? 1. : 1.4142;
                            zg[i][j] += zg[ii][jj] / ( d * d );
                            dist     += d;
                        }
                    }
                }
                if ( dist != 0. )
                    zg[i][j] /= dist;
                else
                    zg[i][j] = zmin;
            }
        }
    }

    plMinMax2dGrid((const PLFLT **) zg, xp, yp, &lzM, &lzm);

    lzm = MIN(lzm, zmin);
    lzM = MAX(lzM, zmax);

    // Increase limits slightly to prevent spurious contours
    // due to rounding errors
    // lzm = lzm - 0.01;
    // lzM = lzM + 0.01;

    clev = (PLFLT*) malloc(nl*sizeof(PLFLT));
    for(i = 0; i < nl; i++) {
        clev[i] = lzm + ( lzM - lzm ) / ( nl - 1 ) * i;
    }

    plcol0(1);
    plenv0(xm, xM, ym, yM, 1, 0);
    plcol0(15);
    pllab("X", "Y", name);
    plshades((const PLFLT **) zg, xp, yp, NULL, xm, xM, ym, yM, clev, nl,
             1, 0, 1, plfill, 1, NULL, NULL);
    // plotting outline of the geometry
    plcol0(15);
    if(line_list != NULL && num_lines > 0) {
        for(i=0; i<num_lines; i++) {
            p1 = line_list[2*i];
            p2 = line_list[2*i+1];
            x_temp[0] = points_list[2*p1];
            y_temp[0] = points_list[2*p1+1];
            x_temp[1] = points_list[2*p2];
            y_temp[1] = points_list[2*p2+1];
            plline(2, x_temp, y_temp);
        }
    }

    plend();

    free_data(x, y, z);
    free_grid( xg, yg );
    free((void *) clev );
    plFree2dGrid( zg, xp, yp );

    return 0;
}

static void cmap1_init()
{
    PLFLT i[2], h[2], l[2], s[2];

    i[0] = 0.0; // left boundary
    i[1] = 1.0; // right boundary

    h[0] = 240; // blue -> green -> yellow ->
    h[1] = 0;   // -> red

    l[0] = 0.6;
    l[1] = 0.6;

    s[0] = 0.8;
    s[1] = 0.8;

    plscmap1n( 256 );
    c_plscmap1l( 0, 2, i, h, l, s, NULL );
}

int mesh_plot(double * const points_list, double * const u, const int num_points,
              int * const line_list, const int num_lines, char * const name,
              char * const filename)
{
    int i, j;
    int ii, jj;
    int argc = 4;
    const char *argv[] = {"-o", filename, "-geometry", "1200x900"};

    PLFLT *x, *y, *z, *clev;
    PLFLT *xg, *yg, **zg;
    PLFLT zmin, zmax, lzm, lzM;
    PLFLT dist, d;

    int p1, p2;
    PLFLT x_temp[2], y_temp[2], z_temp[2];

    // plMergeOpts(options, "interpolate options", NULL);
    plparseopts(&argc, argv, PL_PARSE_NOPROGRAM|PL_PARSE_FULL);

    // Initialize plplot
    plinit();

    get_data(points_list, u, num_points, &x, &y, &z,
             &xm, &xM, &ym, &yM, &zmin, &zmax); // the sampled data
    create_grid(&xg, xp, &yg, yp);  // grid the data at
    plAlloc2dGrid(&zg, xp, yp);     // the output gridded data

    plgriddata( x, y, z, pts, xg, xp, yg, yp, zg, GRID_NNIDW, knn_order);

    for ( i = 0; i < xp; i++ ) {
        for ( j = 0; j < yp; j++ ) {
            if (isnan(zg[i][j])) {
                // average (IDW) over the 8 neighbors
                zg[i][j] = 0.;
                dist = 0.;

                for ( ii = i - 1; ii <= i + 1 && ii < xp; ii++ ) {
                    for ( jj = j - 1; jj <= j + 1 && jj < yp; jj++ ) {
                        if ( ii >= 0 && jj >= 0 && !isnan( zg[ii][jj] ) ) {
                            d         = ( abs( ii - i ) + abs( jj - j ) ) == 1 ? 1. : 1.4142;
                            zg[i][j] += zg[ii][jj] / ( d * d );
                            dist     += d;
                        }
                    }
                }
                if ( dist != 0. )
                    zg[i][j] /= dist;
                else
                    zg[i][j] = zmin;
            }
        }
    }

    plMinMax2dGrid( (const PLFLT **) zg, xp, yp, &lzM, &lzm );

    lzm = MIN(lzm, zmin);
    lzM = MAX(lzM, zmax);

    // Increase limits slightly to prevent spurious contours
    // due to rounding errors
    // lzm = lzm - 0.01;
    // lzM = lzM + 0.01;

    clev = (PLFLT*) malloc(nl*sizeof(PLFLT));
    for(i = 0; i < nl; i++) {
        clev[i] = lzm + ( lzM - lzm ) / ( nl - 1 ) * i;
    }
    cmap1_init();

    plcol0(1);
    pladv(0);
    plvpor(0.0, 1.0, 0.0, 0.9);
    plwind(-1.1, 0.75, -0.65, 1.20);
    //
    // For the comparison to be fair, all plots should have the
    // same z values, but to get the max/min of the data generated
    // by all algorithms would imply two passes. Keep it simple.
    //
    // plw3d(1., 1., 1., xm, xM, ym, yM, zmin, zmax, 30, -60);
    //
    plw3d(1., 1., 1., xm, xM, ym, yM, lzm, lzM, 30, -40);
    plbox3("bntu", "X", 0., 0, "bntu", "Y", 0., 0, "bcdfntu", "", 0., 0);
    plcol0(15);
    pllab("", "", name);
    // plotting outline of the geometry
    if(line_list != NULL && num_lines > 0) {
        z_temp[0] = lzm;
        z_temp[1] = lzm;
        for(i=0; i<num_lines; i++) {
            p1 = line_list[2*i];
            p2 = line_list[2*i+1];
            x_temp[0] = points_list[2*p1];
            y_temp[0] = points_list[2*p1+1];
            x_temp[1] = points_list[2*p2];
            y_temp[1] = points_list[2*p2+1];
            plline3(2, x_temp, y_temp, z_temp);
        }
    }
    plmeshc(xg, yg, (const PLFLT **) zg, xp, yp, DRAW_LINEXY|MAG_COLOR|BASE_CONT, clev, nl);

    plend();

    free_data(x, y, z);
    free_grid( xg, yg );
    free((void *) clev);
    plFree2dGrid( zg, xp, yp );

    return 0;
}
