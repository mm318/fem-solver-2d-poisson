#ifndef _MESH_VIS_H_
#define _MESH_VIS_H_

#include "triangle/triangle.h"
#include "libpng_main/png_writer.h"


typedef struct {
    REAL x;
    REAL y;
} Coordinate;

typedef struct {
    Coordinate *start;
    Coordinate *end;
} Edge;

typedef struct {
    int num_points;
    Coordinate *points;
    int num_edges;
    Edge *edges;
    REAL x_min;
    REAL x_max;
    REAL y_min;
    REAL y_max;
    bitmap_t *bmp;
} MeshVisual;


MeshVisual* create_meshvisual(REAL * const points, const int num_points,
                              int * const edge_list, const int num_edges);

int draw_points(MeshVisual * const mesh, REAL * const points, const int num_points);

int draw_edges(MeshVisual * const mesh);

int delete_meshvisual(MeshVisual * const mesh);


#endif
