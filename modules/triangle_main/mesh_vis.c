#include <stdio.h>
#include <stdlib.h>

#include "mesh_vis.h"


#ifndef NULL
#define NULL 0
#endif


const pixel_t white = {.red = 255, .green = 255, .blue = 255};
const pixel_t black = {.red = 0, .green = 0, .blue = 0};

#define MAX_IMAGE_DIMENSION 2000

MeshVisual* create_meshvisual(REAL * const points, const int num_points,
                              int * const edge_list, const int num_edges)
{
    int i, j;
    MeshVisual *mesh;
    REAL dx, dy;
    int height, width;

    mesh = (MeshVisual*) malloc(sizeof(MeshVisual));
    mesh->num_points = num_points;
    mesh->num_edges = num_edges;

    if(num_points > 0) {
        mesh->points = (Coordinate*) malloc(num_points*sizeof(Coordinate));
    } else {
        mesh->points = NULL;
    }
    if(num_edges > 0) {
        mesh->edges = (Edge*) malloc(num_edges*sizeof(Edge));
    } else {
        mesh->edges = NULL;
    }

    mesh->x_min = 9999999999;
    mesh->x_max = -9999999999;
    mesh->y_min = 9999999999;
    mesh->y_max = -9999999999;

    for(i=0; i<num_points; i++) {
        j = 2*i;
        if(points[j] < mesh->x_min)
            mesh->x_min = points[j];
        else if(points[j] > mesh->x_max)
            mesh->x_max = points[j];

        if(points[j+1] < mesh->y_min)
            mesh->y_min = points[j+1];
        else if(points[j+1] > mesh->y_max)
            mesh->y_max = points[j+1];

        mesh->points[i].x = points[j];
        mesh->points[i].y = points[j+1];
    }

    for(i=0; i<num_edges; i++) {
        j = 2*i;
        mesh->edges[i].start = &mesh->points[edge_list[j]];
        mesh->edges[i].end = &mesh->points[edge_list[j+1]];
    }

    dx = mesh->x_max - mesh->x_min;
    dy = mesh->y_max - mesh->y_min;
    if(dx > dy) {
        width = (int) MAX_IMAGE_DIMENSION;
        height = (int) MAX_IMAGE_DIMENSION*(dy/dx);
    } else {
        height = (int) MAX_IMAGE_DIMENSION;
        width = (int) MAX_IMAGE_DIMENSION*(dx/dy);
    }

    // creating new bitmap and setting background color to white
    mesh->bmp = new_bitmap(height, width);
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) {
            set_pixel(mesh->bmp, j, i, &white);
        }
    }

    return mesh;
}


#define RADIUS  3

int draw_points(MeshVisual * const mesh, REAL * const points, const int num_points)
{
    int i;
    int x, y;
    Coordinate *c;

    const REAL dx = mesh->x_max - mesh->x_min;
    const REAL dy = mesh->y_max - mesh->y_min;
    const REAL height = (REAL) mesh->bmp->height;
    const REAL width = (REAL) mesh->bmp->width;

    if(points == NULL) {
        // use stored coordinates
        for(i=0; i<mesh->num_points; i++) {
            c = &mesh->points[i];
            x = (int) ((c->x - mesh->x_min)/dx*width + 0.5);
            y = (int) ((c->y - mesh->y_min)/dy*height + 0.5);
            fill_circle(mesh->bmp, x, height-y, RADIUS, &black);

            // fprintf(stderr,"drew point at (%d, %d)\n", x, y);
        }
    } else {
        // use provided coordinates
        // for(i=0;i<num_points;i++) {
        // }
    }

    return 0;
}


#define LINE_THICKNESS  2

int draw_edges(MeshVisual * const mesh)
{
    int i;
    int x1, y1, x2, y2;
    Coordinate *s, *e;

    const REAL dx = mesh->x_max - mesh->x_min;
    const REAL dy = mesh->y_max - mesh->y_min;
    const REAL height = (REAL) mesh->bmp->height;
    const REAL width = (REAL) mesh->bmp->width;

    for(i=0; i<mesh->num_edges; i++) {
        s = mesh->edges[i].start;
        e = mesh->edges[i].end;

        x1 = (int) ((s->x - mesh->x_min)/dx*width + 0.5);
        y1 = (int) ((s->y - mesh->y_min)/dy*height + 0.5);
        x2 = (int) ((e->x - mesh->x_min)/dx*width + 0.5);
        y2 = (int) ((e->y - mesh->y_min)/dy*height + 0.5);

        draw_line(mesh->bmp, x1, height-y1, x2, height-y2, LINE_THICKNESS, &black);
    }

    return 0;
}

int delete_meshvisual(MeshVisual * const mesh)
{
    if(mesh->points != NULL) {
        free(mesh->points);
    }

    if(mesh->edges != NULL) {
        free(mesh->edges);
    }

    delete_bitmap(mesh->bmp);

    free(mesh);

    return 0;
}
