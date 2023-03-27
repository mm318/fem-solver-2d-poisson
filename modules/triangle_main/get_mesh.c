#ifndef NULL
#define NULL 0
#endif

#include <stdio.h>
#include <stdlib.h>

#include "node_input.h"
#include "mesh_vis.h"


static double target_area(TriangulateIO *ti)
{
    int i, j;
    double x_min = 9999999999;
    double x_max = -9999999999;
    double y_min = 9999999999;
    double y_max = -9999999999;
    double target;

    for(i=0; i<ti->numberofpoints; i++) {
        j = 2*i;
        if(ti->pointlist[j] < x_min) {
            x_min = ti->pointlist[j];
        } else if(ti->pointlist[j] > x_max) {
            x_max = ti->pointlist[j];
        }

        if(ti->pointlist[j+1] < y_min) {
            y_min = ti->pointlist[j+1];
        } else if(ti->pointlist[j+1] > y_max) {
            y_max = ti->pointlist[j+1];
        }
    }

    target = (x_max - x_min)*(y_max - y_min)/200.0;
    // target = 1000.0;

    return target;
}


static void output_mesh(TriangulateIO *ti, char * const name)
{
    int i;
    FILE *out;
    char temp_buffer[1024];
    double x, y;
    int v1, v2, v3;

    sprintf(temp_buffer, "%s_nodes.txt", name);
    out = fopen(temp_buffer, "wb");
    for(i=0; i<ti->numberofpoints; i++) {
        x = ti->pointlist[2*i];
        y = ti->pointlist[2*i+1];
        fprintf(out, "    %15.7f    %15.7f\n", x, y);
    }
    fclose(out);

    sprintf(temp_buffer, "%s_elements.txt", name);
    out = fopen(temp_buffer, "wb");
    for(i=0; i<ti->numberoftriangles; i++) {
        v1 = ti->trianglelist[3*i] + 1;
        v2 = ti->trianglelist[3*i+1] + 1;
        v3 = ti->trianglelist[3*i+2] + 1;
        fprintf(out, "    %d    %d    %d\n", v1, v2, v3);
    }
    fclose(out);
}


TriangulateIO* get_mesh(char * const name, char * const filename)
{
    // int i;
    char temp_buffer[1024];
    double area;

    TriangulateIO *in = new_TriangulateIO();
    TriangulateIO *out = new_TriangulateIO();

    if(read_polyfile(in, filename) != 0) {
        printf("Unable to read \"%s\"\n", filename);
        return NULL;
    }
    area = target_area(in);
    sprintf(temp_buffer, "pezqa%fV", area);
    triangulate(temp_buffer, in, out, NULL);

    // debug output
    {
        printf("Number of points: %d\n", out->numberofpoints);
        // for(i=0; i<out->numberofpoints; i++) {
        //     printf("Point %d: (x = %f, y = %f)\n", i+1,
        //            out->pointlist[2*i], out->pointlist[2*i+1]);
        // }
        printf("\nNumber of segments: %d\n", out->numberofsegments);
        printf("\nNumber of edges: %d\n", out->numberofedges);
        // for(i=0; i<out->numberofedges; i++) {
        //     printf("Edge %d: (Point 1 = %d, Point 2 = %d)\n", i+1,
        //            out->edgelist[2*i], out->edgelist[2*i+1]);
        // }
        printf("\nNumber of triangles: %d\n", out->numberoftriangles);
        printf("\nNumber of corners: %d\n", out->numberofcorners);

        /* displaying triangulation */
        puts("\nWriting out images...");

        // MeshVisual *initial = create_meshvisual(in->pointlist, in->numberofpoints,
        //                                         in->segmentlist, in->numberofsegments);
        // draw_points(initial, NULL, 0);
        // sprintf(temp_buffer, "%s_initial.png", name);
        // write_bitmap(initial->bmp, temp_buffer); /* Write the image to a file */
        // draw_edges(initial);
        // sprintf(temp_buffer, "%s_initial_edges.png", name);
        // write_bitmap(initial->bmp, temp_buffer);  /* Write the image to a file */

        MeshVisual *final = create_meshvisual(out->pointlist, out->numberofpoints,
                                              out->edgelist, out->numberofedges);
        draw_points(final, NULL, 0);
        sprintf(temp_buffer, "%s_mesh_points.png", name);
        write_bitmap(final->bmp, temp_buffer);  /* Write the image to a file */
        draw_edges(final);
        sprintf(temp_buffer, "%s_mesh_edges.png", name);
        write_bitmap(final->bmp, temp_buffer);  /* Write the image to a file */

        // delete_meshvisual(initial);
        delete_meshvisual(final);

        output_mesh(out, name);
    }
    // end of debug

    // cleanup
    fprintf(stderr, "\nFreeing up resources\n\n");
    delete_TriangulateIO(in);

    return out;
}
