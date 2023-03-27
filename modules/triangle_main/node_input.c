#include <stdio.h>
#include <stdlib.h>

#include "node_input.h"


TriangulateIO* new_TriangulateIO()
{
    TriangulateIO *ti = (TriangulateIO*) malloc(sizeof(TriangulateIO));

    ti->pointlist = NULL;
    ti->pointattributelist = NULL;
    ti->pointmarkerlist = NULL;
    ti->numberofpoints = 0;
    ti->numberofpointattributes = 0;
    ti->trianglelist = NULL;
    ti->triangleattributelist = NULL;
    ti->trianglearealist = NULL;
    ti->neighborlist = NULL;
    ti->numberoftriangles = 0;
    ti->numberofcorners = 0;
    ti->numberoftriangleattributes = 0;
    ti->segmentlist = NULL;
    ti->segmentmarkerlist = NULL;
    ti->numberofsegments = 0;
    ti->holelist = NULL;
    ti->numberofholes = 0;
    ti->regionlist = NULL;
    ti->numberofregions = 0;
    ti->edgelist = NULL;
    ti->edgemarkerlist = NULL;
    ti->normlist = NULL;
    ti->numberofedges = 0;

    return ti;
}


#define INPUTLINESIZE   1024

/*  Read a nonempty line from a file.                                           */
/*  A line is considered "nonempty" if it contains something that looks like    */
/*  a number.  Comments (prefaced by `#') are ignored.                          */
char *readline(char *string, FILE *infile, char *infilename)
{
    char *result;

    /* Search for something that looks like a number. */
    do {
        result = fgets(string, INPUTLINESIZE, infile);
        // fprintf(stderr, "%s\n", result);
        if (result == (char *) NULL) {
            printf("  Error:  Unexpected end of file in %s.\n", infilename);
            return 0;
        }

        /* Skip anything that doesn't look like a number, a comment, */
        /*   or the end of a line.                                   */
        while ((*result != '\0') && (*result != '#')
               && (*result != '.') && (*result != '+') && (*result != '-')
               && ((*result < '0') || (*result > '9'))) {
            result++;
        }
        /* If it's a comment or end of line, read another line and try again. */
    } while ((*result == '#') || (*result == '\0'));
    return result;
}

/*  Find the next field of a string.                                        */
/*  Jumps past the current field by searching for whitespace, then jumps    */
/*  past the whitespace to find the next field.                             */
char *findfield(char *string)
{
    char *result;

    result = string;
    /* Skip the current field.  Stop upon reaching whitespace. */
    while ((*result != '\0') && (*result != '#')
           && (*result != ' ') && (*result != '\t')) {
        result++;
    }
    /* Now skip the whitespace and anything else that doesn't look like a */
    /*   number, a comment, or the end of a line.                         */
    while ((*result != '\0') && (*result != '#')
           && (*result != '.') && (*result != '+') && (*result != '-')
           && ((*result < '0') || (*result > '9'))) {
        result++;
    }
    /* Check for a comment (prefixed with `#'). */
    if (*result == '#') {
        *result = '\0';
    }
    return result;
}

/* Read the vertices from a file, which may be a .node or .poly file */
static int read_nodes(TriangulateIO * const ti, int *index_offset,
                      FILE *infile, char * const filename)
{
    int i, j;

    char inputline[INPUTLINESIZE];
    char *stringptr;
    int mesh_dim;
    int bound_markers;
    REAL x, y;

    /* Read number of vertices, number of dimensions, number of vertex */
    /*   attributes, and number of boundary markers.                   */
    stringptr = readline(inputline, infile, filename);
    ti->numberofpoints = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
        mesh_dim = 2;
    } else {
        mesh_dim = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
        ti->numberofpointattributes = 0;
    } else {
        ti->numberofpointattributes = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
        bound_markers = 0;
    } else {
        bound_markers = (int) strtol(stringptr, &stringptr, 0);
    }

    if (ti->numberofpoints < 3) {
        printf("Error: Input must have at least three input vertices.\n");
        return -1;
    }
    if (mesh_dim != 2) {
        printf("Error: Triangle only works with two-dimensional meshes.\n");
        return -1;
    }

    // allocating memory to store vertex data
    ti->pointlist = (REAL*) malloc(2 * ti->numberofpoints * sizeof(REAL));
    if(ti->numberofpointattributes > 0) {
        ti->pointattributelist = (REAL*) malloc(ti->numberofpointattributes * ti->numberofpoints * sizeof(REAL));
    } else {
        ti->pointattributelist = NULL;
    }
    if(bound_markers) {
        ti->pointmarkerlist = (int*) malloc(ti->numberofpoints * sizeof(int));
    } else {
        ti->pointmarkerlist = NULL;
    }

    /* Read the vertices. */
    for (i = 0; i < ti->numberofpoints; i++) {
        stringptr = readline(inputline, infile, filename);
        if (i == 0) {
            *index_offset = (int) strtol(stringptr, &stringptr, 0);
        }

        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
            printf("Error: Vertex %d has no x coordinate.\n", *index_offset + i);
            return -1;
        }
        x = (REAL) strtod(stringptr, &stringptr);
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
            printf("Error: Vertex %d has no y coordinate.\n", *index_offset + i);
            return -1;
        }
        y = (REAL) strtod(stringptr, &stringptr);

        ti->pointlist[2*i] = x;
        ti->pointlist[2*i + 1] = y;
        // fprintf(stderr, "Point #%d: (x: %f, y: %f)\n", i, x, y);

        /* Read the vertex attributes. */
        for (j = 0; j < ti->numberofpointattributes; j++) {
            stringptr = findfield(stringptr);
            if (*stringptr == '\0') {
                ti->pointlist[ti->numberofpointattributes*i + j] = 0.0;
            } else {
                ti->pointlist[ti->numberofpointattributes*i + j] = (REAL) strtod(stringptr, &stringptr);
            }
        }

        /* Read a vertex marker. */
        if (bound_markers) {
            stringptr = findfield(stringptr);
            if (*stringptr == '\0') {
                ti->pointmarkerlist[i] = 0;
            } else {
                ti->pointmarkerlist[i] = (int) strtol(stringptr, &stringptr, 0);
                /* ti->pointmarkerlist[i] = 1; */
            }
        }
    }

    return 0;
}


// Read the lines from a file, which must be a .poly file
static int read_edges(TriangulateIO * const ti, int *index_offset,
                      FILE *infile, char * const filename)
{
    int i;

    char inputline[INPUTLINESIZE];
    char *stringptr;
    int bound_markers;
    int first_segment_num;
    int index1, index2;

    // Reading number of line edges, number of boundary markers.
    stringptr = readline(inputline, infile, filename);
    ti->numberofsegments = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
        bound_markers = 0;
    } else {
        bound_markers = (int) strtol(stringptr, &stringptr, 0);
    }

    // allocating memory to store line segment data
    ti->segmentlist = (int*) malloc(2 * ti->numberofsegments * sizeof(int));
    if(bound_markers) {
        ti->segmentmarkerlist = (int*) malloc(ti->numberofsegments * sizeof(int));
    } else {
        ti->segmentmarkerlist = NULL;
    }

    // Read the line segments
    for (i = 0; i < ti->numberofsegments; i++) {
        stringptr = readline(inputline, infile, filename);
        if (i == 0) {
            first_segment_num = (int) strtol(stringptr, &stringptr, 0);
        }

        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
            printf("Error: Segment %d is missing a vertex index.\n", first_segment_num + i);
            return -1;
        }
        index1 = (int) strtol(stringptr, &stringptr, 0);
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
            printf("Error: Segment %d is missing a vertex index.\n", first_segment_num + i);
            return -1;
        }
        index2 = (int) strtol(stringptr, &stringptr, 0);

        ti->segmentlist[2*i] = index1 - *index_offset;
        ti->segmentlist[2*i + 1] = index2 - *index_offset;

        // Reading line segment marker
        if (bound_markers) {
            stringptr = findfield(stringptr);
            if (*stringptr == '\0') {
                ti->segmentmarkerlist[i] = 0;
            } else {
                ti->segmentmarkerlist[i] = (int) strtol(stringptr, &stringptr, 0);
                /* ti->pointmarkerlist[i] = 1; */
            }
        }
    }

    return 0;
}

static int read_holes(TriangulateIO * const ti, int *index_offset,
                      FILE *infile, char * const filename)
{
    int i;
    char inputline[INPUTLINESIZE];
    char *stringptr;

    /* Read the holes. */
    stringptr = readline(inputline, infile, filename);
    ti->numberofholes = (int) strtol(stringptr, &stringptr, 0);
    if (ti->numberofholes > 0) {
        ti->holelist = (REAL*) trimalloc(2*ti->numberofholes*sizeof(REAL));
        for (i = 0; i < ti->numberofholes; i++) {
            stringptr = readline(inputline, infile, filename);

            stringptr = findfield(stringptr);
            if (*stringptr == '\0') {
                printf("Error:  Hole %d has no x coordinate.\n", i+1);
                return -1;
            } else {
                ti->holelist[2*i] = (REAL) strtod(stringptr, &stringptr);
            }

            stringptr = findfield(stringptr);
            if (*stringptr == '\0') {
                printf("Error:  Hole %d has no y coordinate.\n", i+1);
                return -1;
            } else {
                ti->holelist[2*i + 1] = (REAL) strtod(stringptr, &stringptr);
            }
        }
    } else {
        ti->holelist = (REAL*) NULL;
    }

    return 0;
}


// getting input from a .poly file
int read_polyfile(TriangulateIO * const ti, char * const filename)
{
    FILE *infile;
    int index_offset;

    printf("Opening %s.\n", filename);
    infile = fopen(filename, "r");
    if (infile == (FILE *) NULL) {
        printf("  Error:  Cannot access file %s.\n", filename);
        return -1;
    }

    read_nodes(ti, &index_offset, infile, filename);
    read_edges(ti, &index_offset, infile, filename);
    read_holes(ti, &index_offset, infile, filename);

    fclose(infile);

    return 0;
}


int delete_TriangulateIO(TriangulateIO* const ti)
{
    if(ti->pointlist != NULL) {
        trifree(ti->pointlist);
    }

    if(ti->pointattributelist != NULL) {
        trifree(ti->pointattributelist);
    }

    if(ti->pointmarkerlist != NULL) {
        trifree(ti->pointmarkerlist);
    }

    if(ti->trianglelist != NULL) {
        trifree(ti->trianglelist);
    }

    if(ti->triangleattributelist != NULL) {
        trifree(ti->triangleattributelist);
    }

    if(ti->trianglearealist != NULL) {
        trifree(ti->trianglearealist);
    }

    if(ti->neighborlist != NULL) {
        trifree(ti->neighborlist);
    }

    if(ti->segmentlist != NULL) {
        trifree(ti->segmentlist);
    }

    if(ti->segmentmarkerlist != NULL) {
        trifree(ti->segmentmarkerlist);
    }

    if(ti->holelist != NULL) {
        trifree(ti->holelist);
    }

    if(ti->regionlist != NULL) {
        trifree(ti->regionlist);
    }

    if(ti->edgelist != NULL) {
        trifree(ti->edgelist);
    }

    if(ti->edgemarkerlist != NULL) {
        trifree(ti->edgemarkerlist);
    }

    if(ti->normlist != NULL) {
        trifree(ti->normlist);
    }

    trifree(ti);

    return 0;
}
