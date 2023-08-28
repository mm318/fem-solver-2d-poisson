#ifndef _NODE_INPUT_H_
#define _NODE_INPUT_H_

#include "triangle/triangle.h"


TriangulateIO* new_TriangulateIO();

/* Read the vertices from a file, which must be a .node */
// int read_nodefile(TriangulateIO * const in, char * const filename);

/* Read the vertices from a file, which must be a .poly */
int read_polyfile(TriangulateIO * const in, char * const filename);

int delete_TriangulateIO(TriangulateIO* const ti);


#endif

