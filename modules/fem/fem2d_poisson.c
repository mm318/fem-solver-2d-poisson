#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "fem2d_poisson.h"
#include "dgb.h"

// for debug functions
#include "main/fem_solver.h"


//****************************************************************************80

char* triangulation_order3_boundary_node(int node_num, int element_num, int *element_node)

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_BOUNDARY_NODE indicates nodes on the boundary.
//
//  Discussion:
//
//    This routine is given a triangulation, an abstract list of triples
//    of nodes.  It is assumed that the nodes in each triangle are listed
//    in a counterclockwise order, although the routine should work
//    if the nodes are consistently listed in a clockwise order as well.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, char TRIANGULATION_ORDER3_BOUNDARY_NODE[NODE_NUM],
//    is TRUE if the node is on a boundary edge.
//
{
    int e1;
    int e2;
    int *edge;
    char equal;
    int i;
    int j;
    int m;
    int n;
    char *node_boundary;

    //
    //  Set up the edge array.
    //
    m = 2;
    n = 3 * element_num;
    edge = (int*) malloc(m*n*sizeof(int));
    for(j = 0; j < element_num; j++) {
        edge[0+(j              )*m] = element_node[0+j*3];
        edge[1+(j              )*m] = element_node[1+j*3];
        edge[0+(j+  element_num)*m] = element_node[1+j*3];
        edge[1+(j+  element_num)*m] = element_node[2+j*3];
        edge[0+(j+2*element_num)*m] = element_node[2+j*3];
        edge[1+(j+2*element_num)*m] = element_node[0+j*3];
    }

    //
    //  In each column, force the smaller entry to appear first.
    //
    for (j = 0; j < n; j++) {
        e1 = i4_min(edge[0+j*m], edge[1+j*m]);
        e2 = i4_max(edge[0+j*m], edge[1+j*m]);
        edge[0+j*m] = e1;
        edge[1+j*m] = e2;
    }

    //
    //  Ascending sort the column array.
    //
    i4col_sort_a(m, n, edge);

    //
    //  Records which appear twice are internal edges and can be ignored.
    //
    node_boundary = (char*) malloc(node_num*sizeof(char));
    for ( i = 0; i < node_num; i++ ) {
        node_boundary[i] = 0;
    }

    j = 0;
    while(j < 3*element_num) {
        j = j + 1;

        if (j == 3*element_num) {
            for(i = 0; i < m; i++) {
                node_boundary[edge[i+(j-1)*m]-1] = 1;
            }
            break;
        }

        equal = 1;
        for(i = 0; i < m; i++) {
            if(edge[i+(j-1)*m] != edge[i+j*m]) {
                equal = 0;
            }
        }

        if(equal) {
            j = j + 1;
        } else {
            for(i = 0; i < m; i++) {
                node_boundary[edge[i+(j-1)*m]-1] = 1;
            }
        }
    }

    return node_boundary;
}
//****************************************************************************80


void quad_rule ( int quad_num, double quad_w[], double quad_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_RULE sets the quadrature rule for assembly.
//
//  Discussion:
//
//    The quadrature rule is given for a reference element.
//
//      0 <= X,
//      0 <= Y, and
//      X + Y <= 1.
//
//      ^
//    1 | *
//      | |.
//    Y | | .
//      | |  .
//    0 | *---*
//      +------->
//        0 X 1
//
//    The rules have the following precision:
//
//    QUAD_NUM  Precision
//
//     1        1
//     3        2
//     4        3
//     6        4
//     7        5
//     9        6
//    13        7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int QUAD_NUM, the number of quadrature nodes.
//
//    Output, double QUAD_W[QUAD_NUM], the quadrature weights.
//
//    Output, double QUAD_XY[2*QUAD_NUM],
//    the coordinates of the quadrature nodes.
//
{
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;
    double h;
    double t;
    double u;
    double v;
    double w;

    if ( quad_num == 1 ) {
        quad_xy[0+0*2] = 1.0 / 3.0;
        quad_xy[1+0*2] = 1.0 / 3.0;

        quad_w[0] = 1.0;
    } else if ( quad_num == 3 ) {
        quad_xy[0+0*2] = 0.5;
        quad_xy[1+0*2] = 0.0;
        quad_xy[0+1*2] = 0.5;
        quad_xy[1+1*2] = 0.5;
        quad_xy[0+2*2] = 0.0;
        quad_xy[1+2*2] = 0.5;

        quad_w[0] = 1.0 / 3.0;
        quad_w[1] = 1.0 / 3.0;
        quad_w[2] = 1.0 / 3.0;
    } else if ( quad_num == 4 ) {
        a =   6.0 / 30.0;
        b =  10.0 / 30.0;
        c =  18.0 / 30.0;

        d =  25.0 / 48.0;
        e = -27.0 / 48.0;

        quad_xy[0+0*2] = b;
        quad_xy[1+0*2] = b;
        quad_xy[0+1*2] = c;
        quad_xy[1+1*2] = a;
        quad_xy[0+2*2] = a;
        quad_xy[1+2*2] = c;
        quad_xy[0+3*2] = a;
        quad_xy[1+3*2] = a;

        quad_w[0] = e;
        quad_w[1] = d;
        quad_w[2] = d;
        quad_w[3] = d;
    } else if ( quad_num == 6 ) {
        a = 0.816847572980459;
        b = 0.091576213509771;
        c = 0.108103018168070;
        d = 0.445948490915965;
        v = 0.109951743655322;
        w = 0.223381589678011;

        quad_xy[0+0*2] = a;
        quad_xy[1+0*2] = b;
        quad_xy[0+1*2] = b;
        quad_xy[1+1*2] = a;
        quad_xy[0+2*2] = b;
        quad_xy[1+2*2] = b;
        quad_xy[0+3*2] = c;
        quad_xy[1+3*2] = d;
        quad_xy[0+4*2] = d;
        quad_xy[1+4*2] = c;
        quad_xy[0+5*2] = d;
        quad_xy[1+5*2] = d;

        quad_w[0] = v;
        quad_w[1] = v;
        quad_w[2] = v;
        quad_w[3] = w;
        quad_w[4] = w;
        quad_w[5] = w;
    } else if ( quad_num == 7 ) {
        a = 1.0 / 3.0;
        b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
        c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
        d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
        e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
        u = 0.225;
        v = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
        w = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;

        quad_xy[0+0*2] = a;
        quad_xy[1+0*2] = a;
        quad_xy[0+1*2] = b;
        quad_xy[1+1*2] = c;
        quad_xy[0+2*2] = c;
        quad_xy[1+2*2] = b;
        quad_xy[0+3*2] = c;
        quad_xy[1+3*2] = c;
        quad_xy[0+4*2] = d;
        quad_xy[1+4*2] = e;
        quad_xy[0+5*2] = e;
        quad_xy[1+5*2] = d;
        quad_xy[0+6*2] = e;
        quad_xy[1+6*2] = e;

        quad_w[0] = u;
        quad_w[1] = v;
        quad_w[2] = v;
        quad_w[3] = v;
        quad_w[4] = w;
        quad_w[5] = w;
        quad_w[6] = w;
    } else if ( quad_num == 9 ) {
        a = 0.124949503233232;
        b = 0.437525248383384;
        c = 0.797112651860071;
        d = 0.165409927389841;
        e = 0.037477420750088;

        u = 0.205950504760887;
        v = 0.063691414286223;

        quad_xy[0+0*2] = a;
        quad_xy[1+0*2] = b;
        quad_xy[0+1*2] = b;
        quad_xy[1+1*2] = a;
        quad_xy[0+2*2] = b;
        quad_xy[1+2*2] = b;
        quad_xy[0+3*2] = c;
        quad_xy[1+3*2] = d;
        quad_xy[0+4*2] = c;
        quad_xy[1+4*2] = e;
        quad_xy[0+5*2] = d;
        quad_xy[1+5*2] = c;
        quad_xy[0+6*2] = d;
        quad_xy[1+6*2] = e;
        quad_xy[0+7*2] = e;
        quad_xy[1+7*2] = c;
        quad_xy[0+8*2] = e;
        quad_xy[1+8*2] = d;

        quad_w[0] = u;
        quad_w[1] = u;
        quad_w[2] = u;
        quad_w[3] = v;
        quad_w[4] = v;
        quad_w[5] = v;
        quad_w[6] = v;
        quad_w[7] = v;
        quad_w[8] = v;
    } else if ( quad_num == 13 ) {
        h = 1.0 / 3.0;
        a = 0.479308067841923;
        b = 0.260345966079038;
        c = 0.869739794195568;
        d = 0.065130102902216;
        e = 0.638444188569809;
        f = 0.312865496004875;
        g = 0.048690315425316;

        w = -0.149570044467670;
        t =  0.175615257433204;
        u =  0.053347235608839;
        v =  0.077113760890257;

        quad_xy[0+ 0*2] = h;
        quad_xy[1+ 0*2] = h;
        quad_xy[0+ 1*2] = a;
        quad_xy[1+ 1*2] = b;
        quad_xy[0+ 2*2] = b;
        quad_xy[1+ 2*2] = a;
        quad_xy[0+ 3*2] = b;
        quad_xy[1+ 3*2] = b;

        quad_xy[0+ 4*2] = c;
        quad_xy[1+ 4*2] = d;
        quad_xy[0+ 5*2] = d;
        quad_xy[1+ 5*2] = c;
        quad_xy[0+ 6*2] = d;
        quad_xy[1+ 6*2] = d;

        quad_xy[0+ 7*2] = e;
        quad_xy[1+ 7*2] = f;
        quad_xy[0+ 8*2] = e;
        quad_xy[1+ 8*2] = g;
        quad_xy[0+ 9*2] = f;
        quad_xy[1+ 9*2] = e;

        quad_xy[0+10*2] = f;
        quad_xy[1+10*2] = g;
        quad_xy[0+11*2] = g;
        quad_xy[1+11*2] = e;
        quad_xy[0+12*2] = g;
        quad_xy[1+12*2] = f;

        quad_w[ 0] = w;
        quad_w[ 1] = t;
        quad_w[ 2] = t;
        quad_w[ 3] = t;
        quad_w[ 4] = u;
        quad_w[ 5] = u;
        quad_w[ 6] = u;
        quad_w[ 7] = v;
        quad_w[ 8] = v;
        quad_w[ 9] = v;
        quad_w[10] = v;
        quad_w[11] = v;
        quad_w[12] = v;
    } else {
        puts("\nQUAD_RULE - Fatal error!");
        printf("  No rule is available of order QUAD_NUM = %d\n", quad_num);
        exit(-1);
    }
    return;
}
//****************************************************************************80


void reference_to_physical_t3(double t[2*3], int n, double *ref, double *phy)

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are literally halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |.
//    |  | .
//    S  |  .
//    |  |   .
//    |  |    .
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of objects to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
    int i;
    int j;

    for ( i = 0; i < 2; i++ ) {
        for ( j = 0; j < n; j++ ) {
            phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                         + t[i+1*2] * ref[0+j*2]
                         + t[i+2*2] * ref[1+j*2];
        }
    }

    return;
}
//****************************************************************************80


int basis_one_t3(double t[2*3], int i, double p[2], double *qi, double *dqidx, double *dqidy)

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_ONE_T3 evaluates basis functions for a linear triangular element.
//
//  Discussion:
//
//    The routine is given the coordinates of the nodes of a triangle.
//
//           3
//          / .
//         /   .
//        /     .
//       1-------2
//
//    It evaluates the linear basis function Q(I)(X,Y) associated with
//    node I, which has the property that it is a linear function
//    which is 1 at node I and zero at the other two nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the nodes.
//
//    Input, int I, the index of the desired basis function.
//    I should be between 1 and 3.
//
//    Input, double P[2], the coordinates of the point where
//    the basis function is to be evaluated.
//
//    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
//    and its X and Y derivatives.
//
{
    double area;
    int ip1;
    int ip2;

    area = t[0+0*2] * (t[1+1*2] - t[1+2*2])
           + t[0+1*2] * (t[1+2*2] - t[1+0*2])
           + t[0+2*2] * (t[1+0*2] - t[1+1*2]);

    if ( area == 0.0 ) {
        puts("\nBASIS_ONE_T3 - Fatal error!");
        puts("  Element has zero area.");
        printf("  Area = %f\n\n", area);
        printf("  Node 1: ( %f, %f )\n", t[0+0*2], t[1+0*2]);
        printf("  Node 2: ( %f, %f )\n", t[0+1*2], t[1+1*2]);
        printf("  Node 3: ( %f, %f )\n", t[0+2*2], t[1+2*2]);
        return -1;
    }

    if ( i < 1 || 3 < i ) {
        puts("\nBASIS_ONE_T3 - Fatal error!");
        puts("  Basis index I is not between 1 and 3.");
        printf("  I = %d\n", i);
        return -2;
    }

    ip1 = i4_wrap(i + 1, 1, 3);
    ip2 = i4_wrap(i + 2, 1, 3);

    *qi = ( ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] )
            * ( p[1]           - t[1+(ip1-1)*2] )
            - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] )
            * ( p[0]           - t[0+(ip1-1)*2] ) ) / area;

    *dqidx = - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) / area;
    *dqidy =   ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) / area;

    return 0;
}
//****************************************************************************80


void assemble_poisson(int node_num, double *node_xy, int element_num, int *element_node,
                      int quad_num, int ib, double *A, double *f)

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_POISSON assembles the system for the Poisson equation.
//
//  Discussion:
//
//    The matrix is known to be banded.  A special matrix storage format
//    is used to reduce the space required.  Details of this format are
//    discussed in the routine DGB_FA.
//
//    Note that a 3 point quadrature rule, which is sometimes used to
//    assemble the matrix and right hand side, is just barely accurate
//    enough for simple problems.  If you want better results, you
//    should use a quadrature rule that is more accurate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Output, double A(3*IB+1,NODE_NUM), the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Output, double F(NODE_NUM), the right hand side.
//
//  Local parameters:
//
//    Local, double BI, DBIDX, DBIDY, the value of some basis function
//    and its first derivatives at a quadrature point.
//
//    Local, double BJ, DBJDX, DBJDY, the value of another basis
//    function and its first derivatives at a quadrature point.
//
{
    int i;
    int j;
    int node;
    int element;
    int quad;
    int ret;

    double *phys_h;
    double *phys_k;
    double *phys_rhs;
    double *phys_xy;
    double *quad_w;
    double *quad_xy;
    double *w;

    double t3[2*3];
    double p[2];

    double bi;
    double bj;
    double dbidx;
    double dbidy;
    double dbjdx;
    double dbjdy;

    int test;
    int basis;
    // double temp_val;
    double area;

    phys_h = (double*) malloc(quad_num*sizeof(double));
    phys_k = (double*) malloc(quad_num*sizeof(double));
    phys_rhs = (double*) malloc(quad_num*sizeof(double));
    phys_xy = (double*) malloc(2*quad_num*sizeof(double));
    quad_w = (double*) malloc(quad_num*sizeof(double));
    quad_xy = (double*) malloc(2*quad_num*sizeof(double));
    w = (double*) malloc(quad_num*sizeof(double));

    //
    //  Initialize the arrays to zero.
    //
    for(node = 0; node < node_num; node++) {
        f[node] = 0.0;
    }
    for(node = 0; node < node_num; node++) {
        for ( i = 0; i < 3*ib+1; i++ ) {
            A[i+node*(3*ib+1)] = 0.0;
        }
    }

    //
    //  Get the quadrature weights and nodes.
    //
    quad_rule(quad_num, quad_w, quad_xy);

    //
    //  Add up all quantities associated with the ELEMENT-th element.
    //
    for(element = 0; element < element_num; element++) {
        for(j = 0; j < 3; j++) {
            for(i = 0; i < 2; i++) {
                //  Make a copy of the element.
                t3[i+j*2] = node_xy[i+(element_node[j+element*3]-1)*2];
            }
        }

        //
        //  Map the quadrature points QUAD_XY to points PHYS_XY in the physical element.
        //
        reference_to_physical_t3(t3, quad_num, quad_xy, phys_xy);

        area = r8_abs(triangle_area_2d(t3));
        for(quad = 0; quad < quad_num; quad++) {
            w[quad] = quad_w[quad]*area;
        }

        rhs(quad_num, phys_xy, phys_rhs);
        h_coef(quad_num, phys_xy, phys_h);
        k_coef(quad_num, phys_xy, phys_k);

        //
        //  Consider the QUAD-th quadrature point.
        //
        for(quad = 0; quad < quad_num; quad++) {
            p[0] = phys_xy[0+quad*2];
            p[1] = phys_xy[1+quad*2];

            //
            //  Consider the TEST-th test function.
            //
            //  We generate an integral for every node associated with an unknown.
            //  But if a node is associated with a boundary condition, we do nothing.
            //
            for(test = 1; test <= 3; test++) {
                i = element_node[test-1+element*3];

                ret = basis_one_t3(t3, test, p, &bi, &dbidx, &dbidy);
                if(ret != 0) {
                    continue;
                }
                f[i-1] += w[quad] * phys_rhs[quad] * bi;

                //
                //  Consider the BASIS-th basis function, which is used to form the
                //  value of the solution function.
                //
                for(basis = 1; basis <= 3; basis++) {
                    j = element_node[basis-1+element*3];

                    ret = basis_one_t3(t3, basis, p, &bj, &dbjdx, &dbjdy);
                    if(ret != 0) {
                        continue;
                    }
                    A[i-j+2*ib+(j-1)*(3*ib+1)] +=
                        w[quad]* (phys_h[quad]*(dbidx*dbjdx + dbidy*dbjdy) + phys_k[quad]*bj*bi);
                }
            }
        }
    }

    free(phys_h);
    free(phys_k);
    free(phys_rhs);
    free(phys_xy);
    free(quad_w);
    free(quad_xy);
    free(w);

    return;
}
//****************************************************************************80


void dirichlet_apply(int node_num, double *node_xy, int *node_condition, int ib,
                     double *A, double *f)

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_APPLY accounts for Dirichlet boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int NODE_CONDITION[NODE_NUM], reports the condition
//    used to set the unknown associated with the node.
//    0, unknown.
//    1, finite element equation.
//    2, Dirichlet condition;
//    3, Neumann condition.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by
//    NODE_NUM coefficient matrix, stored in a compressed format; on output,
//    the matrix has been adjusted for Dirichlet boundary conditions.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//    On output, the right hand side has been adjusted for Dirichlet
//    boundary conditions.
//
{
    int column;
    int column_high;
    int column_low;
    int DIRICHLET = 2;
    int node;
    double *node_bc;

    node_bc = (double*) malloc(node_num*sizeof(double));

    dirichlet_condition(node_num, node_xy, node_bc);

    for(node = 0; node < node_num; node++) {
        if(node_condition[node] == DIRICHLET) {
            column_low = i4_max(node + 1 - ib, 1);
            column_high = i4_min(node + 1 + ib, node_num);

            for(column = column_low; column <= column_high; column++) {
                A[node+1-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
            }
            A[2*ib+node*(3*ib+1)] = 1.0;

            f[node] = node_bc[node];
        }
    }

    free(node_bc);

    return;
}
//****************************************************************************80


double* residual_poisson(int node_num, double *node_xy, int *node_condition,
                         int element_num, int *element_node, int quad_num, int ib,
                         double *A, double *f, double *node_u)

//****************************************************************************80
//
//  Purpose:
//
//    RESIDUAL_POISSON evaluates the residual for the Poisson equation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the
//    coordinates of nodes.
//
//    Input, int NODE_CONDITION[NODE_NUM], reports the condition
//    used to set the unknown associated with the node.
//    0, unknown.
//    1, finite element equation.
//    2, Dirichlet condition;
//    3, Neumann condition.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Workspace, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Workspace, double F[NODE_NUM], the right hand side.
//
//    Input, double NODE_U[NODE_NUM], the value of the solution
//    at each node.
//
//    Output, double NODE_R[NODE_NUM], the finite element
//    residual at each node.
//
//  Local parameters:
//
//    Local, double BI, DBIDX, DBIDY, the value of some basis function
//    and its first derivatives at a quadrature point.
//
//    Local, double BJ, DBJDX, DBJDY, the value of another basis
//    function and its first derivatives at a quadrature point.
//
{
    int i;
    int j;
    int node;
    int element;
    int quad;

    double *phys_h;
    double *phys_k;
    double *phys_rhs;
    double *phys_xy;
    double *quad_w;
    double *quad_xy;
    double *w;

    double t3[2*3];
    double p[2];

    double bi;
    double bj;
    double dbidx;
    double dbidy;
    double dbjdx;
    double dbjdy;

    int test;
    int basis;
    int ret;
    double area;
    // double temp_val;
    double *node_r;

    phys_h = (double*) malloc(quad_num*sizeof(double));
    phys_k = (double*) malloc(quad_num*sizeof(double));
    phys_rhs = (double*) malloc(quad_num*sizeof(double));
    phys_xy = (double*) malloc(2*quad_num*sizeof(double));
    quad_w = (double*) malloc(quad_num*sizeof(double));
    quad_xy = (double*) malloc(2*quad_num*sizeof(double));
    w = (double*) malloc(quad_num*sizeof(double));

    //
    //  Initialize the arrays to zero.
    //
    for(node = 0; node < node_num; node++) {
        f[node] = 0.0;
    }
    for(node = 0; node < node_num; node++) {
        for ( i = 0; i < 3*ib+1; i++ ) {
            A[i+node*(3*ib+1)] = 0.0;
        }
    }
    //
    //  Get the quadrature weights and nodes.
    //
    quad_rule(quad_num, quad_w, quad_xy);

    //
    //  The actual values of A and F are determined by summing up
    //  contributions from all the elements.
    //
    for(element = 0; element < element_num; element++) {
        //
        //  Make a copy of the element.
        //
        for (j = 0; j < 3; j++) {
            for (i = 0; i < 2; i++) {
                //  Make a copy of the element.
                t3[i+j*2] = node_xy[i+(element_node[j+element*3]-1)*2];
            }
        }

        //
        //  Map the quadrature points QUAD_XY to points XY in the physical element.
        //
        reference_to_physical_t3(t3, quad_num, quad_xy, phys_xy);

        area = r8_abs(triangle_area_2d(t3));

        for(quad = 0; quad < quad_num; quad++) {
            w[quad] = quad_w[quad] * area;
        }

        rhs(quad_num, phys_xy, phys_rhs);
        h_coef(quad_num, phys_xy, phys_h);
        k_coef(quad_num, phys_xy, phys_k);

        //
        //  Consider a quadrature point QUAD, with coordinates (X,Y).
        //
        for(quad = 0; quad < quad_num; quad++) {
            p[0] = phys_xy[0+quad*2];
            p[1] = phys_xy[1+quad*2];

            //
            //  Consider one of the basis functions, which will play the
            //  role of test function in the integral.
            //
            //  We generate an integral for every node associated with an unknown.
            //  But if a node is associated with a boundary condition, we do nothing.
            //
            for(test = 1; test <= 3; test++) {
                i = element_node[test-1+element*3];

                ret = basis_one_t3(t3, test, p, &bi, &dbidx, &dbidy);
                if(ret != 0) {
                    continue;
                }
                f[i-1] += w[quad]*phys_rhs[quad]*bi;

                //
                //  Consider another basis function, which is used to form the
                //  value of the solution function.
                //
                for(basis = 1; basis <= 3; basis++) {
                    j = element_node[basis-1+element*3];

                    ret = basis_one_t3(t3, basis, p, &bj, &dbjdx, &dbjdy);
                    if(ret != 0) {
                        continue;
                    }
                    A[i-j+2*ib+(j-1)*(3*ib+1)] +=
                        w[quad]*(phys_h[quad]*(dbidx*dbjdx + dbidy*dbjdy) + phys_k[quad]*bj*bi);
                }
            }
        }
    }

    //
    //  Apply boundary conditions.
    //
    dirichlet_apply(node_num, node_xy, node_condition, ib, A, f);

    //
    //  Compute A*U.
    //
    node_r = dgb_mxv(node_num, node_num, ib, ib, A, node_u);

    //
    //  Set RES = A * U - F.
    //
    for(node = 0; node < node_num; node++) {
        node_r[node] = node_r[node] - f[node];
    }

    free(phys_h);
    free(phys_k);
    free(phys_rhs);
    free(phys_xy);
    free(quad_w);
    free(quad_xy);
    free(w);

    return node_r;
}
//****************************************************************************80


double* fem2d_poisson(TriangulateIO *mesh, double **residual_out)

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_POISSON.
//
//  Discussion:
//
//    FEM2D_POISSON solves the Poisson equation
//
//      -DEL H(X,Y) DEL U(X,Y) + K(X,Y) * U(X,Y) = F(X,Y)
//
//    in a triangulated region in the plane.
//
//    Along the boundary of the region, Dirichlet conditions
//    are imposed:
//
//      U(X,Y) = G(X,Y)
//
//    The code uses continuous piecewise linear basis functions on
//    triangles.
//
//  Problem specification:
//
//    The user defines the geometry by supplying two data files
//    which list the node coordinates, and list the nodes that make up
//    each element.
//
//    The user specifies the right hand side of the Dirichlet boundary
//    conditions by supplying a function
//
//      void dirichlet_condition ( int node_num, double node_xy[2*node_num],
//        double node_bc[node_num] )
//
//    The user specifies the coefficient function H(X,Y) of the Poisson
//    equation by supplying a routine of the form
//
//      void h_coef ( node_num, node_xy, node_h )
//
//    The user specifies the coefficient function K(X,Y) of the Poisson
//    equation by supplying a routine of the form
//
//      void k_coef ( int node_num, double node_xy[2*node_num],
//        double node_k[node_num] )
//
//    The user specifies the right hand side of the Poisson equation
//    by supplying a routine of the form
//
//      void rhs ( int node_num, double node_xy[2*node_num],
//        double node_f[node_num] )
//
//  Usage:
//
//    fem2d_poisson 'prefix'
//
//    where 'prefix' is the common filename prefix so that:
//
//    * prefix_nodes.txt contains the coordinates of the nodes;
//    * prefix_elements.txt contains the indices of nodes forming each element.
//
//    Files created include:
//
//    * prefix_nodes.eps, an image of the nodes;
//    * prefix_elements.eps, an image of the elements;
//    * prefix_solution.txt, the value of the solution at every node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*NODE_NUM], the coefficient matrix.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, integer ELEMENT_ORDER, the element order.
//
//    Local, double F[NODE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, char NODE_BOUNDARY[NODE_NUM], is TRUE if the node is
//    found to lie on the boundary of the region.
//
//    Local, int NODE_CONDITION[NODE_NUM],
//    indicates the condition used to determine the variable at a node.
//    0, there is no condition (and no variable) at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Local, int NODE_NUM, the number of nodes.
//
//    Local, double NODE_R[NODE_NUM], the residual error.
//
//    Local, double NODE_U[NODE_NUM], the finite element coefficients.
//
//    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Local, integer QUAD_NUM, the number of quadrature points used for
//    assembly.  This is currently set to 3, the lowest reasonable value.
//    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
//    of QUAD_NUM greater than 3 may be appropriate.
//
{
    // double temp;
    char debug = 0;
    int node;

    double *node_xy;
    char *node_boundary;
    int *node_condition;
    int node_num;

    int dim_num;
    int *element_node;
    int element_num;
    int element_order;
    // int element_show;

    int quad_num = 7;

    int ib;
    int ierr;
    int job;

    double *A;
    double *f;
    int *pivot;
    double *node_u;
    double *node_r;

    // char * prefix;
    // char * element_eps_filename;
    // char * element_filename;
    // char * node_eps_filename;
    // char * node_filename;
    // char * solution_filename;
    // char node_label;
    // int node_show;

    timestamp();
    puts("");
    puts("FEM2D_POISSON:");
    puts("  C++ version:");
    puts("\n");
    puts("  Solution of the Poisson equation in an arbitrary region in 2 dimensions.\n");
    puts("  - DEL H(x,y) DEL U(x,y) + K(x,y) * U(x,y) = F(x,y) in the region\n");
    puts("                                     U(x,y) = G(x,y) on the boundary.\n");
    puts("  The finite element method is used, with triangular elements,");
    puts("  which must be a 3 node linear triangle.");

    //
    //  Read the node coordinate file.
    //
    dim_num = 2;
    node_num = mesh->numberofpoints;
    puts("");
    printf("  Number of nodes =          %d\n", node_num);
    node_condition = (int*) malloc(node_num*sizeof(int));
    node_xy = mesh->pointlist;
    // r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, 2, 10, "  First 10 nodes");

    //
    //  Read the element description file.
    //
    element_order = mesh->numberofcorners;
    element_num = mesh->numberoftriangles;
    puts("");
    printf("  Element order =            %d\n", element_order);
    printf("  Number of elements =       %d\n", element_num);
    if(element_order != 3) {
        puts("\nFEM2D_POISSON - Fatal error!");
        printf("  The input triangulation has order %d\n", element_order);
        puts("  However, a triangulation of order 3 is required.");
        return NULL;
    }
    // element_node = mesh->trianglelist;
    job = element_order*element_num;
    element_node = (int*) malloc(job*sizeof(int));
    for(node=0; node<job; node++) {
        element_node[node] = mesh->trianglelist[node] + 1;
    }
    // i4mat_transpose_print_some(3, element_num, element_node, 1, 1, 3, 10, "  First 10 elements");

    puts("");
    printf("  Quadrature order =          %d\n", quad_num);

    //
    //  Determine which nodes are boundary nodes and which have a
    //  finite element unknown.  Then set the boundary values.
    //
    node_boundary = triangulation_order3_boundary_node(node_num, element_num, element_node);
    // if(debug) {
    //     lvec_print(node_num, node_boundary, "    Node  Boundary?");
    // }

    //
    //  Determine the node conditions.
    //  For now, we'll just assume all boundary nodes are Dirichlet.
    //
    for(node = 0; node < node_num; node++) {
        if(node_boundary[node]) {
            node_condition[node] = 2;
        } else {
            node_condition[node] = 1;
        }
    }

    //
    //  Determine the bandwidth of the coefficient matrix.
    //
    ib = bandwidth(element_num, element_node);
    puts("");
    printf("  The matrix half bandwidth is %d\n", ib);
    printf("  The matrix bandwidth is      %d\n", 2*ib + 1);
    printf("  The storage bandwidth is     %d\n", 3*ib + 1);

    //
    //  Allocate space for the coefficient matrix A and right hand side F.
    //
    A = (double*) malloc((3*ib+1)*node_num*sizeof(double));
    f = (double*) malloc(node_num*sizeof(double));
    pivot = (int*) malloc(node_num*sizeof(int));

    //
    //  Assemble the finite element coefficient matrix A and the right-hand side F.
    //
    assemble_poisson(node_num, node_xy, element_num, element_node, quad_num, ib, A, f);
    if(debug) {
        //     dgb_print_some(node_num, node_num, ib, ib, A, 1, 1, 10, 10,
        //                    "  Initial block of Finite Element matrix A:");
        //     r8vec_print_some(node_num, f, 1, 10, "  Part of right hand side vector:");
        FILE *out;

        out = fopen("Matrix A Assembled.txt", "wb");
        m_foutput(out, A, node_num, 3*ib+1);
        fclose(out);

        out = fopen("Vector f Assembled.txt", "wb");
        v_foutput(out, f, node_num);
        fclose(out);
    }

    //
    //  Adjust the linear system to account for Dirichlet boundary conditions.
    //
    dirichlet_apply(node_num, node_xy, node_condition, ib, A, f);
    if(debug) {
        //     dgb_print_some(node_num, node_num, ib, ib, A, 1, 1, 10, 10,
        //                    "  Finite Element matrix A after boundary adjustments:");
        //     r8vec_print_some(node_num, f, 1, 10, "  Part of right hand side vector:");

        FILE *out;

        out = fopen("Matrix A Dirichlet.txt", "wb");
        m_foutput(out, A, node_num, 3*ib+1);
        fclose(out);

        out = fopen("Vector f Dirichlet.txt", "wb");
        v_foutput(out, f, node_num);
        fclose(out);
    }

    //
    //  Solve the linear system using a banded solver.
    //
    puts("\nSolving linear system: pivoting...");
    ierr = dgb_fa(node_num, ib, ib, A, pivot);
    if(debug) {
        FILE *out;
        out = fopen("Matrix A Pivot.txt", "wb");
        m_foutput(out, A, node_num, 3*ib+1);
        fclose(out);
    }

    job = 0;
    puts("\nSolving linear system: solving...");
    node_u = dgb_sl(node_num, ib, ib, A, pivot, f, job);
    // r8vec_print_some(node_num, node_u, 1, 10, "  Part of the solution vector U:");
    //
    //  Write an ASCII file that can be read into MATLAB.
    //
    // r8mat_write(solution_filename, 1, node_num, node_u);
    // if ( debug ) {
    //      v_output(node_u);
    //      r8vec_print_some(node_num, node_u, 1, 10, "  Part of the solution vector:");
    // }

    puts("\nCalculating residuals...");
    node_r = residual_poisson(node_num, node_xy, node_condition, element_num, element_node,
                              quad_num, ib, A, f, node_u);
    // temp = v_max(node_r, NULL);
    // puts("");
    // printf("  Maximum absolute residual = %f\n", temp);

    //
    //  Deallocate memory.
    //
    // delete [] node_xy;
    free(node_boundary);
    free(node_condition);
    free(element_node);
    *residual_out = node_r; // free(node_r);
    free(A);
    free(f);
    free(pivot);
    // V_FREE(node_u);

    //
    //  Terminate.
    //
    puts("\nFEM2D_POISSON:");
    puts("  Normal end of execution.\n");
    timestamp();

    return node_u;
}
//****************************************************************************80


int bandwidth(int element_num, int *element_node)

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth of the coefficient matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input,  ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int BANDWIDTH, the half bandwidth of the matrix.
//
{
    int element;
    int global_i;
    int global_j;
    int local_i;
    int local_j;
    int nhba;

    nhba = 0;

    for(element = 0; element < element_num; element++) {
        for(local_i = 0; local_i < 3; local_i++) {
            global_i = element_node[local_i+element*3];
            for(local_j = 0; local_j < 3; local_j++) {
                global_j = element_node[local_j+element*3];
                nhba = i4_max(nhba, abs(global_j - global_i));
            }
        }
    }

    return nhba;
}
//****************************************************************************80


void timestamp(void)

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
#define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
    printf("%s\n", time_buffer);

    return;

#undef TIME_SIZE
}
//****************************************************************************80


int i4_modp(int i, int j)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of integer division.
//
//  Formula:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Discussion:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
    int value;

    if ( j == 0 ) {
        puts("\nI4_MODP - Fatal error!");
        printf("  I4_MODP ( I, J ) called with J = %d\n", j);
        exit(-1);
    }

    value = i % j;
    if(value < 0) {
        value = value + abs(j);
    }

    return value;
}
//****************************************************************************80*


int i4_wrap(int ival, int ilo, int ihi)

//****************************************************************************80*
//
//  Purpose:
//
//    I4_WRAP forces an integer to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I4_WRAP
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
    int jhi;
    int jlo;
    int value;
    int wide;

    jlo = i4_min ( ilo, ihi );
    jhi = i4_max ( ilo, ihi );

    wide = jhi + 1 - jlo;

    if ( wide == 1 ) {
        value = jlo;
    } else {
        value = jlo + i4_modp ( ival - jlo, wide );
    }

    return value;
}
//****************************************************************************80


void i4col_swap(int m, int n, int *a, int icol1, int icol2)

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

    int i;
    int t;

    //
    //  Check.
    //
    if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET ) {
        puts("\nI4COL_SWAP - Fatal error!");
        puts("  ICOL1 is out of range.");
        exit(-1);
    }

    if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET ) {
        puts("\nI4COL_SWAP - Fatal error!");
        puts("  ICOL2 is out of range.");
        exit(-1);
    }

    if ( icol1 == icol2 ) {
        return;
    }
    for ( i = 0; i < m; i++ ) {
        t                     = a[i+(icol1-OFFSET)*m];
        a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
        a[i+(icol2-OFFSET)*m] = t;
    }

    return;
# undef OFFSET
}
//****************************************************************************80


int i4col_compare ( int m, int n, int *a, int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
    int k;

    //
    //  Check.
    //
    if(i < 1) {
        puts("\nI4COL_COMPARE - Fatal error!");
        printf("  Column index I = %d is less than 1.\n", i);
        exit(-1);
    }

    if(n < i) {
        puts("\nI4COL_COMPARE - Fatal error!");
        printf("  N = %d is less than column index I = %d.\n", n, i);
        exit(-1);
    }

    if(j < 1) {
        puts("\nI4COL_COMPARE - Fatal error!");
        printf("  Column index J = %d is less than 1.\n", j);
        exit(-1);
    }

    if(n < j) {
        puts("\nI4COL_COMPARE - Fatal error!");
        printf("  N = %d is less than column index J = %d.\n", n, j);
        exit(-1);
    }

    if(i == j) {
        return 0;
    }

    k = 1;
    while(k <= m) {
        if(a[k-1+(i-1)*m] < a[k-1+(j-1)*m]) {
            return (-1);
        } else if(a[k-1+(j-1)*m] < a[k-1+(i-1)*m]) {
            return 1;
        }
        k = k + 1;
    }

    return 0;
}
//****************************************************************************80


void i4col_sort_a(int m, int n, int *a)

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
    int i;
    int indx;
    int isgn;
    int j;

    //
    //  Initialize.
    //
    i = 0;
    indx = 0;
    isgn = 0;
    j = 0;

    //
    //  Call the external heap sorter.
    //
    for ( ; ; ) {
        sort_heap_external(n, &indx, &i, &j, isgn);

        //
        //  Interchange the I and J objects.
        //
        if ( 0 < indx ) {
            i4col_swap(m, n, a, i, j);
        }

        //
        //  Compare the I and J objects.
        //
        else if ( indx < 0 ) {
            isgn = i4col_compare(m, n, a, i, j);
        } else if ( indx == 0 ) {
            break;
        }

    }

    return;
}
//****************************************************************************80


double r8_abs (double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
    if ( 0.0 <= x ) {
        return x;
    } else {
        return ( -x );
    }
}
//****************************************************************************80


double triangle_area_2d(double t[2*3])

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counterclockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.
//
{
    double area;

    area = 0.5 * (t[0+0*2] * (t[1+1*2] - t[1+2*2]) +
                  t[0+1*2] * (t[1+2*2] - t[1+0*2]) +
                  t[0+2*2] * (t[1+0*2] - t[1+1*2]));

    return area;
}
//****************************************************************************80


void sort_heap_external(int n, int *indx, int *i, int *j, int isgn)

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis and Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
    static int i_save = 0;
    static int j_save = 0;
    static int k = 0;
    static int k1 = 0;
    static int n1 = 0;

    //
    //  INDX = 0: This is the first call.
    //
    if ( *indx == 0 ) {

        i_save = 0;
        j_save = 0;
        k = n / 2;
        k1 = k;
        n1 = n;
    }

    //
    //  INDX < 0: The user is returning the results of a comparison.
    //
    else if ( *indx < 0 ) {
        if ( *indx == -2 ) {
            if ( isgn < 0 ) {
                i_save = i_save + 1;
            }
            j_save = k1;
            k1 = i_save;
            *indx = -1;
            *i = i_save;
            *j = j_save;
            return;
        }

        if ( 0 < isgn ) {
            *indx = 2;
            *i = i_save;
            *j = j_save;
            return;
        }

        if ( k <= 1 ) {
            if ( n1 == 1 ) {
                i_save = 0;
                j_save = 0;
                *indx = 0;
            } else {
                i_save = n1;
                j_save = 1;
                n1 = n1 - 1;
                *indx = 1;
            }
            *i = i_save;
            *j = j_save;
            return;
        }
        k = k - 1;
        k1 = k;
    }

    //
    //  0 < INDX: the user was asked to make an interchange.
    //
    else if ( *indx == 1 ) {
        k1 = k;
    }

    for ( ; ; ) {

        i_save = 2 * k1;

        if ( i_save == n1 ) {
            j_save = k1;
            k1 = i_save;
            *indx = -1;
            *i = i_save;
            *j = j_save;
            return;
        } else if ( i_save <= n1 ) {
            j_save = i_save + 1;
            *indx = -2;
            *i = i_save;
            *j = j_save;
            return;
        }

        if ( k <= 1 ) {
            break;
        }

        k = k - 1;
        k1 = k;
    }

    if ( n1 == 1 ) {
        i_save = 0;
        j_save = 0;
        *indx = 0;
        *i = i_save;
        *j = j_save;
    } else {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
        *i = i_save;
        *j = j_save;
    }

    return;
}
//****************************************************************************80
