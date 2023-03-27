#ifndef _PROB_DEF_H_
#define _PROB_DEF_H_

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


void dirichlet_condition(int node_num, double *node_xy, double *node_g);
void h_coef(int node_num, double *node_xy, double *node_h);
void k_coef(int node_num, double *node_xy, double *node_k);
void rhs(int node_num, double *node_xy, double *node_rhs);

#ifndef LAPLACE_EQUATION
// for figuring out the error in the numerical solution
void u_exact(int node_num, double *node_xy, double *node_u_exact);
#endif

#endif

