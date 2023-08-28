#include <math.h>

#include "fem/prob_def.h"


#ifndef LAPLACE_EQUATION

void dirichlet_condition(int node_num, double *node_xy, double *node_g)
{
    int i;
    double x, y;

    for(i=0; i<node_num; i++) {
        x = node_xy[2*i];
        y = node_xy[2*i+1];

        node_g[i] = (x*x + y*y)/(500.0*500.0);
    }
}

void h_coef(int node_num, double *node_xy, double *node_h)
{
    int i;

    for(i=0; i<node_num; i++) {
        node_h[i] = 1.0;
    }
}

void k_coef(int node_num, double *node_xy, double *node_k)
{
    int i;

    for(i=0; i<node_num; i++) {
        node_k[i] = 1.0;
    }
}

void rhs(int node_num, double *node_xy, double *node_rhs)
{
    int i;
    double x, y;

    for(i=0; i<node_num; i++) {
        x = node_xy[2*i];
        y = node_xy[2*i+1];

        node_rhs[i] = (x*x + y*y - 4.0)/(500.0*500.0);
    }
}

void u_exact(int node_num, double *node_xy, double *node_u_exact)
{
    int i;
    double x, y;

    for(i=0; i<node_num; i++) {
        x = node_xy[2*i];
        y = node_xy[2*i+1];

        node_u_exact[i] = (x*x + y*y)/(500.0*500.0);
    }
}

#else   // do laplace equation here

static void find_min_max(double * const points_list, const int num_points,
                        double * const min_out, double * const max_out)
{
    int i, j;
    double x_min = 9999999999;
    double x_max = -9999999999;

    for(i=0; i<num_points; i++) {
        j = 2*i;
        if(points_list[j] < x_min) {
            x_min = points_list[j];
        } else if(points_list[j] > x_max) {
            x_max = points_list[j];
        }
    }

    *min_out = x_min;
    *max_out = x_max;
}

void dirichlet_condition(int node_num, double *node_xy, double *node_g)
{
    int i;
    double x_min, x_max, x_mid;
    double x, y;

    find_min_max(node_xy, node_num, &x_min, &x_max);
    x_mid = (x_min + x_max)/2.0;

    for(i=0; i<node_num; i++) {
        x = node_xy[2*i];
        y = node_xy[2*i+1];

        if(x < x_mid) {
            node_g[i] = 0;
        } else {
            node_g[i] = 100.0;
        }
    }

}

void h_coef(int node_num, double *node_xy, double *node_h)
{
    int i;
    for(i=0; i<node_num; i++) {
        node_h[i] = 1.0;
    }
}

void k_coef(int node_num, double *node_xy, double *node_k)
{
    int i;
    for(i=0; i<node_num; i++) {
        node_k[i] = 0.0;
    }
}

void rhs(int node_num, double *node_xy, double *node_rhs)
{
    int i;
    for(i=0; i<node_num; i++) {
        node_rhs[i] = 0.0;
    }
}

#endif
