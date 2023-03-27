#ifndef NULL
#define NULL 0
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// for running tests
// #include "tests/test_meschach.h"
#include "triangle_main/tests/test_triangle.h"
#include "plplot_main/tests/test_plplot.h"

// for running actual program
#include "triangle_main/get_mesh.h"
#include "fem/fem2d_poisson.h"  // includes "fem/prob_def.h"
#include "plplot_main/interpolated_plot.h"


void v_foutput(FILE * const out, double * const data, const int num_data)
{
    int i;

    for(i=0; i<num_data; i++) {
        fprintf(out, "%f\n", data[i]);
    }
}

void m_foutput(FILE * const out, double * const data, const int m, const int n)
{
    int i, j;

    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            fprintf(out, "%f\t", data[i*n+j]);
        }
        fputs("\n", out);
    }
}

void u_foutput(FILE * const out, double * const points_list, double * const u_list,
               double * const residuals, const int num_points)
{
    int i;
    double x, y;

    fputs("Solution:\n", out);
    fputs("x:\ty:\tu:\tresidual/error:\n", out);
    for(i=0; i<num_points; i++) {
        x = points_list[2*i];
        y = points_list[2*i+1];
        // u = u_list[i];

        fprintf(out, "%15.7g\t%15.7g\t%15.7g\t%15.7g\n", x, y, u_list[i], residuals[i]);
    }
}

int main(int argc, char *argv[])
{
    // fprintf(stderr, "running tests...\n");
    // solve_random_matrix();  /* testing meschach */
    // triangle_oldtest();     /* testing triangle */
    // test_plplot();          // testing plplot

    if(argc < 2) {
        printf("Usage: %s <geometry file>\n", argv[0]);
        return -1;
    }

    puts("\nRunning FEM Solver\n\n");

    int i;
    char temp_buffer[512], temp_buffer2[512];
    TriangulateIO *mesh;
    double *u, *residual, *exact, *error;
    FILE *out;

    sprintf(temp_buffer, "%s.poly", argv[1]);
    mesh = get_mesh(argv[1], temp_buffer);
    if(mesh == NULL) {
        return -1;
    }

    u = fem2d_poisson(mesh, &residual);
    // sprintf(temp_buffer, "%s_u.txt", argv[1]);
    // out = fopen(temp_buffer, "wb");
    // v_foutput(out, u, mesh->numberofpoints);
    // fclose(out);

    // sprintf(temp_buffer, "%s_residual.txt", argv[1]);
    // out = fopen(temp_buffer, "wb");
    // v_foutput(out, residual, mesh->numberofpoints);
    // fclose(out);

    // outputting solution in text and plots
    sprintf(temp_buffer, "%s_solution.txt", argv[1]);
    out = fopen(temp_buffer, "wb");
    u_foutput(out, mesh->pointlist, u, residual, mesh->numberofpoints);
    fclose(out);

    puts("\nPlotting Results\n");

    sprintf(temp_buffer, "%s_density_plot.png", argv[1]);
    density_plot(mesh->pointlist, u, mesh->numberofpoints, mesh->segmentlist,
                 mesh->numberofsegments, argv[1], temp_buffer);

    sprintf(temp_buffer, "%s_mesh_plot.png", argv[1]);
    mesh_plot(mesh->pointlist, u, mesh->numberofpoints, mesh->segmentlist,
              mesh->numberofsegments, argv[1], temp_buffer);

#ifndef LAPLACE_EQUATION
    // calculating error/deviation from exact solution
    exact = (double*) malloc(mesh->numberofpoints*sizeof(double));
    error = (double*) malloc(mesh->numberofpoints*sizeof(double));
    u_exact(mesh->numberofpoints, mesh->pointlist, exact);
    for(i=0; i<mesh->numberofpoints; i++) {
        // by error, ideally we mean relative error, but...
        // error[i] = fabs((exact[i]-u[i])/exact[i]*100.0);    // ...may run into infinities
        error[i] = fabs(exact[i]-u[i]);
    }

    // making exact solution plots
    sprintf(temp_buffer2, "%s (interpolated exact solution)", argv[1]);

    sprintf(temp_buffer, "%s_exact_density_plot.png", argv[1]);
    density_plot(mesh->pointlist, exact, mesh->numberofpoints, mesh->segmentlist,
                 mesh->numberofsegments, temp_buffer2, temp_buffer);

    sprintf(temp_buffer, "%s_exact_mesh_plot.png", argv[1]);
    mesh_plot(mesh->pointlist, exact, mesh->numberofpoints, mesh->segmentlist,
              mesh->numberofsegments, temp_buffer2, temp_buffer);


    // making error plots
    sprintf(temp_buffer2, "%s error", argv[1]);

    sprintf(temp_buffer, "%s_error_density_plot.png", argv[1]);
    density_plot(mesh->pointlist, error, mesh->numberofpoints, mesh->segmentlist,
                 mesh->numberofsegments, temp_buffer2, temp_buffer);

    sprintf(temp_buffer, "%s_error_mesh_plot.png", argv[1]);
    mesh_plot(mesh->pointlist, error, mesh->numberofpoints, mesh->segmentlist,
              mesh->numberofsegments, temp_buffer2, temp_buffer);


    // outputting exact solution in text
    sprintf(temp_buffer, "%s_exact_error.txt", argv[1]);
    out = fopen(temp_buffer, "wb");
    u_foutput(out, mesh->pointlist, exact, error, mesh->numberofpoints);
    fclose(out);
#endif

    return 0;
}


