#ifndef _INTERPOLATED_PLOT_H_
#define _INTERPOLATED_PLOT_H_


int density_plot(double * const points_list, double * const u, const int num_points,
                 int * const line_list, const int num_lines, char * const name,
                 char * const filename);

int mesh_plot(double * const points_list, double * const u, const int num_points,
              int * const line_list, const int num_lines, char * const name,
              char * const filename);

#endif
