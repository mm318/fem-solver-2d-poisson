#ifndef _FEM2D_POISSON_H_
#define _FEM2D_POISSON_H_

#include "triangle/triangle.h"
// #include "meschach/matrix2.h"


// "prob_def.h" should define the following functions:
// void dirichlet_condition(int node_num, double *node_xy, double *node_rhs);
// void h_coef(int node_num, double *node_xy, double *node_h);
// void k_coef(int node_num, double *node_xy, double *node_k);
// void rhs(int node_num, double *node_xy, double *node_rhs);
#include "prob_def.h"


// FEM related functions
char* triangulation_order3_boundary_node(int node_num, int element_num, int *element_node);
void quad_rule(int quad_num, double *quad_w, double *quad_xy);
void reference_to_physical_t3(double t[2*3], int n, double *ref, double *phy);
int basis_one_t3(double t[2*3], int i, double p[2], double *qi, double *dqidx, double *dqidy);
void assemble_poisson(int node_num, double *node_xy, int element_num, int *element_node,
                      int quad_num, int ib, double *A, double *f);
void dirichlet_apply(int node_num, double *node_xy, int *node_condition, int ib, double *A, double *f);
double* residual_poisson(int node_num, double *node_xy, int *node_condition,
                         int element_num, int *element_node, int quad_num, int ib,
                         double *A, double *f, double *node_u);
double* fem2d_poisson(TriangulateIO *mesh, double **residual_out);


// matrix functions
int bandwidth(int element_num, int *element_node);
// int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
// double *dgb_mxv ( int m, int n, int ml, int mu, double a[], double x[] );
// double *dgb_sl ( int n, int ml, int mu, double a_lu[], int pivot[],
//                  double b[], int job );

// helper functions
void timestamp(void);

int i4_modp(int i, int j);
int i4_wrap(int ival, int ilo, int ihi);
void i4col_swap(int m, int n, int *a, int icol1, int icol2);
int i4col_compare(int m, int n, int *a, int i, int j);
void i4col_sort_a(int m, int n, int *a);

double r8_abs(double x);

double triangle_area_2d(double t[2*3]);
void sort_heap_external(int n, int *indx, int *i, int *j, int isgn);


// debug/output functions
// void lvec_print ( int n, char a[], char * title );
// void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char * title );
// void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
//                                   int ihi, int jhi, char * title );
// void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
//                                   int ihi, int jhi, char * title );
// void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
//                       int jlo, int ihi, int jhi, char * title );


// old deprecated functions

// char ch_cap ( char c );
// char ch_eqi ( char c1, char c2 );
// int ch_to_digit ( char c );
// int file_column_count ( char * input_filename );
// int file_row_count ( char * input_filename );
// int* i4mat_data_read ( char * input_filename, int m, int n );
// void i4mat_header_read ( char * input_filename, int *m, int *n );

// double r8_huge ( void );
// int r8_nint ( double x );
// double *r8mat_data_read ( char * input_filename, int m, int n );
// void r8mat_header_read ( char * input_filename, int *m, int *n );
// void r8mat_write ( char * output_filename, int m, int n, double table[] );
// double r8vec_amax ( int n, double a[] );
// int s_len_trim ( char * s );
// int s_to_i4 ( char * s, int *last, char *error );
// char s_to_i4vec ( char * s, int n, int ivec[] );
// double s_to_r8 ( char * s, int *lchar, char *error );
// char s_to_r8vec ( char * s, int n, double rvec[] );
// int s_word_count ( char * s );
// void solution_evaluate ( double xy[2], double t[2*3], double node_u[3], double *u,
//                          double *dudx, double *dudy );


#endif
