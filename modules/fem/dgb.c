#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dgb.h"


int i4_min(int i1, int i2)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
//
{
    if ( i1 < i2 ) {
        return i1;
    } else {
        return i2;
    }

}
//****************************************************************************80


int i4_max(int i1, int i2)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
//
{
    if ( i2 < i1 ) {
        return i1;
    } else {
        return i2;
    }

}
//****************************************************************************80


int dgb_fa(int n, int ml, int mu, double *a, int *pivot)

//****************************************************************************80
//
//  Purpose:
//
//    DGB_FA performs a LINPACK-style PLU factorization of a DGB matrix.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2003
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.
//    On output, A has been overwritten by the LU factors.
//
//    Output, int PIVOT[N], the pivot vector.
//
//    Output, int DGB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
    int col = 2 * ml + mu + 1;
    int i;
    int i0;
    int j;
    int j0;
    int j1;
    int ju;
    int jz;
    int k;
    int l;
    int lm;
    int m;
    int mm;
    double t;

    m = ml + mu + 1;

    //
    //  Zero out the initial fill-in columns.
    //
    j0 = mu + 2;
    j1 = i4_min ( n, m ) - 1;

    for ( jz = j0; jz <= j1; jz++ ) {
        i0 = m + 1 - jz;
        for ( i = i0; i <= ml; i++ ) {
            a[i-1+(jz-1)*col] = 0.0;
        }
    }

    jz = j1;
    ju = 0;

    for ( k = 1; k <= n-1; k++ ) {
        //
        //  Zero out the next fill-in column.
        //
        jz = jz + 1;
        if ( jz <= n ) {
            for ( i = 1; i <= ml; i++ ) {
                a[i-1+(jz-1)*col] = 0.0;
            }
        }

        //
        //  Find L = pivot index.
        //
        lm = i4_min ( ml, n-k );
        l = m;

        for ( j = m+1; j <= m + lm; j++ ) {
            if ( fabs ( a[l-1+(k-1)*col] ) < fabs ( a[j-1+(k-1)*col] ) ) {
                l = j;
            }
        }

        pivot[k-1] = l + k - m;

        //
        //  Zero pivot implies this column already triangularized.
        //
        if (a[l-1+(k-1)*col] == 0.0) {
            puts("\nDGB_FA - Fatal error!");
            printf("  Zero pivot on step %d\n", k);
            // return k;
            continue;
        }

        //
        //  Interchange if necessary.
        //
        t                = a[l-1+(k-1)*col];
        a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
        a[m-1+(k-1)*col] = t;

        //
        //  Compute multipliers.
        //
        for ( i = m+1; i <= m+lm; i++ ) {
            a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
        }

        //
        //  Row elimination with column indexing.
        //
        ju = i4_max(ju, mu + pivot[k-1]);
        ju = i4_min(ju, n);
        mm = m;

        for ( j = k+1; j <= ju; j++ ) {
            l = l - 1;
            mm = mm - 1;

            if ( l != mm ) {
                t                 = a[l-1+(j-1)*col];
                a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
                a[mm-1+(j-1)*col] = t;
            }
            for ( i = 1; i <= lm; i++ ) {
                a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col]
                                      + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
            }
        }
    }

    pivot[n-1] = n;

    if(a[m-1+(n-1)*col] == 0.0) {
        puts("\nDGB_FA - Fatal error!\n");
        printf("  Zero pivot on step %d\n", n);
        return n;
    }

    return 0;
}
//****************************************************************************80


double* dgb_sl(int n, int ml, int mu, double *a_lu, int *pivot, double *b, int job)

//****************************************************************************80
//
//  Purpose:
//
//    DGB_SL solves a system factored by DGB_FA.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from DGB_FA.
//
//    Input, int PIVOT[N], the pivot vector from DGB_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double DGB_SL[N], the solution.
//
{
    int col = 2 * ml + mu + 1;
    int i;
    int k;
    int l;
    int la;
    int lb;
    int lm;
    int m;
    double t;
    double *x;

    x = (double*) malloc(n*sizeof(double));

    for ( i = 0; i < n; i++ ) {
        x[i] = b[i];
    }
    m = mu + ml + 1;

    //
    //  Solve A * x = b.
    //
    if ( job == 0 ) {
        //
        //  Solve L * Y = B.
        //
        if ( 1 <= ml ) {
            for ( k = 1; k <= n-1; k++ ) {
                lm = i4_min ( ml, n-k );
                l = pivot[k-1];

                if ( l != k ) {
                    t      = x[l-1];
                    x[l-1] = x[k-1];
                    x[k-1] = t;
                }
                for ( i = 1; i <= lm; i++ ) {
                    x[k+i-1] = x[k+i-1] + x[k-1] * a_lu[m+i-1+(k-1)*col];
                }
            }
        }

        //
        //  Solve U * X = Y.
        //
        for ( k = n; 1 <= k; k-- ) {
            x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
            lm = i4_min ( k, m ) - 1;
            la = m - lm;
            lb = k - lm;
            for ( i = 0; i <= lm-1; i++ ) {
                x[lb+i-1] = x[lb+i-1] - x[k-1] * a_lu[la+i-1+(k-1)*col];
            }
        }
    }

    //
    //  Solve A' * X = B.
    //
    else {
        //
        //  Solve U' * Y = B.
        //
        for ( k = 1; k <= n; k++ ) {
            lm = i4_min ( k, m ) - 1;
            la = m - lm;
            lb = k - lm;
            for ( i = 0; i <= lm-1; i++ ) {
                x[k-1] = x[k-1] - x[lb+i-1] * a_lu[la+i-1+(k-1)*col];
            }
            x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
        }
        //
        //  Solve L' * X = Y.
        //
        if ( 1 <= ml ) {
            for ( k = n-1; 1 <= k; k-- ) {
                lm = i4_min ( ml, n-k );
                for ( i = 1; i <= lm; i++ ) {
                    x[k-1] = x[k-1] + x[k+i-1] * a_lu[m+i-1+(k-1)*col];
                }
                l = pivot[k-1];

                if ( l != k ) {
                    t      = x[l-1];
                    x[l-1] = x[k-1];
                    x[k-1] = t;
                }
            }
        }
    }

    return x;
}
//****************************************************************************80


double* dgb_mxv(int m, int n, int ml, int mu, double *a, double *x)

//****************************************************************************80
//
//  Purpose:
//
//    DGB_MXV multiplies a DGB matrix times a vector.
//
//  Discussion:
//
//    The DGB storage format is for an M by N banded matrix, with lower
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
//    extra superdiagonals, which may be required to store nonzero entries
//    generated during Gaussian elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    LINPACK and LAPACK storage of general band matrices requires
//    an extra ML upper diagonals for possible fill in entries during
//    Gauss elimination.  This routine does not access any entries
//    in the fill in diagonals, because it assumes that the matrix
//    has NOT had Gauss elimination applied to it.  If the matrix
//    has been Gauss eliminated, then the routine DGB_MU must be
//    used instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongarra, Bunch, Cleve Moler, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the DGB matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DGB_MXV[M], the product A * x.
//
{
    double *b;
    int i;
    int j;
    int jhi;
    int jlo;

    b = (double*) malloc(m*sizeof(double));

    for ( i = 1; i <= m; i++ ) {
        b[i-1] = 0.0;
        jlo = i4_max ( 1, i - ml );
        jhi = i4_min ( n, i + mu );
        for ( j = jlo; j <= jhi; j++ ) {
            b[i-1] += a[i-j+ml+mu+(j-1)*(2*ml+mu+1)] * x[j-1];
        }
    }

    return b;
}
//****************************************************************************80
