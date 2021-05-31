/*
 * MatrixOps.c
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#include <math.h>
#include <stdlib.h>
#include"MatrixOps.hpp"

void MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2, double** result) {
    /*
     * This function performs matrix multiplication of two matrices A and B
     */

    for(int i = 0; i < rows1; ++i)
        for(int j = 0; j < columns2; ++j) {
            result[i][j] = 0.0;
        }

    for(int i = 0; i < rows1; ++i)
        for(int j = 0; j < columns2; ++j) {
            for(int k = 0; k < columns1; ++k)
                result[i][j] += A[i][k]*B[k][j];
        }
}


void MatTranspose(double** A, int rows, int columns, double** result) {
    /*
     * This function computes the transpose of matrix A
     */

    for(int i = 0; i < columns; ++i)
        for(int j = 0; j < rows; ++j)
            result[i][j] = A[j][i];

}

//************** Inverse matrix start *****************//
void Matinverse(double** num, double** X_inv, int f)
{
    double *b[f];
    for (int br = 0; br < f; br++)
    {
        b[br] = (double*)malloc(f*sizeof(double));
    }
    double *fac[f];
    for (int facr = 0; facr < f; facr++)
    {
        fac[facr] = (double*)malloc(f*sizeof(double));
    }
    int p, q, m, n, i, j;

    for (q = 0; q < f; q++)
    {
        for (p = 0; p < f; p++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++)
            {
                for (j = 0; j < f; j++)
                {
                    b[i][j] = 0;
                    if (i != q && j != p)
                    {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                            n++;
                        else
                        {
                            n = 0;

                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * Matdeterminant(b, f - 1);
        }

    }
    double *b1[f];
    for (int b1r = 0; b1r < f; b1r++)
    {
        b1[b1r] = (double*)malloc(f*sizeof(double));
    }
    double d;

    for (i = 0; i < f; i++)
    {
        for (j = 0; j < f; j++)
        {

            b1[i][j] = fac[j][i];
        }
    }
    double *flk[f];
    for (int flkr = 0; flkr < f; flkr++)
    {
        flk[flkr] = (double*)malloc(f*sizeof(double));
    }
    for (i = 0; i < f; i++)
    {
        for (j = 0; j < f; j++)
        {
            flk[i][j] = num[i][j];
        }
    }
    d = Matdeterminant(flk, f);
    //printf("%f",d);

    for (int i1 = 0; i1 < f; i1++)
    {
        for (int j1 = 0; j1 < f; j1++)
        {
            X_inv[i1][j1] = (b1[i1][j1]) / (d);
        }
    }
}
//************** Inverse matrix over *****************//

//************** Determinant start *****************//
double Matdeterminant(double** a, int k)
{

    static float det;
    double s = 1;
    double *b[k];
    for (int i = 0; i < k; i++)
    {
        b[i] = (double*)malloc(k*sizeof(double));
    }

    int i, j, m, n, c;

    if (k == 1)
    {
        double temp;
        temp = float(a[0][0]);
        det = temp;
        return det;
    }
    else
    {

        det = 0.0;

        for (c = 0; c < k; c++)
        {

            m = 0;

            n = 0;

            for (i = 0; i < k; i++)
            {

                for (j = 0; j < k; j++)
                {

                    b[i][j] = 0;

                    if (i != 0 && j != c)
                    {

                        b[m][n] = a[i][j];

                        if (n < (k - 2))

                            n++;
                        else
                        {

                            n = 0;

                            m++;

                        }

                    }

                }

            }

            det = det + s * (a[0][c] * Matdeterminant(b, k - 1));

            s = -1 * s;

        }

    }

    return (det);

}
//************** Determinant over *****************//
