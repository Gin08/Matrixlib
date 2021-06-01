#include <math.h>
#include"Memory.h"

void MatTranspose(double** A, int rows, int columns, double** result);
void MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2, double** result);
void Matinverse(double** num, double** X_inv, int f);
double Matdeterminant(double** a, int k);


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

void Matinverse(double** num, double** X_inv, int f)
{
    /*
     * This function computes the inverse of matrix num of size f(rows)xf(columns)
     * and returns it by storing into X_inv of size f(rows)xf(columns)
     */

    double **b;
    b = matrix2D(f,f);
    double **fac;
    fac = matrix2D(f,f);
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
    double **b1;
    b1 = matrix2D(f,f);
    double d;

    for (i = 0; i < f; i++)
    {
        for (j = 0; j < f; j++)
        {

            b1[i][j] = fac[j][i];
        }
    }
    double **flk;
    flk = matrix2D(f,f);
    for (i = 0; i < f; i++)
    {
        for (j = 0; j < f; j++)
        {
            flk[i][j] = num[i][j];
        }
    }
    d = Matdeterminant(flk, f);

    for (int i1 = 0; i1 < f; i1++)
    {
        for (int j1 = 0; j1 < f; j1++)
        {
            X_inv[i1][j1] = (b1[i1][j1]) / (d);
        }
    }
    free2D(b,f);
    free2D(fac,f);
    free2D(b1,f);
    free2D(flk,f);
}

double Matdeterminant(double** a, int k)
{
    /*
     * This function returns the determinant of matrix a of size k(rows)xk(columns)
     */

    static float det;
    double s = 1;
    double **b;
    b = matrix2D(k,k);

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
    free2D(b,k);
    return (det);

}
//************** Determinant over *****************//

