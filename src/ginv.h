#include"math.h"
#include"MatrixOps.h"

double sqrt(double x);

// Input : G -> 2 dimensional array/matrix, m rows & n columns
// Output: result -> must be 2 dimensional array/matrix, n rows & m columns

void pinvF(double** G, int m, int n, double** result)
{
    int mi = m; 
    int ni = n;
    double **GTranspose;
    GTranspose = matrix2D(n,m);
    MatTranspose(G,m,n,GTranspose);
    bool transpose = 0;
    double **A;
    int size_A;
    if(m<n)
    {
        size_A = m;
        transpose = 1;
        A = matrix2D(size_A,size_A);
        MatMult(G,m,n,GTranspose,n,m,A);
        n = m;
    }
    else
    {
        size_A = n;
        A = matrix2D(size_A,size_A);
        MatMult(GTranspose,n,m,G,m,n,A);
    }
    double *dA = (double *)malloc(size_A * sizeof(double));
    double minpdA;
    for(int idA = 0; idA<size_A;idA++)
    {
        dA[idA] = A[idA][idA];
        if(dA[idA]>0)
        {
            minpdA = dA[idA];
        }
    }
    for(int idA = 0; idA<size_A;idA++)
    {
        if(dA[idA]>0 && dA[idA]<minpdA)
        {
            minpdA = dA[idA];
        }
    }
    double tol = minpdA * 0.000000001;
    double **L;
    L = matrix2D(size_A,size_A);
    for(int Lr=0;Lr<size_A;Lr++)
    {
        for(int Lc=0;Lc<size_A;Lc++)
        {
            L[Lr][Lc] = 0.0;
        }
    }
    int r = 0;
    for(int k=0;k<n;k++)
    {
        r = r + 1;
        if (r == 1)
        {
            for (int index1 = k; index1 < n; index1++)
            {
                L[index1][(r-1)] = A[index1][k];
            }
        }
        else
        {
            double *Lresult;
            Lresult = (double*)malloc((n-k) * sizeof(double));
            for(int i = 0; i < (n-k); ++i)
            {
                Lresult[i] = 0.0;
            }
            for(int i = 0; i < (n-k); ++i)
            {
                for(int k1 = 0; k1 < r-1; ++k1)
                {
                    Lresult[i] += L[(i+k)][k1]*L[k][k1];
                }
            }
            for (int index1 = k; index1 < n; index1++)
            {
                L[index1][(r-1)] = A[index1][k] - Lresult[(index1-k)];
            }
        }
        if (L[k][(r-1)] > tol)
        {
            L[k][(r-1)] = sqrt(L[k][(r-1)]);
            if (k < (n))
            {
                for (int tempLk1 = k+1; tempLk1 < n; tempLk1++)
                {
                    L[tempLk1][(r-1)] = (L[tempLk1][(r-1)])/(L[k][(r-1)]);
                }
            }
        }
        else
        {
            r = r -1;
        }
    }
    double **Lfinal;
    Lfinal = matrix2D(size_A,r);
    double **LfinalTrans;
    LfinalTrans = matrix2D(r,size_A);
    for (int i = 0; i < size_A; i++)
    {
        for (int Lfinalc = 0; Lfinalc < r; Lfinalc++)
        {
            Lfinal[i][Lfinalc] = L[i][Lfinalc];
            LfinalTrans[Lfinalc][i] = Lfinal[i][Lfinalc];
        }
    }
    double **M1;
    M1 = matrix2D(r,r);
    MatMult(LfinalTrans,r,size_A,Lfinal,size_A,r,M1);
    double **M;
    M = matrix2D(r,r);
    Matinverse(M1,M,r);
    double **Y;
    if (transpose == 1)
    {
        double **GtL;
        double **GtLM;
        double **GtLMM;
        GtL = matrix2D(ni,r);
        GtLM = matrix2D(ni,r);
        GtLMM = matrix2D(ni,r);
        Y = matrix2D(ni,mi);
        MatMult(GTranspose,ni,m,Lfinal,m,r,GtL);
        MatMult(GtL,ni,r,M,r,r,GtLM);
        MatMult(GtLM,ni,r,M,r,r,GtLMM);
        MatMult(GtLMM,ni,r,LfinalTrans,r,m,Y);
        for (int fr = 0; fr < ni; fr++)
        {
            for (int fc = 0; fc < mi; fc++)
            {
                result[fr][fc] = Y[fr][fc];
            }
        }
        free2D(GtL,ni);
        free2D(GtLM,ni);
        free2D(GtLMM,ni);
        free2D(Y,ni);
    }
    else
    {
        double **LM;
        double **LMM;
        double **LMMLt;
        LM = matrix2D(ni,r);
        LMM = matrix2D(ni,r);
        LMMLt = matrix2D(ni,ni);
        Y = matrix2D(ni,m);
        MatMult(Lfinal,ni,r,M,r,r,LM);
        MatMult(LM,ni,r,M,r,r,LMM);
        MatMult(LMM,ni,r,LfinalTrans,r,ni,LMMLt);
        MatMult(LMMLt,ni,ni,GTranspose,ni,m,Y);
        for (int fr = 0; fr < ni; fr++)
        {
            for (int fc = 0; fc < mi; fc++)
            {
                result[fr][fc] = Y[fr][fc];
            }
        }
        free2D(LM,ni);
        free2D(LMM,ni);
        free2D(LMMLt,ni);
        free2D(Y,ni);
    }
    free2D(A,size_A);
    free2D(M1,r);
    free2D(GTranspose,n);
}
