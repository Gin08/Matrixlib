#include"stdio.h"
#include"src\ginv.h"
#include"src\randMatrix.h"

int main()
{
    double G[6][6] = {1,2.3,-2.1,4.5,9.6,6.32,0,2.8,8.2,-2,1,5.02,6.3,5.3,1.1,8.9,0,-0.8,0.2,-9.9,6.9,7.7,5,1,2,9,7.8,-4.7,-2.4,0,0,0,0,0,0,0};
    double **X;
    X = matrix2D(6,6);
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            X[i][j] = G[i][j];
        }
    }
    double **ginv;
    ginv = matrix2D(6,6);
    pinvF(X,6,6,ginv);
    printf("\n ginv matrix : \n");
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            printf("%f,",ginv[i][j]);
        }
        printf("\n");
    }
    double **Xnew = randmatrix(8,9,1.23,10.45);
    printf("\n random matrix :\n");
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            printf("%f, ",Xnew[i][j]);
        }
        printf("\n");
    }
    free2D(Xnew,8);
    free2D(X,6);
    free2D(ginv,6);
    return 0;
}
