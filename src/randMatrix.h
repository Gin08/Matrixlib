#include <stdio.h>
#include <stdlib.h>

double **randmatrix(int rows, int columns, double lowerbound, double upperbound);

double **randmatrix(int rows, int columns, double lowerbound, double upperbound) {
   double **m;

   m = (double **) calloc ( rows, sizeof( double *));
   for (int i = 0; i < rows; i++)
   {
      m[i] = (double *) calloc ( columns, sizeof( double));
      for (int j = 0; j < columns; j++)
      {
          m[i][j] = ((rand() % (int)(upperbound - lowerbound + 1)) + lowerbound);
      }
   }
   return m;
}

