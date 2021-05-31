/*
 * MatrixOps.h
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

void MatTranspose(double** A, int rows, int columns, double** result);
void MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2, double** result);
void Matinverse(double** num, double** X_inv, int f);
double Matdeterminant(double** a, int k);


#endif /* MATRIXOPS_H_ */
