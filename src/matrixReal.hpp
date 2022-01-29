#include<iostream>
#include<math.h>
#include<time.h>
#include<vector>
using namespace std;

class Matrix{
private:
    double **M;
    int num_rows, num_cols;
public:
    Matrix();
    Matrix(int, int);
    Matrix(double**,int,int);
    void input();
    void disp();
    Matrix operator+(Matrix);
    Matrix operator-(Matrix);
    Matrix operator*(Matrix);
    Matrix operator+(double);
    void randMatrix(double, double);

    friend Matrix operator+(double, Matrix);
    friend int row_size(Matrix);
    friend int col_size(Matrix);
    friend ostream& operator<<(ostream&, Matrix&);
    friend istream& operator>>(istream&, Matrix&);
    
    double& operator()(int, int);
    void setVal(int , int , double);
    Matrix copy();
    void clear();
    Matrix transpose();
    double* diag();
    double det();
    int rank();
    Matrix inv();
    Matrix pinv();
    Matrix ginv();
    Matrix eigM();
    double* eigV();
};

Matrix :: Matrix(){
    num_rows = 0;
    num_cols = 0;
}

Matrix  :: Matrix(int m, int n)
{
    num_rows = m;
    num_cols = n;
    M = new double*[num_rows];
    for(int i=0;i<num_rows;i++)
        M[i] = new double[num_cols];
}

Matrix  :: Matrix(double **X,int m, int n)
{
    num_rows = m;
    num_cols = n;
    M = new double*[num_rows];
    for(int i=0;i<num_rows;i++){
        M[i] = new double[num_cols];
        for(int j=0; j<num_cols;j++){
            M[i][j]=X[i][j];
        }
    }
}

double* Matrix::diag(){
    double* D;
    D = new double[num_rows];
    for(int i = 0; i<num_rows; i++)
        D[i] = M[i][i];
    return D;
}

void Matrix ::input()
{
    cout << "Enter the elements row wise" << endl;
    for(int i=0; i<this->num_rows; i++)
    {
        for(int j=0; j<this->num_cols; j++){
            cin >> M[i][j];
        }
    }
}

void Matrix::disp()
{
    for(int i=0;i<num_rows;i++)
    {
        for(int j=0;j<num_cols;j++)
            cout << M[i][j] << "\t";
        cout << endl;
    }
}

Matrix Matrix::operator+(Matrix A)
{
    Matrix P;
    if(num_rows != A.num_rows || num_cols != A.num_cols)
    {
        cout << "Matrix addition is invalid\n";
        exit(0);
    }
    else
    {
        P = Matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<num_cols;j++)
                P.M[i][j] = M[i][j] + A.M[i][j];
        return P;
    }
}

Matrix Matrix::operator+(double a)
{
    Matrix P;
    P = Matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] + a;
    return P;
    
}

Matrix operator+(double a, Matrix A){
    Matrix P;
    P = Matrix(A.num_rows, A.num_cols);
    for(int i=0 ; i < A.num_rows ; i++)
        for(int j=0 ; j < A.num_cols ; j++)
            P.M[i][j] = A.M[i][j] + a;
    return P;
}

Matrix Matrix::operator-(Matrix A)
{
    Matrix P;
    if(num_rows != A.num_rows || num_cols != A.num_cols)
    {
        cout << "Matrix subtraction is invalid\n";
        exit(0);
    }
    else
    {
        P = Matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<num_cols;j++)
                P.M[i][j] = M[i][j] - A.M[i][j];
        return P;
    }
}

void Matrix::randMatrix(double f_min, double f_max)
{
    srand(time(0));
    // Matrix<T> A(n_rows,n_cols);
    for(int i = 0; i<this->num_rows; i++){
        for(int j = 0; j<this->num_cols; j++){
            M[i][j] = f_min + (f_max-f_min)*((double)rand()/INT_MAX);
        }
    }
    // return A;
}

Matrix Matrix::operator*(Matrix A)
{
    Matrix P;
    if (num_cols != A.num_rows)
    {
        cout << "Invalid Matrix multiplication ! \n";
        exit(0);
    }
    else
    {
        P = Matrix(num_rows, A.num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<A.num_cols;j++)
            {
                P.M[i][j]=0.0;
                for(int k=0;k<num_cols;k++)
                    P.M[i][j] = P.M[i][j] + M[i][k]*A.M[k][j];
            }
        return P;
    }
}

int row_size (Matrix A)
{
    return A.num_rows;
}

int col_size (Matrix A)
{
    return A.num_cols;
}

ostream& operator<<(ostream& out, Matrix& A){
   A.disp();
   return out;
}

istream& operator>>(istream& in, Matrix& A){
    A.input();
    return in;
}

double& Matrix :: operator()(int m, int n)
{
    if (m < num_rows && n < num_cols)
        return M[m][n];
    else
    {
        cout << "Index out of bound !" << endl;
        exit(0);
    }
}

void Matrix::setVal(int i, int j, double n){
    this->M[i][j] = n;
}

Matrix Matrix::transpose()
{
    Matrix P(num_cols,num_rows);
    for(int i=0;i<num_cols;i++)
    {
        for(int j=0;j<num_rows;j++){
            P.M[i][j] = M[j][i];
        }
    }
    return P;
}

Matrix Matrix::copy()
{
    Matrix X(num_rows,num_cols);
    for(int i = 0; i< num_rows; i++)
        for(int j = 0; j< num_cols; j++)
            X.setVal(i,j,M[i][j]);
    return X;
}

void Matrix::clear()
{
    for(int i = 0; i< num_rows; i++){
        delete[] M[i];
    }
    delete[] M;
    num_cols = 0;
    num_rows = 0;
}

// LU decomposition-based matrix determinant calculation [*2][*3][*4]
double Matrix::det()
{
    if (num_cols != num_rows)
    {
        cout << "Determinant deosn't exists for non-square matrix ! \n";
        exit(0);
    }
    else
    {
        double detG = 0.0;
        int nSize = num_rows;
        std::vector< vector<double> > matLU;
        std::vector<std::size_t> permuteLU;
        bool changeSign = false;

        for (int i = 0; i < nSize; ++i)
        {
            permuteLU.push_back(i);
        }

        for (int j = 0; j < nSize; ++j)
        {
            double maxv = 0.0;
            for (int i = j; i < nSize; ++i)
            {
                const double currentv = std::abs(M[permuteLU[i]][j]);
                if (currentv > maxv)
                {
                    maxv = currentv;
                    if (permuteLU[i] != permuteLU[j]) // swap rows
                    {
                        changeSign = !changeSign;
                        const int tmp = permuteLU[j];
                        permuteLU[j] = permuteLU[i];
                        permuteLU[i] = tmp;
                    }
                }
            }
        }

        for (int i = 0; i < nSize; ++i)
        {
            std::vector<double> temp_p;
            for (int j = 0; j < nSize; ++j)
            {
                temp_p.push_back(M[i][j]);
            }
            matLU.push_back(temp_p);
        }

        if (matLU[0][0] == 0.0)
        {
            return detG; // Singular matrix, det(G) = 0
        }

        for (int i = 1; i < nSize; ++i)
        {
            matLU[i][0] /= matLU[0][0];
        }

        for (int i = 1; i < nSize; ++i)
        {
            for (int j = i; j < nSize; ++j)
            {
                for (int k = 0; k < i; ++k)
                {
                    matLU[i][j] -= matLU[i][k] * matLU[k][j]; // Calculate U matrix
                }
            }
            if (matLU[i][i] == 0.0)
            {
                return detG; // Singular matrix, det(G) = 0
            }
            for (int k = i + 1; k < nSize; ++k)
            {
                for (int j = 0; j < i; ++j)
                {
                    matLU[k][i] -= matLU[k][j] * matLU[j][i]; // Calculate L matrix
                }
                matLU[k][i] /= matLU[i][i];
            }
        }

        detG = 1.0;
        if (changeSign)
        {
            detG = -1.0; // Change the sign of the determinant
        }
        for (int i = 0; i < nSize; ++i)
        {
            detG *= matLU[i][i]; // det(G) = det(L) * det(U). For triangular matrices, det(L) = prod(diag(L)) = 1, det(U) = prod(diag(U)), so det(G) = prod(diag(U))
        }
        return detG;
    }
}

// Calculate matrix rank (Cholesky decomposition) [*1]
int Matrix::rank()
{
    const double tolerance = 1.0e-9;
    int nSize = num_cols;
    Matrix matA;
    if (num_rows < nSize)
    {
        // A = G * G'
        nSize = num_rows;
        Matrix matB(nSize, nSize);
        for (int i = 0; i < nSize; ++i)
        {
            for (int j = 0; j < nSize; ++j)
            {
                for (int k = 0; k < num_cols; ++k)
                {
                    matB(i,j) += M[i][k] * M[j][k];
                }
            }
        }
        matA = matB.copy();
    }
    else
    {
        // A = G' * G
        Matrix matB(nSize, nSize);
        for (int i = 0; i < nSize; ++i)
        {
            for (int j = 0; j < nSize; ++j)
            {
                for (int k = 0; k < num_cols; ++k)
                {
                    matB(i,j) += M[k][i] * M[k][j];
                }
            }
        }
        matA = matB.copy();
    }
    // Full rank Cholesky decomposition of A
    double tol = std::abs(matA(0,0));
    for (int i = 0; i < nSize; ++i)
    {
        if (matA(i,i) > 0)
        {
            const double temp = std::abs(matA(i,i));
            if (temp < tol)
            {
                tol = temp;
            }
        }
    }
    tol *= tolerance;

    Matrix matL(nSize, nSize);
    int rankA = 0;
    for (int k = 0; k < nSize; ++k)
    {
        for (int i = k; i < nSize; ++i)
        {
            matL(i,rankA) = matA(i,k);
            for (int j = 0; j < rankA; ++j)
            {
                matL(i,rankA) -= matL(i,j) * matL(k,j);
            }
        }
        if (matL(k,rankA) > tol)
        {
            matL(k,rankA) = sqrt(matL(k,rankA));
            if (k < nSize)
            {
                for (int j = k + 1; j < nSize; ++j)
                {
                    matL(j,rankA) /= matL(k,rankA);
                }
            }
            ++rankA;
        }
    }
    return rankA; // rank(G) = rank(A)
}

// LU decomposition-based matrix inversion [*3][*4]
Matrix Matrix::inv()
{
    const bool usePermute = true;
    Matrix matLU;
    if (num_rows != num_cols)
    {
        std::cout << "Error when using inv: matrix is not square.\n";
        exit(0);
    }
    int nSize = num_rows;
    // ******************** Step 1: row permutation (swap diagonal zeros) ********************
    std::vector<int> permuteLU; // Permute vector
    for (int i = 0; i < nSize; ++i)
    {
        permuteLU.push_back(i); // Push back row index
    }

    if (usePermute) // Sort rows by pivot element
    {
        for (int j = 0; j < nSize; ++j)
        {
            double maxv = 0.0;
            for (int i = j; i < nSize; ++i)
            {
                const double currentv = std::abs(M[permuteLU[i]][j]);
                if (currentv > maxv) // Swap rows
                {
                    maxv = currentv;
                    const int tmp = permuteLU[j];
                    permuteLU[j] = permuteLU[i];
                    permuteLU[i] = tmp;
                }
            }
        }
        for (int i = 0; i < nSize; ++i)
        {
            for(int j = 0; j< nSize; ++j){
                matLU.setVal(i,j,M[permuteLU[i]][j]); // Make a permuted matrix with new row order
            }
        }
    }
    else
    {
        matLU = copy(); // Simply duplicate matrix
    }

    // ******************** Step 2: LU decomposition (save both L & U in matLU) ********************
    if (matLU(0,0) == 0.0)
    {
        std::cout << "Warning when using inv: matrix is singular.\n";
        // matLU.clear();
        // return matLU;
        exit(0);
    }
    for (int i = 1; i < nSize; ++i)
    {
        matLU(i,0) /= matLU(0,0); // Initialize first column of L matrix
    }
    for (int i = 1; i < nSize; ++i)
    {
        for (int j = i; j < nSize; ++j)
        {
            for (int k = 0; k < i; ++k)
            {
                matLU(i,j) -= matLU(i,k) * matLU(k,j); // Calculate U matrix
            }
        }
        if (matLU(i,i) == 0.0)
        {
            std::cout << "Warning when using inv: matrix is singular.\n";
            // matLU.clear();
            // return matLU;
            exit(0);
        }
        for (int k = i + 1; k < nSize; ++k)
        {
            for (int j = 0; j < i; ++j)
            {
                matLU(k,i) -= matLU(k,j) * matLU(j,i); // Calculate L matrix
            }
            matLU(k,i) /= matLU(i,i);
        }
    }

    // ******************** Step 3: L & U inversion (save both L^-1 & U^-1 in matLU_inv) ********************
    Matrix matLU_inv(nSize, nSize);

    // matL inverse & matU inverse
    for (int i = 0; i < nSize; ++i)
    {
        // L matrix inverse, omit diagonal ones
        matLU_inv(i,i) = 1.0;
        for (int k = i + 1; k < nSize; ++k)
        {
            for (int j = i; j <= k - 1; ++j)
            {
                matLU_inv(k,i) -= matLU(k,j) * matLU_inv(j,i);
            }
        }
        // U matrix inverse
        matLU_inv(i,i) = 1.0 / matLU(i,i);
        for (int k = i; k > 0; --k)
        {
            for (int j = k; j <= i; ++j)
            {
                matLU_inv(k - 1,i) -= matLU(k - 1,j) * matLU_inv(j,i);
            }
            matLU_inv(k - 1,i) /= matLU(k - 1,k - 1);
        }
    }

    // ******************** Step 4: Calculate G^-1 = U^-1 * L^-1 ********************
    // Lower part product
    for (int i = 1; i < nSize; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            const int jp = permuteLU[j]; // Permute column back
            matLU(i,jp) = 0.0;
            for (int k = i; k < nSize; ++k)
            {
                matLU(i,jp) += matLU_inv(i,k) * matLU_inv(k,j);
            }
        }
    }
    // Upper part product
    for (int i = 0; i < nSize; ++i)
    {
        for (int j = i; j < nSize; ++j)
        {
            const int jp = permuteLU[j]; // Permute column back
            matLU(i,jp) = matLU_inv(i,j);
            for (int k = j + 1; k < nSize; ++k)
            {
                matLU(i,jp) += matLU_inv(i,k) * matLU_inv(k,j);
            }
        }
    }
    return matLU; // Reused matLU as a result container
}