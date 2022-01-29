# C++ Matrix Code
## Features:
* A small c++ real matrix header file with only using standrad library, easy to modify.
* Defined class "Matrix" with private variables num_rows, num_cols & double **M.
* Included all algorithmic operations like +,-,* between two matrices 
---
## Available Functions  
* <b>randMatrix(a,b)</b> - generates random matrix between a & b
* <b>rank</b> - Calculate matrix rank (Cholesky decomposition)
* <b>det</b> - returns the determinant of nxn matrix
* <b>inv</b> - returns inverse of nxn Matrix
---
## Examples
* few code snippets that'll help you understand a little better on how to use .hpp file.
### 1 -> Random Matrix, rank, determinant & inverse
```cpp
#include "src/matrixReal.hpp"
using namespace std;
int main(){
    cout<<"test 0 passed"<<endl;  // Library initialisation test
    return 0;
}
```
---
## References
* Pierre Courrieu, Fast Computation of Moore-Penrose Inverse Matrices, https://arxiv.org/abs/0804.4809
