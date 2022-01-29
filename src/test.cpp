#include<iostream>
#include<ctime>
#include "matrixReal.hpp"
using namespace std;

int main(){
    clock_t start_time,end_time;
    start_time = clock();
    int i = 0;
    for(int k = 0; k<10000; k++){
        i++;
    }
    end_time = clock();
    cout<<"Time taken:\t"<<(end_time-start_time)<<endl;
    
    // Matrix input & output test
    // Matrix A(3,3);
    // cin>>A;
    // cout<<A<<endl;
    // cout<<"test 1 passed"<<endl;
    
    // randomize matrix test
    // srand(time(0));
    // Matrix B(3,3);
    // B.randMatrix((double)rand()/INT_MAX,(double)rand()/INT_MAX);
    // cout<<B<<endl;
    // cout<<"test 2 passed"<<endl;

    // Transpose matrix test
    // srand(time(0));
    // Matrix B(6,3), C;
    // B.randMatrix((double)rand()/INT_MAX,(double)rand()/INT_MAX);
    // cout<<B<<endl;
    // C = B.transpose();
    // cout<<endl;
    // cout<<C<<endl;
    // cout<<"test 3 passed"<<endl;

    // Diagnol of matrix test.
    // srand(time(0));
    // Matrix B(3,5);
    // B.randMatrix((double)rand()/INT_MAX,(double)rand()/INT_MAX);
    // cout<<B<<endl;
    // double *x;
    // x = B.diag();
    // for(int i = 0; i<3; i++)
    //     cout<<x[i]<<"\t";
    // cout<<endl;
    // cout<<"test 4 passed"<<endl;

    // determinant test
    // srand(time(0));
    // Matrix B(3,3);
    // cin>>B;
    // double determinant = B.det();
    // cout<<determinant<<endl;
    // cout<<"test 5 passed"<<endl;

    // rank test
    // srand(time(0));
    // Matrix B(3,3);
    // B.randMatrix((double)rand()/rand(),(double)rand()/rand());
    // int rank = B.rank();
    // cout<<rank<<endl;
    // cout<<"test 5 passed"<<endl;

    // rand determinant test
    // srand(time(0));
    // Matrix B(3,3);
    // B.randMatrix((double)rand()/rand(),(double)rand()/rand());
    // double det3x3 = ((B(0,0)*B(1,1)*B(2,2) + B(0,1)*B(1,2)*B(2,0) + B(0,2)*B(2,1)*B(1,0)) - (B(2,0)*B(1,1)*B(0,2)+B(2,1)*B(1,2)*B(0,0)+B(2,2)*B(0,1)*B(1,0)));
    // cout<<B<<endl;
    // double detB = B.det();
    // cout<<"error % = "<<(double)((det3x3 - detB)/detB)*100<<endl;
    // cout<<"test 6 passed"<<endl;

    // Inverse test

    // pinv test
    return 0;
}

