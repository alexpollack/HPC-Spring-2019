//
//  main.cpp
//  hpc_hw1
//  question 3: solving the Laplace equation by Jacobi Method
//  Created by Alex Pollack on 2/1/19.
//  Copyright © 2019 Alex Pollack. All rights reserved.
//  g++ -std=c++11 -O3 02-memory.cpp && ./a.out

#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include "utils.h"

using namespace std;

void jacobi(double uk[], double f[], int N);
void GSM(double uk[], double f[], int N);

int main()
{
    int N;
    cout << "Enter dimension N: ";  //input dimension size
    cin >> N;
    double uk[N], f[N];
    
    for(int i = 0; i < N; i++)  //f(x) = 1
        f[i] = 1.0;
    for(int i = 0; i < N; i++)  //u_0 is the zero vector
        uk[i] = 0.0;
    
    cout << "Jacobi Method:\n";
    jacobi(uk, f, N);   //call Jacobi Method
    
    for(int i = 0; i < N; i++)  //u_0 is the zero vector
        uk[i] = 0.0;
    cout << "Gauss-Seidel Method:\n";
    GSM(uk, f, N);  //call Gauss-Seidel Method
    
    return 0;
}

void jacobi(double uk[], double f[], int N)
{
    int k = 0;
    Timer t;
    double *A = (double*)malloc(N * N * sizeof(double));
    double r = 1.0, r0= 1.0, sum = 0.0, uk1[N], Au[N], resid[N];
    for (int i = 0; i < N; i++) //create matrix A
    {
        for(int j = 0; j < N; j++)
        {
            if(i == j)
                *(A + i*N + j) = 2.0; //A[i][j] = 2.0;
            else if(j == i + 1 || j == i - 1)
                *(A + i*N + j) = -1.0;//A[i][j] = -1.0;
            else
                *(A + i*N + j) = 0.0; //A[i][j] = 0.0;
        }
    }
    //begin Jacobi Method
    t.tic();
    while(r/r0 > 1e-6)  //end method when residual has been reduced by 1e6 from first try
    {
        for (int i = 0; i < N; i++)
        {
            sum = 0.0;
            for(int j = 0; j < N; j++)
            {
                if(j != i)
                    sum = sum + *(A + i*N + j)*uk[j];
            }
            uk1[i] = (f[i] - sum)/(*(A + i*N + i));
        }
        for(int i = 0; i < N; i++)
            uk[i] = uk1[i]; //u^k gets new u^k+1 value for future computations
        k++;    //iteration count
        for (int i = 0; i < N; i++)
        {
            sum = 0.0;
            for(int j = 0; j < N; j++)
                sum = sum + *(A + i*N + j)*uk[j]; //matrix multiplication Au
            Au[i] = sum;
            resid[i] = pow((Au[i] - f[i]),2);   //residual ||Auk - f||
        }
        sum = 0.0;
        //finish computing residual ||Auk - f||
        for(int i = 0; i < N; i++)
            sum = sum + resid[i];
        r = sqrt(sum);
        if(k == 1)
            r0 = r;
        cout << k << "th iteration residual: " << r << endl;
        if(k > 5000)    //break if exceeds 5000 iterations
        {
            cout << "\n!!Iteration count exceeded 5000!!\n";
            break;
        }
    }
    cout << "\nTime to find approximate solution: " << t.toc() << endl;
    cout << "\nFinal residual after " << k << " iterations: ||Auk - f|| = " << r << endl << endl;
    free(A);
}

void GSM(double uk[], double f[], int N)
{
    int k = 0;
    Timer t;
    double *A = (double*)malloc(N * N * sizeof(double));
    double r = 1.0, r0= 1.0, sum = 0.0, Au[N], resid[N];
    for (int i = 0; i < N; i++) //create matrix A
    {
        for(int j = 0; j < N; j++)
        {
            if(i == j)
                *(A + i*N + j) = 2.0; //A[i][j] = 2.0;
            else if(j == i + 1 || j == i - 1)
                *(A + i*N + j) = -1.0;//A[i][j] = -1.0;
            else
                *(A + i*N + j) = 0.0; //A[i][j] = 0.0;
        }
    }
    //begin Gaus-Seidel Method
    t.tic();
    while(r/r0 > 1e-6)  //end method when residual has been reduced by 1e6 from first try
    {
        for (int i = 0; i < N; i++)
        {
            sum = 0.0;
            for(int j = 0; j < N; j++)
            {
                if(j != i)
                    sum = sum + *(A + i*N + j)*uk[j];
            }
            uk[i] = (f[i] - sum)/(*(A + i*N + i));  //uk is directly upaded as the method is calculated each iteration at i
        }
        k++;    //iteration count
    
        for (int i = 0; i < N; i++)
        {
            sum = 0.0;
            for(int j = 0; j < N; j++)
                sum = sum + *(A + i*N + j)*uk[j]; //matrix multiplication A*u
            Au[i] = sum;
            resid[i] = pow((Au[i] - f[i]),2);   //residual ||Auk - f||
        }
        sum = 0.0;
        
        //finish computing residual ||Auk - f||
        for(int i = 0; i < N; i++)
            sum = sum + resid[i];
        r = sqrt(sum);
        
        if(k == 1)
            r0 = r;
        
        cout << k << "th iteration residual: " << r << endl;
        if(k > 5000)    //break if exceeds 5000 iterations
        {
            cout << "\n!!Iteration count exceeded 5000!!\n";
            break;
        }
    }
    cout << "\nTime to find approximate solution: " << t.toc() << endl;
    cout << "\nFinal residual after " << k << " iterations: ||Auk - f|| = " << r << endl << endl;
    free(A);
}


