//
//  main.cpp
//  hpc_hw1
//  question 3: solving the Laplace equation by Jacobi Method
//  Created by Alex Pollack on 2/1/19.
//  Copyright © 2019 Alex Pollack. All rights reserved.
//  g++ -std=c++11 -O3 jacobi2D-omp.cpp && ./a.out

#include <iostream>
#include <cmath>
#include <omp.h>
#include <math.h>
#include <vector>
#include "utils.h"

using namespace std;

int main()
{
    int N;
    cout << "Enter dimension N: ";  //input dimension size
    cin >> N;
    double h = 1.0/(N+1.0);
    //double *f = (double*)malloc(2*N * 2*N * sizeof(double));
    //double *u = (double*)malloc(2*N * 2*N * sizeof(double));
    double f[N][N], u[N][N];
    for(int i = 0; i <N; i++)//f(x) = 1
    {
        for(int j = 0; j < N; j++)
        {
            f[i][j] = 1.0;
            u[i][j] = 0.0;  //initialize
        }
    }
    
    cout << "Jacobi Method:\n";
    int k = 0;
    Timer t;
    
    double r = 1.0, r0= 1.0, sum = 0.0, uk1[N][N], resid[N];
    for(int i = 1; i < N; i ++)
    {
        uk1[0][i]=0.0;
        uk1[i][0]=0.0;
        uk1[N][i]=1.0;
        uk1[i][N]=1.0;
    }//Boundary conditions
    //begin Jacobi Method
    t.tic();
    int i,j;
    int chunk = 1000;
    while(r/r0 > 1e-6)  //end method when residual has been reduced by 1e6 from first try
    {
        omp_set_num_threads(4);
        #pragma omp parallel for default(none) shared(uk,h,u,f,i,j,chunk) //reduction(+:C_ij)
        {
        for(i = 1; i < N; i ++)
        {
            for(j=1; j < N; j++)
            {
                uk1[i][j] = 0.25 * (h*h*f[i][j] + u[i-1][j] +u[i][j-1]+u[i+1][j]+u[i][j+1]);
            }
        }
        }
        for(i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            u[i][j] = uk1[i][j]; //u^k gets new u^k+1 value for future computations
        }
        sum = 0.0;
        for (i = 0; i < N; i++)
        {
            for(j=0;j<N;j++)
                sum = sum + pow((u[i][j] - f[i][j]),2);   //residual ||uk - f||
            resid[i] = sum;
            sum = 0.0;
        }
        sum = 0.0;
        //finish computing residual ||Auk - f||
        for(int i = 0; i < N; i++)
            sum = sum + resid[i];
        r = sqrt(sum);
        if(k == 1)
            r0 = r;
        cout << k << "th iteration residual: " << r/r0 << endl;
        if(k > 5000)    //break if exceeds 5000 iterations
        {
            cout << "\n!!Iteration count exceeded 5000!!\n";
            break;
        }
        k++;    //iteration count
    }
    cout << "\nTime to find approximate solution: " << t.toc() << endl;
    cout << "\nFinal residual after " << k << " iterations: ||Auk - f|| = " << r/r0 << endl << endl;
    //free(A);
    
    return 0;
}

