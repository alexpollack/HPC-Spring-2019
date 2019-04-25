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
#include <omp.h>
#include <stdio.h>
using namespace std;

void GSM(double uk[], double f[], int N);

int main()
{
    #ifdef _OPENMP
        omp_set_num_threads(4);
    #endif
    int N=10;
    //cout << "Enter dimension N: ";  //input dimension size
    //cin >> N;
    double h = 1/(N+1);
    double f[N+1][N+1], u[N+1][N+1];
    for(int i = 0; i <=N; i++)//f(x) = 1
    {
        for(int j = 0; j <= N; j++)
        {
            f[i][j] = 1.0;
            u[i][j] = 0.0;  //initialize
        }
    }
    cout << "Gauss-Seidel Method:\n";
    int k = 0;
    Timer t;
    double r = 1.0, r0= 1.0, sum = 0.0, uk1[N+1][N+1], resid[N+1];//,ublack[N][N],ured[N][N];
    for(int i = 1; i <= N; i ++)
    {
        uk1[0][i]=0.0;
        uk1[i][0]=0.0;
        uk1[N][i]=1.0;
        uk1[i][N]=1.0;
    }//Boundary conditions



    t.tic();
    int i,j;
    while(r/r0 > 1e-6)  //end method when residual has been reduced by 1e6 from first try
    {

    //#ifdef _OPENMP
        //reduction(+:C_ij)
    //#endif
        #pragma omp parallel 
            #pragma omp for 
            for(i = 1; i <= N; i ++)
            {
                for(j=1; j <= N; j++)
                {
                    if(i+j % 2 == 0)  //i+j is even, do red
                        uk1[i][j] = 0.25 * (h*h*f[i][j] + u[i-1][j] +u[i][j-1]+u[i+1][j]+u[i][j+1]);
                    else if(i+j % 2 == 1) //i+j is odd, do black
                        uk1[i][j] = 0.25 * (h*h*f[i][j] + uk1[i-1][j] +uk1[i][j-1]+uk1[i+1][j]+uk1[i][j+1]);
                }
            }
        
        for(i = 0; i <=N; i++)
        {
            for(int j = 0; j <=N; j++)
                u[i][j] = uk1[i][j]; //u^k gets new u^k+1 value for future computations
        }
        sum = 0.0;
        for (i = 0; i <=N; i++)
        {
            for(j=0;j<=N;j++)
                sum = sum + pow((u[i][j] - f[i][j]),2);   //residual ||uk - f||
            resid[i] = sum;
            sum = 0.0;
        }
        sum = 0.0;
        //finish computing residual ||Auk - f||
        for(int i = 0; i <=N; i++)
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
        k++;
    }
    return 0;
}
    

