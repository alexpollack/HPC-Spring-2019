#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = 0;
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i-1];
  }
}

void scan_omp(long* prefix_sum, const long* A, long n) {
    
    prefix_sum[0]  = 0;
    int p = 4;//, thread; sum = 0
    long* partial;
    int nthr;
    omp_set_num_threads(4);
    
    #pragma omp parallel //private (thread)
    {
        int i;
        //printf("number threads %d\n",omp_get_num_threads());
        #pragma omp single
        {
            nthr = omp_get_num_threads();
            partial = (long*) malloc(p* sizeof(long));
            partial[0] = 0;
        }
        int thread = omp_get_thread_num();
        int sum = 0;
        //printf("Hello from thread: %d\n",thread);
        
        #pragma omp for schedule(static) //nowait //private(i, sum)
        //for(i = thread*(n/p)+1; i < (thread+1)*(n/p); i++)// +1? //for(i = 0; i < n; i++)
        for( i = 1; i<n; i++)
        {
            //partial[thread] += A[i-1];
            sum += A[i-1];//A[i-1];
            //partial[thread] += A[i-1];
            prefix_sum[i] = sum;
            //printf("i: %d\n", thread, i);
        }
        //printf("\nThread %d has sum %d\n", thread, sum);
        partial[thread] = sum;
        #pragma omp barrier
        
        int offset = 0;//,j;
        for(i=1;i<thread+1;i++)//for(i = 0; i < p; i++)
            offset += partial[i-1];
        
        //printf("offset %d from thread: %d\n",offset,thread);
        #pragma omp for schedule(static) //private(j)
        //for(j = thread*(n/p) + 1; j < (thread+1)*(n/p); j++)
        for( i = 1; i<n; i++)
            prefix_sum[i] += offset; //partial[thread];
        
    }
    free(partial);
    
  // TODO: implement multi-threaded OpenMP scan
    
}

int main() {
    long N = 100000000;
  long* A = (long*) malloc(N * sizeof(long));
  long* B0 = (long*) malloc(N * sizeof(long));
  long* B1 = (long*) malloc(N * sizeof(long));
  for (long i = 0; i < N; i++) A[i] = rand();
  //  for (long i = 0; i < N; i++) A[i] = i;
    
    /*printf("A:\n");
    for(int i = 0; i< N; i++)
        printf("%ld\n", A[i]);
    printf("\n");*/
    
  double tt = omp_get_wtime();
  scan_seq(B0, A, N);
  printf("sequential-scan = %fs\n", omp_get_wtime() - tt);
    
    /*printf("B0:\n");
    for(int i = 0; i< N; i++)
        printf("%ld\n", B0[i]);
    printf("\n");*/
    
  tt = omp_get_wtime();
  scan_omp(B1, A, N);
  printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);
    
    /*printf("B1:\n");
    for(int i = 0; i< N; i++)
        printf("%ld\n", B1[i]);
    printf("\n");*/
    
  long err = 0;
  for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
  printf("error = %ld\n", err);

  free(A);
  free(B0);
  free(B1);
  return 0;
}
