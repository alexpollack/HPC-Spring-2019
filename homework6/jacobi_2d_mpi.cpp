#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq)
{
  int i, j;
  double tmp, gres = 0.0, lres = 0.0, sum = 0.0;;

  for (i = 1; i <= lN; i++)
  {
    for (i = 1; i <= lN; i++)
    {
      tmp = ((4.0*lu[i*lN + j] - lu[(j-1)+i*lN] - lu[j+(i-1)*lN] - lu[(j+1)+i*lN] - lu[j+(i+1)*lN]) * invhsq-1); // - 1 ??
      sum += tmp;
    }
    lres += sum * sum;
  }
  /* use allreduce for convenience; a reduce would also be sufficient */
  MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(gres);
}


int main(int argc, char * argv[])
{
  int mpirank, i, j, p, N, NN, lN, iter, max_iters, index;
  MPI_Status status, status1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /* get name of host running MPI process */
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

  //N = 4; 
  //max_iters = 10;
  //NN = N*N;

  /* compute number of unknowns handled by each process */
  //lN = N / p;
  lN = N / p;
  if ((N % p != 0) && mpirank == 0 ) 
  {
    printf("N: %d, local N: %d\n", N, lN);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  if ( mpirank == 0 ) 
  {
    printf("N: %d, local N: %d\n", N, lN);
    printf("Exiting. NN must be a multiple of p\n");
  }
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  double tt = MPI_Wtime();

  /* Allocation of vectors, including left/upper and right/lower ghost points */
  double * lu    = (double *) calloc(sizeof(double), (lN + 2)*(lN + 2));
  double * lunew = (double *) calloc(sizeof(double), (lN + 2)*(lN + 2));
  double * lutemp = (double *) calloc(sizeof(double), (lN + 2)*(lN + 2));

  double * sendright = (double *) calloc(sizeof(double), lN);
  double * sendleft = (double *) calloc(sizeof(double), lN);
  double * sendup = (double *) calloc(sizeof(double), lN);
  double * senddown = (double *) calloc(sizeof(double), lN);

  double * recright = (double *) calloc(sizeof(double), lN);
  double * recleft = (double *) calloc(sizeof(double), lN);
  double * recup = (double *) calloc(sizeof(double), lN);
  double * recdown = (double *) calloc(sizeof(double), lN);

  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  double gres, gres0, tol = 1e-5;

  /* initial residual */
  gres0 = compute_residual(lu, lN, invhsq);
  gres = gres0;

  // Note: matrices are stored in column major order; i.e. the array elements in
  // the (m x n) matrix C are stored in the sequence: {C_00, C_10, ..., C_m0,
  // C_01, C_11, ..., C_m1, C_02, ..., C_0n, C_1n, ..., C_mn}

  for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) 
  {

    /* Jacobi step for local points */
    for(i = 1; i <= lN; i++)
    {
      for(j = 1; j <= lN; j++)
        lunew[i*lN + j]  = 0.25 * (hsq+lu[(j-1)+i*lN]+lu[j+(i-1)*lN]+lu[(j+1)+i*lN]+lu[j+(i+1)*lN]);
        //lunew[i]  = 0.25 * (hsq*f[i+j*lN]+lu[(i-1)+j*lN]+lu[i+(j-1)*lN]+lu[(i+1)+j*lN]+lu[i+(j+1)*lN]);
    }
    //set up values to send
    for(int k = 0; k < lN; k++) //rows
    {
      sendup[k] = lunew[k]; // values to send up is just first row
      senddown[k] = lunew[(lN - 1)*lN + k]; // values to send down is just bottom (lN'th) row
      sendleft[k] = lunew[k*lN]; // values to send left is just right most column
      sendright[k] = lunew[(k+1)*lN - 1]; // values to send right is just left most column
    }
    /* communicate ghost values */
    if (mpirank < (p - lN))//(mpirank < p - 1) // (lN-1)*lN
    {
      //set up values to send

      /* If not the top row of processes, send top row, rec row above */
      MPI_Send(&(sendup), lN, MPI_DOUBLE, mpirank+lN, 124, MPI_COMM_WORLD);
      MPI_Recv((recup), lN, MPI_DOUBLE, mpirank+lN, 123, MPI_COMM_WORLD, &status);
    }
    if (mpirank >= lN) 
    {
      /* If not the first row of process, send bottom row, rec row below */
      MPI_Send(&(senddown), lN, MPI_DOUBLE, mpirank-lN, 123, MPI_COMM_WORLD);
      MPI_Recv((recdown), lN, MPI_DOUBLE, mpirank-lN, 124, MPI_COMM_WORLD, &status);
    }
    if (mpirank % lN != 0)
    {
      // if not in left most column of processors, send left most column, rec col to right
      MPI_Send(&(sendleft), lN, MPI_DOUBLE, mpirank-1, 122, MPI_COMM_WORLD);
      MPI_Recv((recleft), lN, MPI_DOUBLE, mpirank-1, 121, MPI_COMM_WORLD, &status);
    }
    if ((mpirank+1) % lN != 0)
    { 
      // if not in right most column of processors, send right most column, rec col to right
      MPI_Send(&(sendright), lN, MPI_DOUBLE, mpirank+1, 121, MPI_COMM_WORLD);
      MPI_Recv((recright), lN, MPI_DOUBLE, mpirank+1, 122, MPI_COMM_WORLD, &status);
    }


    /* copy newu to u using pointer flipping
     * i.e. lutemp = lu; lu = lunew; lunew = lutemp */
    for(i = 1; i <= lN; i++)
    {
      for(j = 1; j <= lN; j++)
      {
        lutemp[i*lN + j] = lu[i*lN + j];
        lu[i*lN + j] = lunew[i*lN + j];
        lunew[i*lN + j] = lutemp[i*lN + j];
      }
    }
    if (0 == (iter % 10)) 
    {
      gres = compute_residual(lu, lN, invhsq);
      if (0 == mpirank) 
      {
        printf("Iter %d: Residual: %g\n", iter, gres);
      }
    }
  }

  if (0 == mpirank) 
  {
    printf("Iter %d: Residual: %g\n", iter, gres);
  }

  /* Clean up */
  free(lu);
  free(lunew);
  free(lutemp);

  free(sendright);
  free(sendleft);
  free(sendup);
  free(senddown);

  free(recright);
  free(recleft);
  free(recup);
  free(recdown);

  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  double elapsed = MPI_Wtime() - tt;
  if (0 == mpirank) 
  {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }
  MPI_Finalize();
  return 0;
}



