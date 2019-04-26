// make clean
// make
// mpirun -np 2 int_ring

#include <stdio.h>
#include <cstdlib>
#include <mpi.h>
#include <iostream>


int main(int argc, char **argv)
{
  int world_rank, world_size, N = 100000; // N = number of repeats
  long Nsize = 1000000; // 2MBytes, will b int, from pingpong code.

  char* data = (char*) malloc(Nsize);
  for (long i = 0; i < Nsize; i++) data[i] = 1;

  MPI_Init(&argc , &argv);
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_rank == 0) {
    printf("\nNsize = %d\n", Nsize);
    printf("\nN = %d\n", N);
  }

  double tt = MPI_Wtime();
  for(int i = 0; i < N; i++)
  {
    if (world_rank == 0) 
    {
      MPI_Send( &data, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv( &data, 1 , MPI_INT, world_size-1, 0, MPI_COMM_WORLD, &status);
    }
    else
    {
      MPI_Recv( &data, 1 , MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, &status);
      if ( world_rank < world_size - 1) {
        MPI_Send( &data, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
      }
      if ( world_rank == world_size - 1) 
        MPI_Send( &data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
    }
    /*if (world_rank == 0)
      printf("Process %d got %d from %d\n", world_rank, data, world_size-1);
    else
      printf("Process %d got %d from %d\n", world_rank, data, world_rank-1);*/
  }

  
  tt = MPI_Wtime() - tt;
  if (!world_rank) printf("\nRing bandwidth: %e GB/s\n", (Nsize*N)/tt/1e9);
  MPI_Finalize();  

  if (world_rank == 0)
    printf("\nDone.\n");
  return 0;
}






