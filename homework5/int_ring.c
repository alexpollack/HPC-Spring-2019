// make clean
// make
// mpirun -np 2 int_ring

#include <stdio.h>
#include <cstdlib>
#include <mpi.h>
#include <iostream>


int main(int argc, char **argv)
{
  int world_rank, world_size, data = 0, N = 10000; // N = number of repeats

  MPI_Init(&argc , &argv);
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_rank == 0)
    printf("\nN = %d\n", N);

  double tt = MPI_Wtime();
  for(int i = 0; i < N; i++)
  {
    if (world_rank == 0) 
    {
      MPI_Send( &data, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv( &data, 1 , MPI_INT, world_size-1, 0, MPI_COMM_WORLD, &status);
      data += world_size - 1;
    }
    else
    {
      MPI_Recv( &data, 1 , MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, &status);
      data += world_rank - 1;
      if ( world_rank < world_size - 1) {
        MPI_Send( &data, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
      }
      if ( world_rank == world_size - 1) 
        MPI_Send( &data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
    }
    // Print out what gets what from where (only ideal if # of repeats isn't large)
    /*if (world_rank == 0)
      printf("Process %d got %d from %d\n", world_rank, data, world_size-1);
    else
      printf("Process %d got %d from %d\n", world_rank, data, world_rank-1);*/
  }
  tt = MPI_Wtime() - tt;
  if (!world_rank) printf("\nRing latency: %e ms\n", tt/N * 1000);
  MPI_Finalize();  

  if (world_rank == 0)
    printf("\nFinal: Process 0 ended with %d\n\n", data); // Print what the root process 0 ended up with as thats the end of the ring
  return 0;
}








