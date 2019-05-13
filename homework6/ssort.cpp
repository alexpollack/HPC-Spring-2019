// Parallel sample sort
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <algorithm>

int main( int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  //MPI_Init(&argc, &argv);

  int rank, p, s;// rootSamples, splitters, sdispls;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm comm = MPI_COMM_WORLD;
  
  s = p - 1; // number of samples

  int* rootSplitters = (int *)malloc(p*s*sizeof(int)); // root to store and sort all sample sets
  int* splitters = (int *)malloc(s*sizeof(int)); // splitters

  // Number of random numbers per processor (this should be increased
  // for actual tests or could be passed in through the command line
  int N = 100;

  int* vec = (int*)malloc(N*sizeof(int));
  int* sdispls = (int*)malloc(s*sizeof(int));
  int* rdispls = (int*)malloc(s*sizeof(int));

  // seed random number generator differently on every core
  srand((unsigned int) (rank + 393919));

  // fill vector with random integers
  for (int i = 0; i < N; ++i)
  {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  // sort locally
  std::sort(vec, vec+N);

  // sample p-1 entries from vector as the local splitters, i.e.,
  // every N/P-th entry of the sorted vector
  for(int i = 1; i < s; i++)
  {
    splitters[i] =  vec[(N/p)*i];
    //splitters[i] = vec[N/p * (i+1)];
  }

  // every process communicates the selected entries to the root
  // process; use for instance an MPI_Gather
  MPI_Gather(splitters, s, MPI_INT, rootSplitters, s, MPI_INT, 0, comm);

  // root process does a sort and picks (p-1) splitters (from the
  // p(p-1) received elements)
  if (rank == 0)
  {
    std::sort(rootSplitters, rootSplitters+s);
    for(int i = 1; i < s; i++)
      splitters[i-1] = rootSplitters[i*s - 1];
      //splitters[i] = rootSplitters[s*(i+1)];
    
    // root process broadcasts splitters to all other processes
  }
  MPI_Bcast(splitters, s, MPI_INT, 0, comm);

  // every process uses the obtained splitters to decide which
  // integers need to be sent to which other process (local bins).
  // Note that the vector is already locally sorted and so are the
  // splitters; therefore, we can use std::lower_bound function to
  // determine the bins efficiently.
  //
  // Hint: the MPI_Alltoallv exchange in the next step requires
  // send-counts and send-displacements to each process. Determining the
  // bins for an already sorted array just means to determine these
  // counts and displacements. For a splitter s[i], the corresponding
  // send-displacement for the message to process (i+1) is then given by,
  // sdispls[i+1] = std::lower_bound(vec, vec+N, s[i]) - vec;
  int * bcounts = (int*)malloc(p*sizeof(int));
  int * allct = (int*)malloc(p*sizeof(int));

  sdispls[0] = 0;
  int sum = 0, bucketSize;
  for(int i = 0; i < s; i++)
  {
    sdispls[i+1] = std::lower_bound(vec, vec+N, splitters[i]) - vec;
    bcounts[i+1] =  sdispls[i+1] - sdispls[i];
  }
  bcounts[0] = sdispls[1] - sdispls[0];

  // send and receive: first use an MPI_Alltoall to share with every
  // process how many integers it should expect, and then use
  // MPI_Alltoallv to exchange the data
  MPI_Alltoall(bcounts,1,MPI_INT,allct,1,MPI_INT,comm);
  rdispls[0] = 0;
  for(int i = 1 ; i < p; i++)
  {
    rdispls[i] = rdispls[i-1] + allct[i-1]; 
    bucketSize += allct[i];
  }
  MPI_Barrier(comm);
  MPI_Alltoallv(vec, bcounts, sdispls, MPI_INT, buckets,allct, rdispls, MPI_INT, comm);
  MPI_Barrier(comm);
  // do a local sort of the received data
  //std::sort(buckets, buckets+bucketSize);
  // every process writes its result to a file
/*
  { // Write output to a file
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd) {
      printf("Error opening file \n");
      return 1;
    }
    for(int i = 0; i < N; i++)
      fprintf(fd, "  %f\n", buckets[i]);

    fclose(fd);
  }*/


  // free allocated memory
  free(bucket_sizes);
  bucket_sizes = NULL;
  free(rootSplitters);
  rootSplitters= NULL;
  free(splitters);
  splitters= NULL;
  free(sdispls);
  sdispls= NULL;
  free(rdispls);
  rdispls= NULL;
  free(bcounts);
  bcounts= NULL;
  free(allct);
  allct= NULL;
  free(buckets);
  free(vec);
  vec= NULL;


  MPI_Finalize();
  return 0;
}














