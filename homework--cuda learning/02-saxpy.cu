#include <stdio.h>

#define N 2048 * 2048 // Number of elements in each vector

/*
 * Optimize this already-accelerated codebase. Work iteratively,
 * and use nsys to support your work.
 *
 * Aim to profile `saxpy` (without modifying `N`) running under
 * 20us.
 *
 * Some bugs have been placed in this codebase for your edification.
 */

__global__ void saxpy(int * a, int * b, int * c)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
   int gridnum=blockDim.x*gridDim.x;
    for(;tid<N;tid+=gridnum)
        c[tid] = 2 * a[tid] + b[tid];
}

int main()
{
    int *a, *b, *c;
    int deviceId;
    int numberofSMs;
    cudaGetDevice(&deviceId);
    cudaDeviceGetAttribute(&numberofSMs,cudaDevAttrMultiProcessorCount, deviceId);

    size_t size = N * sizeof (int); // The total number of bytes per vector

    cudaMallocManaged(&a, size);
    cudaMallocManaged(&b, size);
    cudaMallocManaged(&c, size);

     // Initialize memory
  for( int i = 0; i < N; ++i )
  {
    a[i] = 2;
    b[i] = 1;
    c[i] = 0;
  }
    

    int threads_per_block = 1024;
    int number_of_blocks = numberofSMs*32;

   cudaMemPrefetchAsync(a, size, deviceId);
  cudaMemPrefetchAsync(b, size, deviceId);
  cudaMemPrefetchAsync(c, size, deviceId);
     // Initialize memory
     
     
    saxpy <<< number_of_blocks, threads_per_block >>> ( a, b, c );
    cudaDeviceSynchronize();

    // Print out the first and last 5 values of c for a quality check
    for( int i = 0; i < 5; ++i )
        printf("c[%d] = %d, ", i, c[i]);
    printf ("\n");
    for( int i = N-5; i < N; ++i )
        printf("c[%d] = %d, ", i, c[i]);
    printf ("\n");

    cudaFree( a ); cudaFree( b ); cudaFree( c );
}
