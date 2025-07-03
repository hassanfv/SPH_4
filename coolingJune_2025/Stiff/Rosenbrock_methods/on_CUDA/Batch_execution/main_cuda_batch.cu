%%writefile stiff_test.cu
#include <iostream>
#include <cstdio>
#include "diff_jacob_batch_01.h"
#include "stiff_libs_hfv_batch_01.h"
#include <cuda_runtime.h>
#include <chrono>

using namespace std;

int main()
{
  const int N_part = 100000; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  const int N_con = 5000; //!!!!!!!!!!! number of concurrent access of the matrices a and dfdy !// Select efficiently. The higher the better !
  int *h_slot_status = new int[N_con];
  for (int i = 0; i < N_con; i++)
    h_slot_status[i] = 0;
  
  const int n = 3; // !!!!!!!!!!!!!!!!!! This is nvar !!!!!!!!!!!!!!!!! Also modify the value for vv in ludcmp function !!!!!!

  const float x1 = 0.0f; // Initial time !!!!!!!!!!!!!!!!!!!!!!!
	const float x2 = 50.0f; // Final time !!!!!!!!!!!!!!!!!!!!!!!!!
	float eps = 1e-4f;   // THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float htry = 2.9e-4f; // Step size !!!!!!!!!!!!!!!!!!!!!!!!
	float hmin = 1e-6f;  //!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

  float *h_y0 = new float[N_part * n];
  //!!!!!!!!!!!!! initial values !!!!!!!!!!!!!!!!!
  // These will be the initial ionization fractions of each SPH particle!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for (int i = 0; i < N_part; i++)
  {
    h_y0[i * n + 0] = 1.0f; // e.g. 1.0f will be replaced by nHII or nHeI of particle i !
    h_y0[i * n + 1] = 1.0f;
    h_y0[i * n + 2] = 0.0f;
  }
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  // Allocate device memory-- except "d_y0", all other device variables are only needed in the device and they do not need host counterparts !
  float *d_y0, *d_dydx, *d_yscal; // Note that the final result will also be copied back to d_y0 !
  float *d_a, *d_dfdy, *d_dfdx, *d_dysav, *d_err, *d_ysav; // a and dfdy are each a matrix !
  float *d_g1, *d_g2, *d_g3, *d_g4;
  int *d_indx, *d_slot_status;

  CHECK_CUDA(cudaMalloc(&d_slot_status, N_con * sizeof(int)));

  CHECK_CUDA(cudaMalloc(&d_y0, N_part * n * sizeof(float))); // Only y0 has N_part*n others have N_con*n.
  CHECK_CUDA(cudaMalloc(&d_dydx, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_yscal, N_con * n * sizeof(float)));

  CHECK_CUDA(cudaMalloc(&d_a, N_con * n * n * sizeof(float)));  // 3D !
  CHECK_CUDA(cudaMalloc(&d_dfdy, N_con * n * n * sizeof(float))); // 3D !
  CHECK_CUDA(cudaMalloc(&d_indx, N_con * n * sizeof(int)));
  CHECK_CUDA(cudaMalloc(&d_dfdx, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_dysav, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_err, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_ysav, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_g1, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_g2, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_g3, N_con * n * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_g4, N_con * n * sizeof(float)));

  // Copy to device
  CHECK_CUDA(cudaMemcpy(d_slot_status, h_slot_status, N_con * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_y0, h_y0, N_part * n * sizeof(float), cudaMemcpyHostToDevice));



  int nThreadsPerBlock = 256;
  //int nBlocksPerGrid = (N_part + nThreadsPerBlock - 1) / nThreadsPerBlock; // automate it !!
  
  auto start = chrono::high_resolution_clock::now();
  
  //************** Batch running of the Kernel ****************
  for (int i = 0; i < N_part; i += N_con)
  {
    int N_batch = min(N_con, N_part - i); // last batch might be smaller.
    
    odeint_d<<<(N_batch + nThreadsPerBlock - 1)/nThreadsPerBlock, nThreadsPerBlock>>>(d_y0 + i * n, x1, x2, eps, htry, hmin, n, d_dydx, d_yscal, d_a,
                                                                                      d_dfdy, d_indx, d_dfdx, d_dysav, d_err, d_ysav, d_g1, d_g2, d_g3,
                                                                                      d_g4, d_slot_status, N_batch); // n is nvar !
  }
  //***********************************************************
  cudaDeviceSynchronize();
  
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = end - start;
  cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;
  
  cudaMemcpy(h_y0, d_y0, N_part * n * sizeof(float), cudaMemcpyDeviceToHost);
  
  int i = 0; // i is an SPH particle index.
  
  float y0 = h_y0[i * n + 0];
  float y1 = h_y0[i * n + 1];
  float y2 = h_y0[i * n + 2];
  
  cout << y0 << ", " << y1 << ", " << y2 << endl;
  
}


