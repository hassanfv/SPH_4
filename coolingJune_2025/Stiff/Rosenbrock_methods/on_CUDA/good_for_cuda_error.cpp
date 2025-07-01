#include <iostream>
#include <cuda_runtime.h>
#include <cstdio>

#define CHECK_CUDA_ERROR(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << ": " \
                  << cudaGetErrorString(err) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

using namespace std;

__global__ void test_h(int *a)
{
    *a += 1;
}

int main()
{
    int a = 0;
    int *d_a;

    CHECK_CUDA_ERROR(cudaMalloc(&d_a, sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpy(d_a, &a, sizeof(int), cudaMemcpyHostToDevice));

    test_h<<<1, 1>>>(d_a);
    
    // Check for kernel launch errors
    CHECK_CUDA_ERROR(cudaGetLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    CHECK_CUDA_ERROR(cudaMemcpy(&a, d_a, sizeof(int), cudaMemcpyDeviceToHost));

    cout << "a = " << a << endl;

    cudaFree(d_a);
    return 0;
}
