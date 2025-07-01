



N_con *


int main()
{

  const float x1 = 0.0f; // Initial time !!!!!!!!!!!!!!!!!!!!!!!
	const float x2 = 50.0f; // Final time !!!!!!!!!!!!!!!!!!!!!!!!!
	float eps = 1e-4f;   // THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float htry = 2.9e-4f; // Step size !!!!!!!!!!!!!!!!!!!!!!!!
	float hmin = 1e-6f;  //!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

  float y_host[N] = {1.0f, 0.0f, 0.0f, 0.0f};  // initial values
  float dydx_host[N], yscal_host[N], hdid_host, hnext_host;


  // Allocate device memory
  float *d_y, *d_dydx, *d_x, *d_yscal, *d_hdid, *d_hnext;
  float *d_a, *d_dfdy, *d_dfdx, *d_dysav, *d_err, *d_ysav;
  float *d_g1, *d_g2, *d_g3, *d_g4;
  int *d_indx;

  CHECK_CUDA(cudaMalloc(&d_x, N_con * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_hdid, N_con * sizeof(float)));
  CHECK_CUDA(cudaMalloc(&d_hnext, N_con * sizeof(float)));

  CHECK_CUDA(cudaMalloc(&d_y, N_con * n * sizeof(float)));
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

  // Copy inputs to device
  CHECK_CUDA(cudaMemcpy(d_x, &h_x, N_con * sizeof(float), cudaMemcpyHostToDevice));
  
  CHECK_CUDA(cudaMemcpy(d_y, h_y, N_con * n * sizeof(float), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_dydx, h_dydx, N_con * n * sizeof(float), cudaMemcpyHostToDevice));
  CHECK_CUDA(cudaMemcpy(d_yscal, h_yscal, N_con * n * sizeof(float), cudaMemcpyHostToDevice));
  
//!!!!!!!!!!!!!!! TRY TO GET RID of h_x and d_x

}
