// original
#include "utils.hpp"
// std
#include <stdio.h>
#include <time.h>
// lib
#include <cuda_runtime.h>

void initial_data(float *ip, int size) {
  // create random seed
  time_t t;
  srand((unsigned int) time(&t));

  for (int i = 0; i < size; i++)
    ip[i] = (float)(rand() & 0xFF) / 10.0f;
  return;
}

void sum_array_on_host(float* a, float* b, float* c, const int n)
{
  for (int i = 0; i < n; i++)
    c[i] = a[i] + b[i];
}

// this function is called by all thread parallely
__global__ void sum_array_on_gpu(float* a, float* b, float* c)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  c[i] = a[i] + b[i];
}

bool check_result(float* host_ref, float* gpu_ref, const int n)
{
  double epsilon = 1.0e-8;
  for (int i = 0; i < n; i++)
    if (abs(host_ref[i] - gpu_ref[i]) > epsilon)
      return false;
  return true;
}

void setup_device()
{
  int dev = 0;
  cudaDeviceProp device_prop;
  CHECK(cudaGetDeviceProperties(&device_prop, dev));
  std::cout << "Using Device : " << device_prop.name << std::endl;
  CHECK(cudaSetDevice(dev));
}

int main(int argc, char *argv[]) {
  setup_device();

  if (argc != 2)
    throw std::runtime_error("usage : binary_filename block_size");

  float *data_a, *data_b, *data_c;
  float *host_a, *host_b, *host_ref, *gpu_ref;
  int n_element = 1 << 24;
  std::cout << "Vector size : " << n_element << std::endl;
  size_t n_bytes = n_element * sizeof(float); 
  
  // allocate memory for host
  host_a   = (float *)malloc(n_bytes);
  host_b   = (float *)malloc(n_bytes);
  host_ref = (float *)malloc(n_bytes);
  gpu_ref  = (float *)malloc(n_bytes);

  // allocate memory for cuda device
  cudaMalloc((float**)&data_a, n_bytes);
  cudaMalloc((float**)&data_b, n_bytes);
  cudaMalloc((float**)&data_c, n_bytes);
  
  // create data
  initial_data(host_a, n_element);
  initial_data(host_b, n_element);

  // copy data to the gpu
  cudaMemcpy(data_a, host_a, n_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(data_b, host_b, n_bytes, cudaMemcpyHostToDevice);

  {
    scope_timer timer("sum_array_on_host");
    sum_array_on_host(host_a, host_b, host_ref, n_element);
  }

  // execute kernel
  int  block_size = atoi(argv[1]);
  dim3 block(block_size);
  dim3 grid((n_element + block.x - 1) / block.x);
  {
    scope_timer timer("sum_array_on_gpu");
    sum_array_on_gpu<<<grid, block>>>(data_a, data_b, data_c);

    CHECK(cudaDeviceSynchronize());
    // copy the result
    cudaMemcpy(gpu_ref, data_c, n_bytes, cudaMemcpyDeviceToHost);
  } 

  auto result = check_result(host_ref, gpu_ref, n_element);
  if (result)
    printf("correct answer.\n");
  else
    printf("wrong answer.\n");

  // free host memory
  free(host_a);
  free(host_b);
  free(host_ref);
  free(gpu_ref);

  // free cuda memory
  cudaFree(data_a);
  cudaFree(data_b);
  cudaFree(data_c);
}