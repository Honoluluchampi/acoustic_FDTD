#include <chrono>
#include <iostream>

#include <cuda_runtime.h>

// error checking macro
#define CHECK(call)\
{\
  const cudaError_t error = call;\
  if (error != cudaSuccess) {\
    printf("Error:%s:%d, ", __FILE__, __LINE__);\
    printf("code:%d, reason:%s\n", error, cudaGetErrorString(error));\
    exit(1);\
  }\
}

struct scope_timer {
  // start timer by ctor
  scope_timer(const std::string& _entry = "") {
    entry = _entry;
    start = std::chrono::system_clock::now();
  }

  // stop and output elapsed time by dtor
  ~scope_timer() {
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "\x1b[35m" << "SCOPE TIMER : " << entry << std::endl;
    std::cout << "\t elapsed time : " << elapsed << " ms" << "\x1b[0m" << std::endl;
  }

  std::string entry;
  std::chrono::system_clock::time_point start;
};