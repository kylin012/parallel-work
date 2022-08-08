#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "cuda_runtime_api.h"
#include<stdlib.h>
const int n = 1000;
float m1[n*n];
void init(int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			m1[i*n + j] = 0;
		m1[i*n + i] = 1.0;
		for (int j = i + 1; j < n; j++)
			m1[i*n + j] = rand() % 1000 + 1;
	}
	for (int k = 0; k < n; k++)
		for (int i = k + 1; i < n; i++)
			for (int j = 0; j < n; j++)
				m1[i*n + j] = int((m1[i*n + j] + m1[k*n + j])) % 1000 + 1.0;
}
__global__ void division_kernel(float* data, int k, int N) {	//除法核函数
	int tid = blockDim.x * blockIdx.x + threadIdx.x;//计算线程索引
	if (tid > k&&tid < n)
		data[k*N + tid] /= data[k*N + k];
	return;
}
__global__ void eliminate_kernel(float* data, int k, int N) {	//消去核函数
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	if (tx == 0)  data[k*N + k] = 1.0;//对角线元素设为 1 
	int row = k + 1 + blockIdx.x;//每个块负责一行
	while (row > k&&row < N) {
		int tid = threadIdx.x;
		while (k + 1 + tid < N) {
			int col = k + 1 + tid;
			float temp_1 = data[(row*N) + col];
			float temp_2 = data[(row*N) + k];
			float temp_3 = data[k*N + col];
			data[(row*N) + col] = temp_1 - temp_2 * temp_3;
			tid = tid + blockDim.x;
		}
		__syncthreads();//块内同步
		if (threadIdx.x == 0) {
			data[row * N + k] = 0;
		}
		row += gridDim.x;
	}
	return;
}

int main() {
	init(n);
	float* gdata;
	int size = n * n * sizeof(float);
	cudaError_t ret;

	ret = cudaMalloc(&gdata, size);	//分配gpu空间
	if (ret != cudaSuccess) {
		printf("cudaMalloc gpudata failed!\n");
	}

	ret = cudaMemcpy(gdata, m1, size, cudaMemcpyHostToDevice);//将数据传输至 GPU 端

	if (ret != cudaSuccess) {
		printf("cudaMemcpyHostToDevice failed!\n");
	}
	dim3 grid(1024, 1);//线程块
	dim3 block(1024, 1);//线程网格
	cudaEvent_t start, stop;//计时器
	float etime = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);//开始计时

	for (int k = 0; k < n; k++) {
		division_kernel << <grid, block >> > (gdata, k, n);//负责除法任务的核函数
		cudaDeviceSynchronize();//CPU 与 GPU 之间的同步函数
		ret = cudaGetLastError();
		if (ret != cudaSuccess) {
			printf("division_kernel failed, %s\n", cudaGetErrorString(ret));
		}
		eliminate_kernel << <grid, block >> > (gdata, k, n);//负责消去任务的核函数
		cudaDeviceSynchronize();
		ret = cudaGetLastError();
		if (ret != cudaSuccess) {
			printf("eliminate_kernel failed, %s\n", cudaGetErrorString(ret));
		}
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);//停止计时
	cudaEventElapsedTime(&etime, start, stop);
	printf("GPU_LU:%f ms\n", etime);

	ret = cudaMemcpy(m1, gdata, size, cudaMemcpyDeviceToHost);//将数据传回 CPU 端
	/*for(int i=0;i<n;i++){
		for (int j = 0; j < n; j++)
			printf("%.2f ", m1[i*n+j]);
		puts("");
	}*/
	if (ret != cudaSuccess) {
		printf("cudaMemcpyDeviceToHost failed!\n");
	}
	cudaFree(gdata);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
}