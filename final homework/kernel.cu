#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "cuda_runtime_api.h"

#include<fstream>
#include<sstream>

const int col = 254, elinenum = 53; //����������Ԫ����
const int bytenum = (col - 1) / 32 + 1;   //ÿ��ʵ���е�byte��������

class bitmatrix {
public:
	int mycol;    //����
	int *mybyte;
	bitmatrix() {    //��ʼ��
		mycol = -1;
		mybyte = new int[bytenum];
		for (int i = 0; i < bytenum; i++)
			mybyte[i] = 0;
	}
	void insert(int x) { //���ݶ���
		if (mycol == -1)mycol = x;
		int a = x / 32, b = x % 32;
		mybyte[a] |= (1 << b);
	}
};


bitmatrix *eliminer = new bitmatrix[col], *eline = new bitmatrix[elinenum];
void readdata() {
	using namespace std;
	ifstream ifs;
	ifs.open("D:\\VS��Ŀ\\cuda\\cuda\\eliminer1.txt");  //��Ԫ��
	string temp;
	while (getline(ifs, temp)) {
		istringstream ss(temp);
		int x;
		int trow = 0;
		while (ss >> x) {
			if (!trow)trow = x;    //��һ������Ԫ�ش����к�
			eliminer[trow].insert(x);
		}
	}
	ifs.close();
	ifstream ifs2;
	ifs2.open("D:\\VS��Ŀ\\cuda\\cuda\\eline1.txt");     //����Ԫ��,���뷽ʽ����Ԫ�Ӳ�ͬ
	int trow = 0;
	while (getline(ifs2, temp)) {
		istringstream ss(temp);
		int x;
		while (ss >> x) {
			eline[trow].insert(x);
		}
		trow++;
	}
	ifs2.close();
}
/*void dowork() {  //������Ԫ--����Ԫ��->��Ԫ��
	for (int i = 0; i < elinenum; i++) {
		while (!eline[i].isnull()) {  //ֻҪ����Ԫ�зǿգ�ѭ������
			int tcol = eline[i].mycol;  //����Ԫ�е�����
			if (!eliminer[tcol].isnull())    //������ڶ�Ӧ��Ԫ��
				eline[i].doxor(eliminer[tcol]);
			else {
				eliminer[tcol] = eline[i];    //���ڱ���Ԫ������Ϊ��Ԫ�Ӻ󲻲��������������ֱ����=��ǳ����
				break;
			}
		}
	}
}*/

__global__ void dowork1(int **gelinebyte,int *gelinecol,int **geliminerbyte,int *geliminercol,int i,int elinenum1,int bytenum1) {  //������Ԫ--��Ԫ��->����Ԫ��
	if (geliminercol[i] == -1)return;
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	while (bid < elinenum1) {
		if (gelinecol[bid] == i) {
			int temp = tid;
			while (temp < bytenum1) {
				gelinebyte[bid][temp] ^= geliminerbyte[i][temp];
				temp += blockDim.x;
			}
			__syncthreads();
			if (tid == 0) {
				bool f = 0;
				for (int k = bytenum1 - 1; k >= 0&&!f; k--)
					for (int j = 31; j >= 0&&!f; j--)
						if ((gelinebyte[bid][k] & (1 << j)) != 0) {
							gelinecol[bid] = k * 32 + j;
							f = 1;
						}
				if(!f)gelinecol[bid] = -1;
			}
		}
		bid += gridDim.x;
	}
}
__global__ void dowork2(int **gelinebyte, int *gelinecol, int **geliminerbyte, int *geliminercol, int i, int elinenum1) {
	if (blockIdx.x==0&&threadIdx.x == 0) {
		for (int j = 0; j < elinenum1; j++)
			if (gelinecol[j] == i) {
				geliminerbyte[i] = gelinebyte[j];
				geliminercol[i] = gelinecol[j];
				return;
			}
	}
}

void printres(int** celinebyte,int *celinecol) { //��ӡ���
	for (int i = 0; i < elinenum; i++) {
		if (celinecol[i]==-1) { puts(""); continue; }   //���е��������
		for (int j = bytenum - 1; j >= 0; j--) {
			for (int k = 31; k >= 0; k--)
				if ((celinebyte[i][j] & (1 << k)) != 0) {     //һ��������˰�Сʱ�����ǵ���λΪ1ʱ>>�����ڳ�����
					printf("%d ", j * 32 + k);
				}
		}
		puts("");
	}
}
int main() {
	readdata();
	int** geliminerbyte;
	int* geliminercol;
	int** gelinebyte;
	int* gelinecol;
	int** celinebyte = new int*[elinenum];
	int *celinecol=new int[elinenum];
	int** celiminerbyte = new int*[col];
	int *celiminercol=new int[col];
	for (int i = 0; i < elinenum; i++) {
		int* host_1d = new int[bytenum];
		for (int j = 0; j < bytenum; j++)
			host_1d[j] = eline[i].mybyte[j];
		celinecol[i] = eline[i].mycol;
		int* dev_1d;
		cudaMalloc((void**)&dev_1d, sizeof(int)*bytenum);// ��ʱ dev_ld ָ�� һƬ�Դ�ռ�
		cudaMemcpy(dev_1d, host_1d, sizeof(int)*bytenum, cudaMemcpyHostToDevice);
		celinebyte[i] = dev_1d;
	}
	cudaMalloc((void**)&gelinebyte, sizeof(int*)*elinenum); // ����һ���Դ��ϵĶ���ָ��
	cudaMemcpy(gelinebyte, celinebyte, sizeof(int*)*elinenum, cudaMemcpyHostToDevice);
	cudaMalloc(&gelinecol, sizeof(int)*elinenum);
	cudaMemcpy(gelinecol, celinecol, sizeof(int)*elinenum, cudaMemcpyHostToDevice);
	
	for (int i = 0; i < col; i++) {
		int* host_1d = new int[bytenum];
		for (int j = 0; j < bytenum; j++)
			host_1d[j] = eliminer[i].mybyte[j];
		celiminercol[i] = eliminer[i].mycol;
		int* dev_1d;
		cudaMalloc((void**)&dev_1d, sizeof(int)*bytenum);// ��ʱ dev_ld ָ�� һƬ�Դ�ռ�
		cudaMemcpy(dev_1d, host_1d, sizeof(int)*bytenum, cudaMemcpyHostToDevice);
		celiminerbyte[i] = dev_1d;
	}
	cudaMalloc((void**)&geliminerbyte, sizeof(int*)*col); // ����һ���Դ��ϵĶ���ָ��
	cudaMemcpy(geliminerbyte, celiminerbyte, sizeof(int*)*col, cudaMemcpyHostToDevice);
	cudaMalloc(&geliminercol, sizeof(int)*col);
	cudaMemcpy(geliminercol, celiminercol, sizeof(int)*col, cudaMemcpyHostToDevice);
	cudaEvent_t start, stop;//��ʱ��
	float etime = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);//��ʼ��ʱ
	for (int i = col - 1; i >= 0; i--) {
		if (eliminer[i].mycol == -1)
			dowork2 << <1, 1>> > (gelinebyte,gelinecol,geliminerbyte,geliminercol, i, elinenum);
		cudaDeviceSynchronize();
		dowork1 << <1024, 1024 >> > (gelinebyte, gelinecol,geliminerbyte,geliminercol, i, elinenum, bytenum);
		cudaDeviceSynchronize();
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);//ֹͣ��ʱ
	cudaEventElapsedTime(&etime, start, stop);
	printf("GPU_LU:%f ms\n", etime);

	//cudaMemcpy(celinebyte, gelinebyte, sizeof(int*)*elinenum, cudaMemcpyDeviceToHost);
	cudaMemcpy(celinecol, gelinecol, sizeof(int)*elinenum, cudaMemcpyDeviceToHost);
	//for (int i = 0; i < elinenum; i++)
	//	printf("%d\n", celinecol[i]);
	//printres(celinebyte,celinecol);
	cudaFree((void*)gelinebyte);
	cudaFree(gelinecol);
	cudaFree((void*)geliminerbyte);
	cudaFree(geliminercol);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	return 0;
}