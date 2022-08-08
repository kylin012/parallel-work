#include "mpi.h" 
#include <iostream> 
#include <emmintrin.h>
#include <immintrin.h>
#include<algorithm>
using namespace std;
const int n = 1000,mpisize = 16;
float m1[n][n], time1;
bool parallel = 1;
void printres(int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			std::cout << m1[i][j] << " ";
		std::cout << std::endl;
	}
}
int main(int argc, char* argv[])
{
	int rank;
	double st, ed;
	MPI_Status status;
	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//��ȡ��ǰ���̺�
	if (rank == 0) {	//0�Ž���ʵ�ֳ�ʼ��
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				m1[i][j] = 0;
			m1[i][i] = 1.0;
			for (int j = i + 1; j < n; j++)
				m1[i][j] = rand() % 1000 + 1;
		}
		for (int k = 0; k < n; k++)
			for (int i = k + 1; i < n; i++)
				for (int j = 0; j < n; j++)
					m1[i][j] = int((m1[i][j] + m1[k][j])) % 1000 + 1.0;
		for (int j = 1; j < mpisize; j++) {
			for(int i= j;i<n;i+=mpisize-1)
			MPI_Send(&m1[i][0], n, MPI_FLOAT, j, i, MPI_COMM_WORLD);
		}
	}
	else 
		for(int i=rank;i<n;i+=mpisize-1)
			MPI_Recv(&m1[i][0], n, MPI_FLOAT,0, i, MPI_COMM_WORLD, &status);	
	MPI_Barrier(MPI_COMM_WORLD);	//������ͬ��
	st = MPI_Wtime();	 //��ʱ
	for (int k = 0; k < n; k++) {	
		if (rank == 0) {	//0�Ž��̸����������
			for (int j = k + 1; j < n; j++)
				m1[k][j] /= m1[k][k];
			m1[k][k] = 1.0;
			for (int j = 1; j < mpisize; j++)
				MPI_Send(&m1[k][0], n, MPI_FLOAT, j,n+k + 1, MPI_COMM_WORLD);
		}
		else
			MPI_Recv(&m1[k][0], n, MPI_FLOAT, 0, n+k + 1, MPI_COMM_WORLD, &status);
		if (rank != 0) {
			int r2 = rank;
			while (r2 < k + 1)r2 += mpisize - 1;
			for (int i = r2; i < n; i += mpisize - 1) {//����k+1��֮��ĸ����̲������м���
				for (int j = k + 1; j < n; j++)
					m1[i][j] -= m1[i][k] * m1[k][j];
				m1[i][k] = 0;
				if (i == k + 1) 
					MPI_Send(&m1[i][0], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	//��������㲥��0�Ž���
			}
		}
		else if(k<n-1)
			MPI_Recv(&m1[k + 1][0], n, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);	//������ͬ��
	ed = MPI_Wtime();	 //��ʱ
	MPI_Finalize();
	if (rank == 0) {	//ֻ��0�Ž����������ս��
		printf("cost time:%.4lf\n", ed - st);
		//printres(n);
	}
	return 0;
}

