#include<iostream>
#include<algorithm>
#include<fstream>
#include<sstream>
#include<omp.h>
#include"mpi.h"
using namespace std;
const int col = 3799, elinenum = 1953, num_thread = 8,mpisize=8;//����������Ԫ����
bool parallel = 1;
int isupgrade;
int tmp = 0;
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
	bool isnull() {  //�жϵ�ǰ���Ƿ�Ϊ����
		if (mycol == -1)return 1;
		return 0;
	}
	void insert(int x) { //���ݶ���
		if (mycol == -1)mycol = x;
		int a = x / 32, b = x % 32;
		mybyte[a] |= (1 << b);
	}
	void doxor(bitmatrix b) {  //�����������������ڽ�����ڱ�ʵ���У�ֻ�б���Ԫ����ִ����һ����,����������Ҫ��������
		for (int i = 0; i < bytenum; i++)
			mybyte[i] ^= b.mybyte[i];
		for (int i = bytenum - 1; i >= 0; i--)
			for (int j = 31; j >= 0; j--)
				if ((mybyte[i] & (1 << j)) != 0) {
					mycol = i * 32 + j;
					return;
				}
		mycol = -1;
	}
};
bitmatrix *eliminer = new bitmatrix[col], *eline = new bitmatrix[elinenum], *pass = new bitmatrix[mpisize];
void readdata() {
	ifstream ifs;
	ifs.open("D:\\VS��Ŀ\\mpi1\\x64\\Debug\\eliminer1.txt");  //��Ԫ��
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
	ifs2.open("D:\\VS��Ŀ\\mpi1\\x64\\Debug\\eline1.txt");     //����Ԫ��,���뷽ʽ����Ԫ�Ӳ�ͬ
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
void printres() { //��ӡ���
	for (int i = 0; i < elinenum; i++) {
		if (eline[i].isnull()) { puts(""); continue; }   //���е��������
		for (int j = bytenum - 1; j >= 0; j--) {
			for (int k = 31; k >= 0; k--)
				if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //һ��������˰�Сʱ�����ǵ���λΪ1ʱ>>�����ڳ�����
					printf("%d ", j * 32 + k);
				}
		}
		puts("");
	}
}
/*void dowork() {  //������Ԫ--����Ԫ��->��Ԫ��
	int rank;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//��ȡ��ǰ���̺�
	double st = MPI_Wtime();	 //��ʱ
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
	double ed = MPI_Wtime();	 //��ʱ
	if (rank == 0) {	//ֻ��0�Ž����������ս��
		printf("cost time:%.4lf\n", ed - st);
	}
}*/

 void dowork(){  //������Ԫ--��Ԫ��->����Ԫ��
	 int rank;
	 MPI_Status status;
	 MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//��ȡ��ǰ���̺�
	 int r1 = rank * (elinenum / mpisize), r2 = (rank == mpisize - 1) ? elinenum - 1 : (rank + 1)*(elinenum / mpisize) - 1;
	 double st = MPI_Wtime();	 //��ʱ��ʼ
	 int i, j,k;
	 for (i = col - 1; i >= 0; i--) {
		 if (!eliminer[i].isnull()) {
			 for (j = r1; j <= r2; j++) {
				 if (eline[j].mycol == i)
					 eline[j].doxor(eliminer[i]);
			 }
		 }
		 else {
			 //cout << rank<<" "<<2 << endl;
			 isupgrade = -1;
			 int t = -1;
			 if (rank != 0) {
				 for (k = r1; k <= r2; k++)
					 if (eline[k].mycol == i) {
						 pass[rank] = eline[k];
						 t = k;
						 MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
						 MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);	//����������Ԫ��ϵı���Ԫ�д���0�Ž���
						 MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
						 break;
					 }
				 if (k > r2) {
					 MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
					 MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);	//����������Ԫ��ϵı���Ԫ�д���0�Ž���
					 MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
				 }
			 }
			 else {
				 for (k = mpisize - 1; k > 0; k--) {
					 MPI_Recv(&t, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					 MPI_Recv(&pass[k].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					 MPI_Recv(&pass[k].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					 if (t != -1) {
						 eliminer[i] = pass[k];
						 isupgrade = t;
					 }
				 }
				 for (k = r1; k <= r2; k++)
					 if (eline[k].mycol == i) {
						 eliminer[i] = eline[k];
						 isupgrade = k;
						 break;
					 }
				 for (k = 1; k < mpisize; k++) {
					 int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1)*(elinenum / mpisize) - 1;
					 MPI_Send(&isupgrade, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
					 if (isupgrade != -1 && t2 >= isupgrade) {
						 MPI_Send(&eliminer[i].mybyte[0], bytenum, MPI_INT, k, 1, MPI_COMM_WORLD);
						 MPI_Send(&eliminer[i].mycol, 1, MPI_INT, k, 2, MPI_COMM_WORLD);
					 }
				 }
			 }

			 //MPI_Barrier(MPI_COMM_WORLD);
			 if (rank != 0) {
				 MPI_Recv(&isupgrade, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				 if (isupgrade != -1 && r2 >= isupgrade) {
					 MPI_Recv(&eliminer[i].mybyte[0], bytenum, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
					 MPI_Recv(&eliminer[i].mycol, bytenum, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
				 }
			 }
			 if (isupgrade != -1 && r2 >= isupgrade)
				 for (j = r1; j <= r2; j++) {
					 if (eline[j].mycol == i && j != isupgrade)
						 eline[j].doxor(eliminer[i]);
				 }
		 }
	 }
	 
	 if (rank != 0) 
		 for ( k = r1; k <= r2; k++) {
			 MPI_Send(&eline[k].mybyte[0], bytenum, MPI_INT, 0, k+3+elinenum*2, MPI_COMM_WORLD);	//����������Ԫ��ϵı���Ԫ�д���0�Ž���
			 MPI_Send(&eline[k].mycol, 1, MPI_INT, 0, k+3+elinenum*3, MPI_COMM_WORLD);
		 }
	 else 
		 for ( k = 1; k < mpisize; k++) {
			 int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1)*(elinenum / mpisize) - 1;
			 for (int q = t1; q <= t2; q++) {
				 MPI_Recv(&eline[q].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				 MPI_Recv(&eline[q].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			 }
		 }
	 double ed = MPI_Wtime();	 //��ʱ����
	 if (rank == 0) {	//ֻ��0�Ž����������ս��
		 printf("cost time:%.4lf\n", ed - st);
		 //printres();
	 }
 }

int main() {
	MPI_Init(0, 0);
	readdata();
	dowork();
	MPI_Finalize();
	return 0;
}
