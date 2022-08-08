#include<omp.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include<iostream>
//#include<arm_neon.h>
const int n=500,num_thread=8;
float m1[n][n],ini[n][n],time1,time2,time3;
bool parallel=1;
void init(int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            ini[i][j]=0;
        ini[i][i]=1.0;
        for(int j=i+1;j<n;j++)
            ini[i][j]=rand()%1000+1;
    }
    for(int k=0;k<n;k++)
        for(int i=k+1;i<n;i++)
            for(int j=0;j<n;j++)
                ini[i][j]=int((ini[i][j]+ini[k][j]))%1000+1.0;
}
void reset(int n){
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            m1[i][j]=ini[i][j];
}
void serial(int n){ //串行
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          for(int j=k+1;j<n;j++)
              m1[k][j]/=m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
              for(int j=k+1;j<n;j++)
                  m1[i][j]-=m1[i][k]*m1[k][j];
              m1[i][k]=0;
          }
      }
      num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time1=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time1<<" num: "<<num<<" serial: "<<time1/num<<std::endl;
    time1=time1/num;
    //  for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++)
    //         std::cout<<m1[i][j]<<" ";
    //     std::cout<<std::endl;
    //     }
}
void openmp1(int n){ //openmp实验手册版本
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int i,j,k;
    #pragma omp parallel if(parallel),num_threads(num_thread),private(i,j,k)
    for(k=0;k<n;k++){
        #pragma omp single    //串行部分
        for(j=k+1;j<n;j++)
            m1[k][j]/=m1[k][k];
        m1[k][k]=1.0;
        #pragma omp for     //并行部分
        for(i=k+1;i<n;i++){
            for(j=k+1;j<n;j++)
                m1[i][j]-=m1[i][k]*m1[k][j];
            m1[i][k]=0;
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time2=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"openmp: "<<time2<<std::endl;
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++)
    //         std::cout<<m1[i][j]<<" ";
    //     std::cout<<std::endl;
    //     }
}

void openmp2(int n){ //openmp+simd
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int i,j,k;
    #pragma omp parallel if(parallel),num_threads(num_thread),private(i,j,k)
    for(k=0;k<n;k++){
        //#pragma omp simd  1.自动simd
        //#pragma omp master    2.master
        #pragma omp single
        for(j=k+1;j<n;j++)
            m1[k][j]/=m1[k][k];
        m1[k][k]=1.0;
        //#pragma omp for simd     //并行部分 1.自动simd
        //#pragma omp for schedule(guided)  3.guided
        //#pragma omp for schedule(dynamic) 4.dynamic
        //#pragma omp for schedule(runtime) 5.runtime
        //#pragma omp for schedule(auto) 6.auto
        #pragma omp for schedule(auto)
        for(i=k+1;i<n;i++){
            for(j=k+1;j<n;j++)
                m1[i][j]-=m1[i][k]*m1[k][j];
            m1[i][k]=0;
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time3=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"openmp+simd: "<<time3<<std::endl;
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++)
    //         std::cout<<m1[i][j]<<" ";
    //     std::cout<<std::endl;
    //     }
}
int main(){
    using namespace std::chrono;
    init(n);
    reset(n);
    serial(n);
    float ttime=time1;
    reset(n);
    openmp1(n);
    std::cout<<ttime/time2<<std::endl;
    reset(n);
    openmp2(n);
    std::cout<<ttime/time3<<std::endl;
    system("pause");
    return 0;
}