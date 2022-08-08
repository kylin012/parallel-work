#include<omp.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include<iostream>
#include <emmintrin.h>
#include <immintrin.h>
const int n=10,num_thread=8;
float m1[n][n],ini[n][n],time1,time2,time3,time4,time5;
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
     for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<m1[i][j]<<" ";
        std::cout<<std::endl;
        }
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
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<m1[i][j]<<" ";
        std::cout<<std::endl;
        }
}

void openmp2(int n){ //openmp+sse
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int i,j,k;
    #pragma omp parallel if(parallel),num_threads(num_thread),private(i,j,k)
      for(int k=0;k<n;k++){
          #pragma omp single
          for(int i=0;i<2;i++){
            if(i==0)  
            for(j=k+1;j+4<=n;j+=4){
                __m128 vt= _mm_set1_ps(m1[k][k]);
                __m128 va =_mm_loadu_ps(m1[k]+j);
                va=_mm_div_ps(va,vt);
                _mm_storeu_ps(m1[k]+j,va);
            }
            else
            for(;j<n;j++)
                m1[k][j]=m1[k][j]/m1[k][k];
          }
          m1[k][k]=1.0;
          #pragma omp for
          for(i=k+1;i<n;i++){
              __m128 vaik=_mm_set1_ps(m1[i][k]);
              int j;
              for(j=k+1;j+4<=n;j+=4){
                  __m128 vakj=_mm_loadu_ps(m1[k]+j);
                  __m128 vaij=_mm_loadu_ps(m1[i]+j);
                  __m128 vx= _mm_mul_ps(vaik,vakj);
                  vaij = _mm_sub_ps(vaij,vx);
                  _mm_storeu_ps(m1[i]+j,vaij); 
              }
              for(;j<n;j++)
                  m1[i][j]=m1[i][j]-m1[i][k]*m1[k][j];
              m1[i][k]=0;
          }
        }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time3=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"openmp+sse: "<<time3<<std::endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<m1[i][j]<<" ";
        std::cout<<std::endl;
        }
}

void openmp3(int n){    //openmp+avx
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int i,j,k;
    #pragma omp parallel if(parallel),num_threads(num_thread),private(i,j,k)
      for(int k=0;k<n;k++){
          #pragma omp single
          for(int i=0;i<2;i++){
            if(i==0)  
            for(j=k+1;j+8<=n;j+=8){
                __m256 vt = _mm256_set1_ps(m1[k][k]);
                __m256 va =_mm256_loadu_ps(m1[k]+j);
                va=_mm256_div_ps(va,vt);
                _mm256_storeu_ps(m1[k]+j,va);
            }
            else
            for(;j<n;j++)
                m1[k][j]=m1[k][j]/m1[k][k];
            }
          m1[k][k]=1.0;
          #pragma omp for
          for(i=k+1;i<n;i++){
              __m256 vaik=_mm256_set1_ps(m1[i][k]);
              int j;
              for(j=k+1;j+8<=n;j+=8){
                  __m256 vakj=_mm256_loadu_ps(m1[k]+j);
                  __m256 vaij=_mm256_loadu_ps(m1[i]+j);
                  __m256 vx= _mm256_mul_ps(vaik,vakj);
                  vaij = _mm256_sub_ps(vaij,vx);
                  _mm256_storeu_ps(m1[i]+j,vaij); 
              }
              for(;j<n;j++)
                  m1[i][j]=m1[i][j]-m1[i][k]*m1[k][j];
              m1[i][k]=0;
          }
        }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time4=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"openmp+avx: "<<time4<<std::endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<m1[i][j]<<" ";
        std::cout<<std::endl;
        }
}

void openmp4(int n){ //openmp+avx512
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int i,j,k;
    #pragma omp parallel if(parallel),num_threads(num_thread),private(i,j,k)
      for(int k=0;k<n;k++){
          #pragma omp single
          for(int i=0;i<2;i++){
            if(i==0) 
            for(j=k+1;j+16<=n;j+=16){
                __m512 vt= _mm512_set1_ps(m1[k][k]);
                __m512 va =_mm512_loadu_ps(m1[k]+j);
                va=_mm512_div_ps(va,vt);
                _mm512_storeu_ps(m1[k]+j,va);
            }
            else
            for(;j<n;j++)
                m1[k][j]=m1[k][j]/m1[k][k];
          }
          m1[k][k]=1.0;
          #pragma omp for
          for(i=k+1;i<n;i++){
              __m512 vaik=_mm512_set1_ps(m1[i][k]);
              int j;
              for(j=k+1;j+16<=n;j+=16){
                  __m512 vakj=_mm512_loadu_ps(m1[k]+j);
                  __m512 vaij=_mm512_loadu_ps(m1[i]+j);
                  __m512 vx= _mm512_mul_ps(vaik,vakj);
                  vaij = _mm512_sub_ps(vaij,vx);
                  _mm512_storeu_ps(m1[i]+j,vaij); 
              }
              for(;j<n;j++)
                  m1[i][j]=m1[i][j]-m1[i][k]*m1[k][j];
              m1[i][k]=0;
          }
        }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time5=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"openmp+avx512: "<<time5<<std::endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<m1[i][j]<<" ";
        std::cout<<std::endl;
        }
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
    openmp3(n);
    std::cout<<ttime/time4<<std::endl;
    //openmp4(n);
    //std::cout<<ttime/time5<<std::endl;
    system("pause");
    return 0;
}
