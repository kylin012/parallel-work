#include<iostream>
#include <emmintrin.h>
#include <immintrin.h>
#include <ctime>
#include <ratio>
#include <chrono>
 const int n=128;
float m1[n][n],time1,time2,time3,time4; //time1、2、3、4对应串行算法、sse算法、avx算法、avx512算法
void reset(int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            m1[i][j]=0;
        m1[i][i]=1.0;
        for(int j=i+1;j<n;j++)
            m1[i][j]=rand()%1000+1;
    }
    for(int k=0;k<n;k++)
        for(int i=k+1;i<n;i++)
            for(int j=0;j<n;j++)
                m1[i][j]+=m1[k][j];
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
}
void sse(int n){ //普通并行
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          float temp[4]={m1[k][k],m1[k][k],m1[k][k],m1[k][k]};
          __m128 vt= _mm_loadu_ps(temp);
          int j;
          for(j=k;j+4<=n;j+=4){
              __m128 va =_mm_loadu_ps(m1[k]+j);
              va=_mm_div_ps(va,vt);
              _mm_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
              float temp2[4]={m1[i][k],m1[i][k],m1[i][k],m1[i][k]};
              __m128 vaik=_mm_loadu_ps(temp2);
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
        num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time2=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time2<<" num: "<<num<<" sse: "<<time2/num<<std::endl;
    time2=time2/num;
}
void avx(int n){
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          float temp[8]={m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k]};
          __m256 vt= _mm256_loadu_ps(temp);
          int j;
          for(j=k;j+8<=n;j+=8){
              __m256 va =_mm256_loadu_ps(m1[k]+j);
              va=_mm256_div_ps(va,vt);
              _mm256_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
              float temp2[8]={m1[i][k],m1[i][k],m1[i][k],m1[i][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k]};
              __m256 vaik=_mm256_loadu_ps(temp2);
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
        num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time3=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time3<<" num: "<<num<<" avx: "<<time3/num<<std::endl;
    time3=time3/num;
}
// void avx512(int n){
//     using namespace std::chrono;
//     int num=0;
//     high_resolution_clock::time_point t1 = high_resolution_clock::now();
//     while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
//       for(int k=0;k<n;k++){
//           float temp[16]={m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k]};
//           __m512 vt= _mm512_loadu_ps(temp);
//           int j;
//           for(j=k;j+16<=n;j+=16){
//               __m512 va =_mm512_loadu_ps(m1[k]+j);
//               va=_mm512_div_ps(va,vt);
//               _mm512_storeu_ps(m1[k]+j,va);
//           }
//           for(;j<n;j++)
//               m1[k][j]=m1[k][j]/m1[k][k];
//           m1[k][k]=1.0;
//           for(int i=k+1;i<n;i++){
//               float temp2[16]={m1[i][k],m1[i][k],m1[i][k],m1[i][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k],m1[i][k],m1[i][k],m1[i][k],m1[i][k],m1[k][k],m1[k][k],m1[k][k],m1[k][k]};
//               __m512 vaik=_mm512_loadu_ps(temp2);
//               int j;
//               for(j=k+1;j+16<=n;j+=16){
//                   __m512 vakj=_mm512_loadu_ps(m1[k]+j);
//                   __m512 vaij=_mm512_loadu_ps(m1[i]+j);
//                   __m512 vx= _mm512_mul_ps(vaik,vakj);
//                   vaij = _mm512_sub_ps(vaij,vx);
//                   _mm512_storeu_ps(m1[i]+j,vaij); 
//               }
//               for(;j<n;j++)
//                   m1[i][j]=m1[i][j]-m1[i][k]*m1[k][j];
//               m1[i][k]=0;
//           }
//         }
//         num++;
//     }
//     high_resolution_clock::time_point t2 = high_resolution_clock::now();
//     time4=duration_cast<duration<double>>(t2-t1).count();
//     std::cout<<"total: "<<time4<<" num: "<<num<<" avx512: "<<time4/num<<std::endl;
//     time4=time4/num;
// }
int main(){
    using namespace std::chrono;
    reset(n);
    serial(n);
    reset(n);
    sse(n);
    std::cout<<"serial/sse: "<<time1/time2<<std::endl;
    double ttime = time1;
    reset(n);
    avx(n);
    std::cout<<"serial/avx: "<<ttime/time3<<std::endl;
    system("pause");
    //reset(n);
    // std::cout<<"serial/avx512: "<<time1/time4<<std::endl;
    // avx512(n);
}