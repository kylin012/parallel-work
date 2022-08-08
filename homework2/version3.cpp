#include<iostream>
#include<arm_neon.h>
#include <ctime>
#include <ratio>
#include <chrono>
 const int n=100;
float m1[n][n],time1,time3; //time1、3对应并行算法、串行算法
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
void parallel(int n){ //普通并行
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          float32x4_t vt=vld1q_dup_f32(m1[k]+k);
          int j;
          if(k%4!=0){
              for(j=k;j<(k+4)/4*4;j++){
                  m1[k][j]/=m1[k][k];
              }
          }
          for(;j+4<=n;j+=4){
              float32x4_t va =vld1q_f32(m1[k]+j);
              va=vdivq_f32(va,vt);
              vst1q_f32(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
              float32x4_t vaik=vld1q_dup_f32(m1[i]+k);
              if((k+1)%4!=0){
                  for(j=k+1;j<(k+5)/4*4;j++){
                       m1[i][j]-=m1[i][k]*m1[k][j];
                  }
              }
              for(;j+4<=n;j+=4){
                  float32x4_t vakj=vld1q_f32(m1[k]+j);
                  float32x4_t vaij=vld1q_f32(m1[i]+j);
                  float32x4_t vx= vmulq_f32(vaik,vakj);
                  vaij = vsubq_f32(vaij,vx);
                  vst1q_f32(m1[i]+j,vaij); 
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
    std::cout<<"total: "<<time3<<" num: "<<num<<" parallel: "<<time3/num<<std::endl;
    time3=time3/num;
}
int main(){
    using namespace std::chrono;
    reset(n);
    serial(n);
    reset(n);
    parallel(n);
    std::cout<<time1/time3<<std::endl;
}