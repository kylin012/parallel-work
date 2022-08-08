#include<iostream>
#include<arm_neon.h>
#include <ctime>
#include <ratio>
#include <chrono>
const int N=5001;
float m1[N][N],m2[N][N],m3[N][N],m4[i][j]; //对应普通算法、特殊算法、普通串行算法、特殊串行算法
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
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            m3[i][j]=m1[i][j];
}
void serial(int n){ //串行
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int k=0;k<n;k++){
        for(int j=k+1;j<n;j++)
            m3[k][j]/=m3[k][k];
        m3[k][k]=1.0;
        for(int i=k+1;i<n;i++){
            for(int j=k+1;j<n;j++)
                m3[i][j]-=m3[i][k]*m3[k][j];
            m3[i][k]=0;
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    std::cout<<"serial: "<<duration_cast<duration<double>>(t2-t1).count()<<std::endl;
}
void simple(int n){ //普通并行
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int k=0;k<n;k++){
        float32x4_t vt=vld1q_dup_f32(m1[k]+k);
        int j;
        for(j=k;j+4<=n;j+=4){
            float32x4_t va =vld1q_f32(m1[k]+j);
            va=vdivq_f32(va,vt);
            vst1q_f32(m1[k]+j,va);
        }
        for(;j<n;j++)
            m1[k][j]=m1[k][j]/m1[k][k];
        m1[k][k]=1.0;
        for(int i=k+1;i<n;i++){
            float32x4_t vaik=vld1q_dup_f32(m1[i]+k);
            int j;
            for(j=k+1;j+4<=n;j+=4){
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
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    std::cout<<"simple: "<<duration_cast<duration<double>>(t2-t1).count()<<std::endl;
}
void special(int n){

}
int main(){
    using namespace std::chrono;
    std::cout<<"input n"<<std::endl;
    int n;
    std::cin>>n;
    reset(n);
    serial(n);
    simple(n);
    
}