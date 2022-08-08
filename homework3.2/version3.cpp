#include<iostream>
#include <emmintrin.h>
#include <immintrin.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include<pthread.h>
#include<semaphore.h>
 const int n=100,num_thread=8;
float m1[n][n],ini[n][n],time1,time2,time3,time4,time5; //time1、3对应并行算法、串行算法
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
}

void sse(int n){    //sse
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          __m128 vt= _mm_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+4<=n;j+=4){
              __m128 va =_mm_loadu_ps(m1[k]+j);
              va=_mm_div_ps(va,vt);
              _mm_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
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
        num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time2=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time2<<" num: "<<num<<" sse: "<<time2/num<<std::endl;
    time2=time2/num;
}

void avx(int n){    //avx
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          __m256 vt = _mm256_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+8<=n;j+=8){
              __m256 va =_mm256_loadu_ps(m1[k]+j);
              va=_mm256_div_ps(va,vt);
              _mm256_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
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
        num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time3=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time3<<" num: "<<num<<" avx: "<<time3/num<<std::endl;
    time3=time3/num;
}

void avx512(int n){ //avx512
    using namespace std::chrono;
    int num=0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
      for(int k=0;k<n;k++){
          __m512 vt= _mm512_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+16<=n;j+=16){
              __m512 va =_mm512_loadu_ps(m1[k]+j);
              va=_mm512_div_ps(va,vt);
              _mm512_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;
          for(int i=k+1;i<n;i++){
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
        num++;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time4=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"total: "<<time4<<" num: "<<num<<" avx512: "<<time4/num<<std::endl;
    time4=time4/num;
}
typedef struct{
    int t_id;   //线程id
}threadparam_t;
sem_t sem_div;  //信号量定义
sem_t sem_elimination;
pthread_barrier_t barrier;
void *threadfunc(void *param){      //pthread
    threadparam_t *p=(threadparam_t *)param;
    int t_id=p->t_id;
    for(int k=0;k<n;k++){   //除法操作
        if(t_id==0){
            for(int j=k+1;j<n;j++)
                m1[k][j]/=m1[k][k];
            m1[k][k]=1.0;  
            }
        else sem_wait(&sem_div);    //阻塞，等待完成除法操作
        if(t_id==0){
            for(int i=1;i<num_thread;i++)
                sem_post(&sem_div);     //唤醒，准备进行消去操作
        }
        for(int i=k+1+t_id;i<n;i+=num_thread){  //处理消去操作时，每个线程以8行为间隔处理某一行
             for(int j=k+1;j<n;j++)
                m1[i][j]-=m1[i][k]*m1[k][j];
            m1[i][k]=0;
        }
        pthread_barrier_wait(&barrier);
    }   
    pthread_exit(NULL);
}

void *threadfunc2(void *param){      //pthread+sse
    threadparam_t *p=(threadparam_t *)param;
    int t_id=p->t_id;
    for(int k=0;k<n;k++){   //除法操作
        if(t_id==0){
            __m128 vt= _mm_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+4<=n;j+=4){
              __m128 va =_mm_loadu_ps(m1[k]+j);
              va=_mm_div_ps(va,vt);
              _mm_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0; 
            }
        else sem_wait(&sem_div);    //阻塞，等待完成除法操作
        if(t_id==0){
            for(int i=1;i<num_thread;i++)
                sem_post(&sem_div);     //唤醒，准备进行消去操作
        }
        for(int i=k+1+t_id;i<n;i+=num_thread){  //处理消去操作时，每个线程以8行为间隔处理某一行
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
        pthread_barrier_wait(&barrier);
    }   
    pthread_exit(NULL);
}

void *threadfunc3(void *param){      //pthread+avx
    threadparam_t *p=(threadparam_t *)param;
    int t_id=p->t_id;
    for(int k=0;k<n;k++){   //除法操作
        if(t_id==0){
            __m256 vt = _mm256_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+8<=n;j+=8){
              __m256 va =_mm256_loadu_ps(m1[k]+j);
              va=_mm256_div_ps(va,vt);
              _mm256_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0;  
            }
        else sem_wait(&sem_div);    //阻塞，等待完成除法操作
        if(t_id==0){
            for(int i=1;i<num_thread;i++)
                sem_post(&sem_div);     //唤醒，准备进行消去操作
        }
        for(int i=k+1+t_id;i<n;i+=num_thread){  //处理消去操作时，每个线程以8行为间隔处理某一行
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
        pthread_barrier_wait(&barrier);
    }   
    pthread_exit(NULL);
}

void *threadfunc4(void *param){      //pthread+avx512
    threadparam_t *p=(threadparam_t *)param;
    int t_id=p->t_id;
    for(int k=0;k<n;k++){   //除法操作
        if(t_id==0){
            __m512 vt= _mm512_set1_ps(m1[k][k]);
          int j;
          for(j=k+1;j+16<=n;j+=16){
              __m512 va =_mm512_loadu_ps(m1[k]+j);
              va=_mm512_div_ps(va,vt);
              _mm512_storeu_ps(m1[k]+j,va);
          }
          for(;j<n;j++)
              m1[k][j]=m1[k][j]/m1[k][k];
          m1[k][k]=1.0; 
            }
        else sem_wait(&sem_div);    //阻塞，等待完成除法操作
        if(t_id==0){
            for(int i=1;i<num_thread;i++)
                sem_post(&sem_div);     //唤醒，准备进行消去操作
        }
        for(int i=k+1+t_id;i<n;i+=num_thread){  //处理消去操作时，每个线程以8行为间隔处理某一行
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
        pthread_barrier_wait(&barrier);
    }   
    pthread_exit(NULL);
}
void pthread(int n){
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    sem_init(&sem_div,0,0);     //初始化信号量
    sem_init(&sem_elimination,0,0);
    pthread_barrier_init(&barrier,NULL,num_thread);
    //创建线程
    pthread_t handles[num_thread];
    threadparam_t param[num_thread];
    for(int t_id=0;t_id<num_thread;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadfunc,(void*)&param[t_id]);
    }
    for(int i=0;i<num_thread;i++)
        pthread_join(handles[i],NULL);
    sem_destroy(&sem_div);
    sem_destroy(&sem_elimination);
    pthread_barrier_destroy(&barrier);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    time5=duration_cast<duration<double>>(t2-t1).count();
    std::cout<<"pthread: "<<time5<<std::endl;
}
int main(){
    using namespace std::chrono;
    init(n);
    reset(n);
    serial(n);
    double ttime=time1;
    reset(n);
    sse(n);
    std::cout<<ttime/time2<<std::endl;
    reset(n);
    avx(n);
    std::cout<<ttime/time3<<std::endl;
    reset(n);
    avx512(n);
    std::cout<<ttime/time4<<std::endl;
    reset(n);
    pthread(n);
    std::cout<<ttime/time5<<std::endl;
    return 0;
}