#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
int main()
{
    using namespace std::chrono;
    int n=10;
    while(n<=1000){
        double *sum = new double[n],*a = new double[n],**b = new double*[n];
        for(int i=0;i<n;i++){
            sum[i]=0,a[i]=i;
            b[i]=new double[n];
            for(int j=0;j<n;j++)
                b[i][j]=i+j;
        }
        int counter=0;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
            for(int i = 0; i < n; i++)
                sum[i] = 0;
            counter++;
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++)
                sum[i] += b[j][i] * a[j];
            }
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        std::cout <<"n= "<<n<<" counter= "<<counter<<" time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        if(n<100)n+=10;
        else n+=100;
    }
    std::cout<<std::endl<<"cache:"<<std::endl;
    n=10;
    while(n<=1000){
        double *sum = new double[n],*a = new double[n],**b = new double*[n];
        for(int i=0;i<n;i++){
            sum[i]=0,a[i]=i;
            b[i]=new double[n];
            for(int j=0;j<n;j++)
                b[i][j]=i+j;
        }
        int counter=0;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
            for(int i = 0; i <n; i++)
                sum[i] = 0;
            counter++;
            for(int j = 0; j < n; j++)
                for(int i = 0; i < n; i++)
                    sum[i] += b[j][i] * a[j];
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        std::cout <<"n= "<<n<<" counter= "<<counter<<" total time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        if(n<100)n+=10;
        else n+=100;
    }
    system("pause");   
}