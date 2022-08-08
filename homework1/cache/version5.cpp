#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
double res[2][40];
int main()
{
    using namespace std::chrono;
    std::cout<<std::endl<<"common:"<<std::endl;
    int n=10,k=0;
    while(n<=10000){
        int *sum = new int[n],*a = new int[n],**b = new int*[n];
        for(int i=0;i<n;i++){
            sum[i]=0,a[i]=i%100;
            b[i]=new int[n];
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
        res[0][k++]=duration_cast<duration<double>>(t2 - t1).count()/counter;
        std::cout <<"n= "<<n<<" counter= "<<counter<<" time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        if(n<100)n+=10;
        else if(n<1000)n+=100;
        else n+=1000;
    }
    std::cout<<std::endl<<"cache:"<<std::endl;
    n=10,k=0;
    while(n<=10000){
        int *sum = new int[n],*a = new int[n],**b = new int*[n];
        for(int i=0;i<n;i++){
            sum[i]=0,a[i]=i%100;
            b[i]=new int[n];
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
        res[1][k++]=duration_cast<duration<double>>(t2 - t1).count()/counter;
        std::cout <<"n= "<<n<<" counter= "<<counter<<" total time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        if(n<100)n+=10;
        else if(n<1000)n+=100;
        else n+=1000;
    }
    for(int i=0;i<k;i++)
        std::cout<<res[0][i]/res[1][i]<<std::endl;
    system("pause");   
}