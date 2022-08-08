#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
int main()
{
    using namespace std::chrono;
    std::cout<<"common:"<<std::endl;
    int n=4;
    while(n<=131072){
        double *a = new double[n];
        for(int i=0;i<n;i++)
            a[i]=i;
        int counter=0;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
            counter++;
            double sum=0;
            for (int i = 0; i < n; i++)
                sum += a[i];
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        std::cout <<"n= "<<n<<" counter= "<<counter<<" time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        n*=2;
    }
    std::cout<<std::endl<<"scalar:"<<std::endl;
    n=4;
    while(n<=131072){
        double *a = new double[n];
        for(int i=0;i<n;i++)
            a[i]=i;
        int counter=0;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<1){
            counter++;
            double sum=0;
            double sum1 = 0, sum2 = 0;
        for (int i = 0;i < n; i += 2) {
            sum1 += a[i];
            sum2 += a[i + 1];
        }
            sum = sum1 + sum2;
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        std::cout <<"n= "<<n<<" counter= "<<counter<<" time: "<< duration_cast<duration<double>>(t2 - t1).count()<<" single time:"<<
        duration_cast<duration<double>>(t2 - t1).count()/counter<<std::endl;
        n*=2;
    }
    
    system("pause");
}
