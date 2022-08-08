#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
double a[10001];
int main()
{
    using namespace std::chrono;
    for(int i=1;i<=10000;i++)
        a[i]=i;
    std::cout<<"common:"<<std::endl;
    int n=4;
    while(n<=8192){
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
    system("pause");
}
