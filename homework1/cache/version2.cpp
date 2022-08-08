#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
const int N = 10240; // matrix size
double b[1001][1001],a[1001],sum[1001];
int main()
{
    using namespace std::chrono;
    for(int i=1;i<=1000;i++){
        a[i]=i;
        for(int j=1;j<=1000;j++)
            b[i][j]=i+j;
    }
    int n=10,step = 10;
    while(n<=1000){
        int counter=0;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        while(duration_cast<duration<double>>(high_resolution_clock::now() - t1).count()<0.1){
            counter++;
            for(int i = 0; i < n; i++){
                sum[i] = 0.0;
                for(int j = 0; j < n; j++)
                sum[i] += b[j][i] * a[j];
            }
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        std::cout <<"n= "<<n<<" counter= "<<counter<<" time: "<< duration_cast<duration<double>>(t2 - t1).count()<<std::endl;
        if(n<100)n+=10;
        else n+=100;
    }
}