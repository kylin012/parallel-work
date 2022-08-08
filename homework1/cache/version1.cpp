#include <iostream>
#include<stdlib.h>
#include<time.h>
using namespace std ;
const int N = 10240; // matrix size
double b[1001][1001],a[1001],sum[1001];
int main()
{
    for(int i=1;i<=1000;i++){
        a[i]=i;
        for(int j=1;j<=1000;j++)
            b[i][j]=i+j;
    }
    int n=10,step = 10;
    clock_t start , finish ;
    while(n<=1000){
        int counter=0;
        start=clock();
        while(clock()-start<10){
            counter++;
            for(int i = 0; i < n; i++){
                sum[i] = 0.0;
                for(int j = 0; j < n; j++)
                sum[i] += b[j][i] * a[j];
            }
        }
        finish=clock();
        cout <<"n= "<<n<<" counter= "<<counter<<" time: "<<  ( finish - start)/float (CLOCKS_PER_SEC)<<endl;
        if(n<100)n+=10;
        else n+=100;
    }
}