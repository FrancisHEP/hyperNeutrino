#include <iostream>
#include "TH2.h"
using namespace std;

int main(){
    
    //generate a matrix
    int N = 3;
    double a[N][N];
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            a[i][j] = i+j;
    
    //do
    double eps = 1e-8;
    
    TH2D* fun = new TH2D("","",N,0,1,N,0,1);
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            fun->SetBinContent(i+1,j+1,a[i][j]);
    double b = fun->Interpolate(x,y);
    
    return 0;
}



