#include <iostream>
using namespace std;


double **interp(double **M,int N){
    double **m;
    //distribute bytes to newly created m
    m = (double **)malloc(N*sizeof(double *));
    for (int i=0;i<N;i++)
        m[i] = (double *)malloc(N*sizeof(double));
    
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            m[i][j] = *( (double *)M + N*i + j ) + 1;
    
    return m;
}

int main(){
    
    double a[3][3];
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            a[i][j] = i + j;
    double **b = interp((double **)a,3);
    
    return 0;
}
