#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double** fun(){
    double **a;
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++) a[i][j] = i+j+10;
    }
    cout<<a[2][2]<<endl;

    return a;
}

int main(){
    double **b = fun();
    ofstream ofile;
    ofile.open("example2",ios::out|ios::app);    //it means if not created, create now
    //cout<<b[0][0]<<endl;
    /*
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            ofile<<b[i][j]<<endl;
        }
    }*/
    return 0;
}

