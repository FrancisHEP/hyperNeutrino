#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double fun(double x,double y){
    return pow(x,2)+y;
}

int main(){
    
    ofstream ofile;
    ofile.open("example",ios::out|ios::app);    //it means if not created, create now
    
    for (int i=0;i<10;i++){
        for (int j=0;j<3;j++){
            ofile<<fun(i,j)<<endl;
        }
    }
    return 0;
}

