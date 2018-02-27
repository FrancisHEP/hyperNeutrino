#include <iostream>
#include <cmath>
using namespace std;

double funcB(double x,double y,double z)
{
    return 100*x+10*y+z;
}

double funcA(double x,double y)
{
    return funcB(x,y,z);
}

int main()
{
    
    return 0;
}
