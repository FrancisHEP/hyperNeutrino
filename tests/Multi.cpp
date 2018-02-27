#include <iostream>
#include <vector>
using namespace std;

// b = ax for 1D vectors

vector<double> multi(vector<double> a, double x)
{
    vector<double> b(a);
    for (int i=0;i<a.size();i++) b[i] = a[i]*x;
    return b;
}

int main ()
{
    double ini[3] = {1,2,3};
    vector<double> a(ini,ini+3);
    vector<double> b;
    b = multi(a,2);
    for (int i=0;i<b.size();i++) printf("%f\n",b[i]);
    
    return 0;
}
