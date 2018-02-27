#include <iostream>
#include <vector>
#include <TRandom.h>
using namespace std;
void PrintVector(vector<double> a);                 //show vector
vector<double> multi(vector<double> a, double x);   //multi in Vectors
static vector<double> I(5,1);   //Initialization vector, eg. P=multi(I,1);
int GetPoisson(double mean);    //generate poisson randoms

vector<double> multi(vector<double> a, double x)
{   // b = ax for 1D vectors
    vector<double> b(a);
    for (int i=0;i<a.size();i++) b[i] = a[i]*x;
    return b;
}

void PrintVector(vector<double> a)
{   //show vector
    if (a.empty()) {
        cout<<"NULL"<<endl;
    } else {
        if (a.size()>100){
            for (int i=0;i<100;i++) cout<<a[i]<<" ";
            cout<<endl<<a.size()<<" elements in the vector. "<<endl;
        } else {
            for (int i=0;i<a.size();i++) cout<<a[i]<<" ";
            cout<<endl; }
    }
}

int GetPoisson(double mean)
{
    TRandom* r = new TRandom();
    r -> SetSeed();
    int expect = r -> Poisson(mean);
    return expect;
}

vector<double> GetPoisson(vector<double> mean)
{
    vector<double> expect(mean.size(),0);
    for (int i=0;i<mean.size();i++)
    {
        TRandom* r = new TRandom();
        r -> SetSeed();
        expect[i] = r -> Poisson(mean[i]);
    }
    return expect;
}

/*
 sin(2*theta)^2 = 0.006000, dm^2 = 40.000000
 muTrue:
 95.8257 218.041 196.977 147.627 109.599
 n:
 89 220 202 172 109
*/
double a[5] = {95.8257,218.041,196.977,147.627,109.599};
vector<double> muTrue(a,a+5);

int main()
{
    vector<double> n(5,0);
    for (int i=0;i<5;i++) n[i] = GetPoisson(muTrue[i]);
    PrintVector(n);
    
    vector<double> N(5,0);
    N = GetPoisson(muTrue);
    PrintVector(N);
    
    return 0;
}
