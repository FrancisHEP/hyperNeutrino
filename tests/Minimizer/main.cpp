//    The main function for the simulation of 2 dimensional neutrino
//    oscillation.
#include <iostream>
#include <vector>
#include <algorithm>
#include "TF2.h"
#include "TMinuit.h"
#include "TRandom.h"
using namespace std;
static vector<double> I(5,1);   //Initialization vector, eg. P=multi(I,1);
static double b = 100;
extern vector<double> nOrigin(5,0);     //MUST INITIALIZE IT!!!(in case of FATAL error)
void PrintVector(vector<double> a);     //show vector
double SumVector(vector<double> a);     //find sum for the vector
vector<double> multi(vector<double> a, double x);       //multi in Vectors
int GetPoisson(double mean);                            //generate poisson randoms
vector<double> GetPoisson(vector<double> mean);         //generate poisson randoms V
vector<double> FastIntegral(double x1, double x2);


void PrintVector(vector<double> a)
{   // mathematic tools
    // show vector
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

double SumVector(vector<double> a)
{   // mathematic tools
    double sum = 0;
    for (int i=0;i<a.size();i++) sum = sum + a[i];
    return sum;
}

vector<double> multi(vector<double> a, double x)
{   // mathematic tools
    // b = ax for 1D vectors
    vector<double> b(a);
    for (int i=0;i<a.size();i++) b[i] = a[i]*x;
    return b;
}

int GetPoisson(double mean)
{   // mathematic tools
    TRandom* r = new TRandom();
    r -> SetSeed();
    int expect = r -> Poisson(mean);
    return expect;
}

vector<double> GetPoisson(vector<double> mean)
{   // mathematic tools
    vector<double> expect(mean.size(),0);
    for (int i=0;i<mean.size();i++)
    {
        TRandom* r = new TRandom();
        r -> SetSeed();
        expect[i] = r -> Poisson(mean[i]);
    }
    return expect;
}

vector<double> GetPoisson(vector<double> mean, double b)
{   // poisson with background
    vector<double> expect(mean.size(),0);
    for (int i=0;i<mean.size();i++)
    {
        TRandom* r = new TRandom();
        r -> SetSeed();
        expect[i] = r -> Poisson(mean[i]+b);
    }
    return expect;
}

vector<double> FastIntegral(double x1, double x2)
{   //  integral tool
    vector<double> result(5,0);
    TF2 * func1 = new TF2("func1","[0]*pow(sin(1.27*x/y*[1]),2)");
    func1->SetParameters(x1,x2);
    double L1 = 0.6;
    double L2 = 1.0;
    for (int i=0;i<5;i++)
    {
        double E1 = 10 + 10*i;
        double E2 = 20 + 10*i;
        result[i] = func1->Integral(L1,L2,E1,E2,1e-5);
    }
    return result;
}

vector<double> L(vector<double> n, vector<double> mu, double b)
{   //  tools for dchi2, CAUTION: L is Chi, Chi is L.
    vector<double> L(n.size(),0);
    for (int i=0;i<n.size();i++) L[i] = pow(n[i]-mu[i]-b,2)/(mu[i]+b);
    return L;
}

class PhasePoint
{
    //    (x1,x2) = (sin2, dm2), is the position of phase point
    //    (sin(2*theta)^2,dm^2)
    //    for example,
    //    sin2 = pow(10,-1.33), dm2  = pow(10,2);
public:
    void SetValue (double x1, double x2);
    void PrintValue();
    void FindDchi2();
    
    double sin2;
    double dm2;
    double dchi2;
    
    vector<double> P;       //Probability in 5 bands, [1 2 3 4 5]
    vector<double> muTrue;  //P[5]*10000
};

void PhasePoint::SetValue(double x1, double x2)
{   //SET value for x1,x2 and INITIALIZE P,muTrue
    sin2 = x1;
    dm2  = x2;
    P = FastIntegral(sin2,dm2);
    muTrue = multi(P,10000);
}

void PhasePoint::PrintValue()
{
    printf("sin(2*theta)^2 = %f, dm^2 = %f\n",sin2,dm2);
    printf("muTrue = ");PrintVector(muTrue);
}

//  tools for dchi2
double funcL (double x1, double x2)
{
    PhasePoint M;
    M.SetValue(x1,x2);
    vector<double> mu = M.muTrue;
    return SumVector(L(nOrigin,mu,b));
}

void newFCN(int& nDim, double* gout, double& result, double p[], int flg)
{   result = funcL(p[0],p[1]);  }

vector<double> FindMinimum(double x1, double x2)
{
    //int flag = 0;
    TMinuit John(2);
    John.mninit(5,6,7); //initialization?
    John.SetPrintLevel(-1);
    John.SetFCN(newFCN);
    John.DefineParameter(0,"x1",x1,1e-6,0,1);
    John.DefineParameter(1,"x2",x2,1e-3,0,1e3);
    John.SetErrorDef(1);
    John.SetMaxIterations(100000);
    int err1 = John.Migrad();   //NEVER DELETE IT!!!!!!
    while (0) err1 = 0;
    
    double x1_,x1_err;
    double x2_,x2_err;
    John.GetParameter(0,x1_,x1_err);
    John.GetParameter(1,x2_,x2_err);
    
    vector<double> vec(2,0);
    vec[0] = x1_;
    vec[1] = x2_;
    return vec;
}

void PhasePoint::FindDchi2()
{
    int N = 1000;  //measures
    int lim = int(0.9*N);
    vector<double> chi(N,0);
    for (int i=0;i<N;i++)
    {
        nOrigin = GetPoisson(muTrue,b);
        vector<double> Pnew = FindMinimum(sin2,dm2);
        vector<double> muBest = multi(FastIntegral(Pnew[0],Pnew[1]),10000);
        double L1 = SumVector(L(nOrigin,muTrue,b));
        double L2 = SumVector(L(nOrigin,muBest,b));
        chi[i] = L1 - L2;
    }
    //PrintVector(chi);
    sort(chi.begin(),chi.end());
    dchi2 = chi[lim];
}

int main()
{
    /*
     vector<double> dchi2V;
    for (int i=0;i<100;i++){
        for (int j=0;j<100;j++){
            double x1 = 0.01*i;
            double x2 = 10*j;
             PhasePoint M;
             M.SetValue(x1,x2);
             M.FindDchi2();
             dchi2V.push_back(M.dchi2);
        }
    }*/
    double x1 = pow(10,-3.5);
    double x2 = pow(10,0.5);
    printf("M:\n");
    PhasePoint M;
    M.SetValue(x1,x2);
    M.PrintValue();
    M.FindDchi2();
    printf("dchi2 = %f\n",M.dchi2);

    return 0;
}
