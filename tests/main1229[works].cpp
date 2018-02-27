//    The main function for the simulation of 2 dimensional neutrino
//    oscillation.
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "TF2.h"
#include "TH2.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
static vector<double> I(5,1);   //Initialization vector, eg. P=multi(I,1);
static double b = 100;
extern vector<double> nOrigin(5,0);     //MUST INITIALIZE!!!in case of FATAL error
void PrintVector(vector<double> a);     //show vector
double SumVector(vector<double> a);     //find sum for the vector
vector<double> multi(vector<double> a, double x);       //multi in Vectors
int GetPoisson(double mean);                            //generate poisson randoms
vector<double> GetPoisson(vector<double> mean);
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

double SimpFunc(double x2,double E)
{
    double P;
    P = 0.2-E/(4*1.27*x2)*(sin(2*1.27*x2/E)-sin(1.2*1.27*x2/E));
    return P;
}

vector<double> SimptonIntegral(double x1,double x2)
{   //  integral hinted by dana
    vector<double> result(5,0);
    int N = 10;
    double h = 10/double(N);
    double E = 0;
    for (int i=0;i<5;i++){
        double sum = 0;
        double s = 0;
        for (int j=1;j<=N+1;j++){
            E = 10*(i+1) + 10*(j-1)/double(N);
            s = SimpFunc(x2,E)*h/3;
            if (j==1||j==N+1)   s = s*1;
            else {
                if (j%2==0)     s = s*4;
                else            s = s*2;
            }
            sum += s;
        }
        result[i] = sum * x1;
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
    void FindDchi2(int samples);
    
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
    P = SimptonIntegral(sin2,dm2);
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

void PhasePoint::FindDchi2(int samples)
{
    int N = samples;  //measures
    int lim = int(0.9*N);
    vector<double> chi(N,0);
    for (int i=0;i<N;i++)
    {
        nOrigin = GetPoisson(muTrue,b);
        vector<double> Pnew = FindMinimum(sin2,dm2);
        vector<double> muBest = multi(SimptonIntegral(Pnew[0],Pnew[1]),10000);
        double L1 = SumVector(L(nOrigin,muTrue,b));
        double L2 = SumVector(L(nOrigin,muBest,b));
        chi[i] = L1 - L2;
    }
    //PrintVector(chi);
    sort(chi.begin(),chi.end());
    dchi2 = chi[lim];
}

void PlotDchi2(int N,int samples)   //points,samples
{   int count = 0; printf("%d*%d points in TOTAL.\n",N,N);
    vector<double> dchi2V;
    double Matrix[N][N];
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){  count += 1; printf("count = %d\n",count);
            double x1 = 4*(double(i)+0.5)/double(N)-4;
            double x2 = 3*(double(j)+0.5)/double(N);
            PhasePoint M;
            M.SetValue(pow(10,x1),pow(10,x2));
            //printf("M:\n"); M.PrintValue();
            M.FindDchi2(samples);                   ////////samples
            dchi2V.push_back(M.dchi2);
            Matrix[i][j] = M.dchi2;
        }
    }
    //PrintVector(dchi2V);
    
    TCanvas cvs;
    cvs.SetWindowSize(200,200);
    TH2D* hist = new TH2D("","",N,-4,0,N,0,3);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            hist->SetBinContent(i+1,j+1,Matrix[i][j]);
        }
    }
    gStyle->SetOptStat(0);  //hide notes
    TAxis *xa = hist->GetXaxis();
    TAxis *ya = hist->GetYaxis();
    
    xa->SetTitle("log_{10}(sin^{2}(2#theta))");
    ya->SetTitle("log_{10}(#Delta m^{2})");
    hist->Draw("surf3");    cvs.Print("5.png");    //good
    hist->Draw("LEGO2Z");   cvs.Print("6.png");    //good
    hist->Draw("COLZ");     cvs.Print("7.png");    //good
    
    
    //do interpolation
    
}

void SinglePoint(double x1,double x2,int samples)
{
    printf("M:\n");
    PhasePoint M;
    M.SetValue(x1,x2);
    M.PrintValue();
    M.FindDchi2(samples);
    printf("dchi2 = %f\n",M.dchi2);
}

int main(int null, char* in[])
{
    //SinglePoint();
    PlotDchi2(10,10000);
    
    

    /*
     double x1=0;
     double x2=0;
     sscanf(in[2],"%lf%lf",&x1,&x2);
     */
    
    return 0;
}
