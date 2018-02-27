//	The main function for the simulation of 2 dimensional neutrino
//	oscillation.
#include <iostream>
#include <vector>
using namespace std;
void PrintVector(vector<double> a);                 //show vector
vector<double> multi(vector<double> a, double x);   //multi in Vectors
static vector<double> I(5,1);

class PhasePoint
{   //  emmm...the declaration of class must be on the top of the code.
    //    (x1,x2) or (sin2, dm2) is the position of phase point
    //    (sin(2*theta)^2,dm^2)
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
    P = multi(I,0);
    muTrue = multi(I,100);
}

void PhasePoint::PrintValue()   {printf("sin(2*theta)^2 = %f, dm^2 = %f\n",sin2,dm2);}

void PhasePoint::FindDchi2()
{
    //functions declaration
    vector<double> FastIntegral(double x1, double x2);
    //    P = FastIntegral(sin2,dm2);
    muTrue = multi(P,10000);
    PrintVector(muTrue);
    /*
     in 10000 measures,
     find n;
     for each n:
     calculate mubest;
     calculate likelihood for mubest and mu;
     calculate Ratio;
     plot a histogram for Ratio;
     find dchi2 of this phase point.
     */
}

/*
 vector<double> FastIntegral(double x1, double x2)
 {
 ///////////////////
 return P;
 }
 */

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


int main()
{
    PhasePoint M;
    M.SetValue(1,3);
    M.PrintValue();
    
    PrintVector(M.P);
    PrintVector(M.muTrue);
    
    M.FindDchi2();
    PrintVector(M.muTrue);

}






