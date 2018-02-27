//	The main function for the simulation of 2 dimensional neutrino
//	oscillation.
#include <iostream>
#include <vector>
#include <cmath>
#include "TF2.h"
using namespace std;
double sin2 = pow(10,-1.33);
double dm2  = pow(10,2);
// x:L y:E
// L[0.6 1] E[10 20]
// (10^-1.33,10^2) -> 0.0939963
int main(){
    TF2 * func1 = new TF2("func1","[0]*sin(1.27*[1]*x/y)**2");
    func1->SetParameters(sin2,dm2);
    double result = func1->Integral(0.6,1,30,40,1e-7);
    cout<<result<<endl;
    return 0;
    
}

/*
 #include <iostream>
 #include <vector>
 #include <cmath>
 #include <time.h>
 #include "TF2.h"
 using namespace std;
 double sin2 = pow(10,-1.33);
 double dm2  = pow(10,2);
 // x:L y:E
 // L[0.6 1] E[10 20]
 int main(){
 TF2 * func1 = new TF2("func1","[0]*sin(1.27*[1]*x/y)**2");
 func1->SetParameters(sin2,dm2);
 time_t start, end;
 time(&start);
 for (int i=0;i<100;i++)
 {   double result = func1->Integral(0.6,1,10,20,1e-7);
 cout<<result<<endl;
 }
 time(&end);
 double cost;
 cost = difftime(end,start);
 cout<<cost<<endl;
 //delete func1;
 return 0;
 
 }

 */
