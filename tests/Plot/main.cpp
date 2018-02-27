#include <iostream>
#include <vector>
#include <cmath>
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;

double z(double x,double y)
{
    double z = pow(x,3) + pow(y,2);
    return z;
}

int main()
{
    TCanvas cvs;
    cvs.SetWindowSize(200,200);
    int N = 20;
    TH2D* hist = new TH2D("","",N,-1,1,N,-1,1);
    for (int i=0;i<N;i++){
        double x = 2*double(i)/N-1;
        for (int j=0;j<N;j++){
            double y = 2*double(j)/N-1;
            hist->SetBinContent(i+1,j+1,z(x,y));
        }
    }
    gStyle->SetOptStat(0);  //hide notes
    
    TAxis* xa = hist->GetXaxis();
    TAxis* ya = hist->GetYaxis();
    xa->SetTitle("x");
    ya->SetTitle("y");
    hist->Draw("CONTZ"); cvs.Print("1.png");
    hist->Draw("CONT3"); cvs.Print("2.png");    //contour
    hist->Draw("LEGO"); cvs.Print("3.png");    //good~
    hist->Draw("surf2"); cvs.Print("4.png");    //the most beautiful fig :)
    hist->Draw("surf3"); cvs.Print("5.png");    //ref
    hist->Draw("LEGO2Z"); cvs.Print("6.png");    //good~
    hist->Draw("COLZ"); cvs.Print("7.png");    //good~


    return 0;
}
