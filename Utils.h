#include <iostream>

#include <sstream>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TRandom3.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>



void addText(double x1, double x2, double y1, double y2, TString TEXT, Color_t color, Float_t angle = 0)
{
    TPaveText* T = new TPaveText(x1,y1,x2,y2, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    T->SetTextColor(color);
    TText *text = T->AddText(TEXT);
    text->SetTextAngle(angle);
    text->SetTextAlign(22);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);
};

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (double number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

TString Convert2TString (int number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (float number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (double number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TStyle * createTdrStyle() {
    TStyle * tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    //// For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(600); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);

    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    //tdrStyle->SetPadGridX(1);
    //tdrStyle->SetPadGridY(1);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);

    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);

    // For the histo:
    // tdrStyle->SetHistFillColor(1);
    // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    // tdrStyle->SetLegoInnerR(Float_t rad = 0->5);
    // tdrStyle->SetNumberContours(Int_t number = 20);

    tdrStyle->SetEndErrorSize(2);
    //tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(.5);

    tdrStyle->SetMarkerStyle(20);

    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);

    //For the date:
    tdrStyle->SetOptDate(0);
    // tdrStyle->SetDateX(Float_t x = 0->01);
    // tdrStyle->SetDateY(Float_t y = 0->01);

    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(43);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);


    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.16);
    tdrStyle->SetPadRightMargin(0.05);

    // For the Global title:
    tdrStyle->SetTitleFont(43);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);

    // For the axis titles:
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");

    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);

    // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);

    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);


    tdrStyle->cd();
    return tdrStyle;
}

void applyStyleToGraph(TGraph * h_g1){
    h_g1->SetLineWidth(4);
    h_g1->SetMarkerStyle(20);
    h_g1->SetLineStyle(1);
    h_g1->SetMarkerSize(0.8);
    h_g1->SetMarkerColor(2);
    h_g1->SetLineColor(2);
}
TGraph * sortGraph(TGraph * h){
    double x(0), y(0);
    TGraph * h_sorted = new TGraph();

    double last_x(0),last_y(0), this_x(0), this_y(0), new_x(0), new_y(0);
    int n = 0;
    int new_i(0);
    h->GetPoint(0, last_x, last_y );
    h_sorted->SetPoint(n++,last_x,last_y);
    h->RemovePoint(1);

    double distance = 0;
    int last_i=0;
    int npoints = h->GetN();
    while(n < npoints) {
        double min_distance = 9999999;
        for( int i = 0; i<h->GetN();i++) {
            h->GetPoint(i,this_x,this_y);
            distance = sqrt(pow(this_x-last_x,2) + pow(this_y-last_y,2));
            if( distance < min_distance ) {
                min_distance = distance;
                new_x = this_x;
                new_y = this_y;
                new_i = i;
            }
        }
        //~ std::cout << last_x << " " << last_y << " " << min_distance << " " << std::endl;
        h_sorted->SetPoint(n,new_x,new_y);
        last_x = new_x;
        last_y = new_y;
        h->RemovePoint(new_i);
        //~ std::cout << n << std::endl;
        n++;
    }
    return h_sorted;
}
TGraph *InterpolateDM(TString tag, TGraph2D* h_Limit, double sqrt_gxgq=1.0, double xmin=1, double xmax=1000, double ymin=1, double ymax=200)
{
    bool debug = false; 
    bool verboseDebug = false; 

    std::cout << "Start interpolation." << std::endl;
    TGraph *h_ =  new TGraph();
    TRandom3 *r0 = new TRandom3();

    int maxTries = 5000;
    int nnpts=0;
    int Nsteps = 300; //100;
    //int stepsWithoutHit = 0;
    //double x(0.), y(0.),muval(0.);

    int cnt=0; 
    // Let's try "diagonal" interpolation 
    double x0 = 80.; 
    double y0 = 1.; 
    double xL = xmax - x0; 
    double yL = ymax - y0; 
    double alphamin = 0.001; 
    double alphamax = 1.0*TMath::Pi(); 
    double deltaAlpha = (alphamax-alphamin)/Nsteps; 
    bool notyetchanged=true; 
    //for(int j=0; j<Nsteps; j++) { // to be changed, this is just to skip alpha=0 
    for(double alpha=alphamin; alpha<alphamax; alpha+=deltaAlpha) { 
      //x = j*(xmax-xmin)/Nsteps;
      //double alpha = alphamin + j*(alphamax-alphamin)/Nsteps;
      if(notyetchanged && alpha>0.5*TMath::Pi()) { 
	deltaAlpha *= 30.; 
	xL = x0 - xmin; 
	notyetchanged = false; 
      } 
      double xsel(0.), ysel(0.), muvalratiosel(999999.); 
      for(int i=1; i<maxTries; i++) {
	//y=TMath::Exp(r0->Uniform(TMath::Log(ymin),TMath::Log(ymax)));
	//y=r0->Uniform(ymin,ymax);
	//y=ymin + r0->Exp((ymax-ymin)/10);
	//double rho = r0->Exp(10.); 
	double rhoL = std::sqrt( std::pow(xL*std::cos(alpha),2 ) + std::pow(yL*std::sin(alpha), 2) ); 
	double rho = (-1.)*rhoL*std::log10(1-i*1./maxTries); 
	double x = x0 + rho*std::cos(alpha); 
	double y = y0 + rho*std::sin(alpha); 
	if(x>=xmax || y>=ymax || x<=xmin || y<=ymin) continue;  
	double muval = h_Limit->Interpolate(x,y);
	double muvalratio = fabs(muval-sqrt_gxgq)/sqrt_gxgq; 
	//double muval = h_Limit->
	//std::cout << muval << " - "; 
	if(muvalratio<0.01) {
	  if(debug && verboseDebug) std::cout << " --- Candidate: (x, y;   mu) = (" << x << ", " << y << ";   " << muval << ")"; 
	  if(muvalratio<muvalratiosel) {
	    if(debug && verboseDebug) std::cout << "  --->  selected!" << std::endl; 
	    xsel = x; 
	    ysel = y; 
	    muvalratiosel = muvalratio; 
	    //h_->SetPoint(nnpts, x, y);
	    //nnpts++; 
	    //break; 
	  }
	  else {
	    if(debug && verboseDebug) std::cout << std::endl; 
	  }
	}
	//if( i == maxTries - 1 ) stepsWithoutHit++;
      }
      if(muvalratiosel<999990.) {
	h_->SetPoint(nnpts++, xsel, ysel);
	if(debug) std::cout << " *** Selected: (x, y;   mu ratio) = (" << xsel << ", " << ysel << ";   " << muvalratiosel << ")" << std::endl; 
      }
      if(debug && verboseDebug) std::cout << std::endl; 
      //~ if( stepsWithoutHit > 10 ) break;
    }

    /*
    for(int j=0; j<Nsteps; j++) {
        x = j*(xmax-xmin)/Nsteps;
        for(int i=0; i<maxTries; i++) {
            //~ y=TMath::Exp(r0->Uniform(TMath::Log(ymin),TMath::Log(ymax)));
            //~ y=r0->Uniform(ymin,ymax);
            y=ymin + r0->Exp((ymax-ymin)/10);
            muval = h_Limit->Interpolate(x,y);
            //~ std::cout << x << " " << y << std::endl;
            if(fabs(muval-sqrt_gxgq)/sqrt_gxgq<3e-2) {
                h_->SetPoint(nnpts,x,y);
                nnpts++;
            }
            //if( i == maxTries - 1 ) stepsWithoutHit++;
        }
        //~ if( stepsWithoutHit > 10 ) break;
    }

    //stepsWithoutHit = 0;
    for(int j=0; j<Nsteps; j++) {
        y = j*(ymax-ymin)/Nsteps;
        for(int i=0; i<maxTries; i++) {
            //~ x=r0->Uniform(xmin,xmax);
            x=xmin + r0->Exp((xmax-xmin)/10);
            muval = h_Limit->Interpolate(x,y);
            //~ std::cout << x << " " << y << std::endl;
            if(fabs(muval-sqrt_gxgq)/sqrt_gxgq<3e-2) {
                h_->SetPoint(nnpts,x,y);
                nnpts++;
            }
            //if( i == maxTries - 1 ) stepsWithoutHit++;
        }
        //~ if( stepsWithoutHit > 10 ) break;
    }
    */

    applyStyleToGraph(h_);
    std::cout << "Finish interpolation." << std::endl;
    //h_->Sort();

    
    //return sortGraph(h_);
    return h_; 
}



void dumpGraphToFile( TGraph* g, TString filename ) {
    ofstream f;
    f.open(filename);
    double x(0), y(0);
    for( int i=0; i<g->GetN(); i++ ) {
        g->GetPoint(i,x,y);
        f << x << " " << y << endl;
    }
    f.close();
}
