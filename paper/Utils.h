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



void addText(double x1, double x2, double y1, double y2, TString TEXT, Color_t color, Float_t angle = 0, Float_t tsize=-1.)
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
    if(tsize>0.) T->SetTextSize(tsize);
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
TGraph *InterpolateDM(TString tag, TGraph2D* h_Limit, double sqrt_gxgq=1.0, double xmin=1, double xmax=1000, double ymin=1, double ymax=200, double xctr=100, double yctr=1, int nsteps=100)
{
    bool debug = true;
    bool verboseDebug = false;

    std::cout << "Start interpolation: " << h_Limit->GetName() << std::endl;
    TGraph *h_ =  new TGraph();
    TRandom3 *r0 = new TRandom3();

    //int maxTries = 5000;
    int nnpts=0;
    int Nsteps = nsteps; // 100 for MV-gQ=1.0 and MA-gQ=1.0, 300 for MV-gQ=0.25, 100 for  MA-gQ=0.25

    int cnt=0;
    // Let's try "diagonal" interpolation
    double x0 = xctr; // 200. for MV=gQ=1, 120. for MA-gQ=1, 120. for MV-gQ=0.25
    double y0 = yctr;
    double xpre = -1.;
    double ypre = -1.;
    double deltarpre = 999.;
    double skippedDeltar = 0.;
    //double deltarmax = 2000./Nsteps;
    //double xL = xmax - x0;
    //double yL = ymax - y0;
    //double cosmin = xL/std::sqrt(xL*xL + yL*yL);
    double alphamin = 0.001;
    double alphamax = 1.0*TMath::Pi();
    //double alphamin = 0.8*TMath::Pi();
    //double alphamax = 0.9*TMath::Pi();
    double deltaAlpha = (alphamax-alphamin)/Nsteps;
    double oldDeltaAlpha = deltaAlpha;
    bool notyetchanged_1(true), notyetchanged_2(true), notyetchanged_3(true);

    double muvalratiomax = 0.01;
    double maxstepdist = 1.5;
    double deltarhored = 0.01;

    for(double alpha=alphamin; alpha<alphamax; alpha+=deltaAlpha) {
	double cosa = std::cos(alpha);
	double sina = std::sin(alpha);
	//double rhoL = cosa>cosmin ? xL/cosa : yL/sina;
	double rhomax = 1000.;
	double deltarho = 0.001 * std::sqrt(cosa*cosa*xmax*xmax + sina*sina*ymax*ymax); // Adjust step with angle
	bool islargestep = true;

        double xsel(0.), ysel(0.), muvalratiosel(999999.), deltarsel(999.);
        //for(int i=1; i<maxTries; i++) {
	  //double rhoL = std::sqrt( std::pow(xL*std::cos(alpha),2 ) + std::pow(yL*std::sin(alpha), 2) );
	  //double rho = (-1.)*rhoL*std::log10(1-i*1./maxTries);

        for(double rho=deltarho; rho<rhomax; rho+=deltarho) {
            double x = x0 + rho*std::cos(alpha);
            double y = y0 + rho*std::sin(alpha);
	    double deltar = xpre>0. ? sqrt(std::pow(x-xpre, 2) + std::pow(y-ypre, 2)) : 999.;
            if(x>=xmax || y>=ymax || x<=xmin || y<=ymin) break; // out of boundary
            double muval = h_Limit->Interpolate(x,y);
	    if(muval>(1.+muvalratiomax)*sqrt_gxgq) break;
	    if(islargestep && muval>(1.-2.*muvalratiomax)*sqrt_gxgq) {
	      deltarho *= deltarhored;
	      islargestep = false;
	    }
            double muvalratio = fabs(muval-sqrt_gxgq)/sqrt_gxgq;
            if(muvalratio<muvalratiomax && deltar<maxstepdist*deltarpre) {
	        if(debug && verboseDebug) printf(" --- Candidate: (x, y;  mu; delta-r) = (%7.3f, %7.3f;  %8.6f;  %7.3f [prev. %7.3f])", x, y, muval, deltar, deltarpre);

	        //if( (deltar>998. && muvalratio<muvalratiosel) || ( deltar<998. && deltar<deltarsel) ) {
                if(muvalratio<muvalratiosel) {
		    if(debug && verboseDebug) printf("  --->  selected!\n");
                    xsel = x;
                    ysel = y;
                    muvalratiosel = muvalratio;
		    deltarsel = deltar;
                }
                else {
		    if(debug && verboseDebug) printf("\n");
                }
            }
        } // end for(int i=1; i<maxTries; i++)
        if(muvalratiosel<999990.) {
	    //if(alpha<0.78 || (alpha>0.79 && alpha<0.94) || alpha>0.95) // for MA-gQ=1
	    //if(alpha<0.60 || (alpha>0.61 && alpha<0.65) || alpha>0.67) // for MV-gQ=0.25
	    h_->SetPoint(nnpts++, xsel, ysel);
	    //if(deltarsel<deltarmax || deltarsel>998.) h_->SetPoint(nnpts++, xsel, ysel);
            if(debug) printf(" *** Selected: (x, y;  mu; delta-r) = (%7.3f, %7.3f;  %8.6f;  %7.3f [prev. %7.3f])\n", xsel, ysel, muvalratiosel, deltarsel, deltarpre);
	    xpre = xsel;
	    ypre = ysel;
	    deltarpre = deltarsel;
	    skippedDeltar = 0.;
        }
	else {
	    if(skippedDeltar==0.) skippedDeltar = deltarpre;
	    deltarpre += skippedDeltar;
	}
        if(debug && verboseDebug) std::cout << std::endl;

	// Change step
        if(notyetchanged_1 && alpha<0.3*TMath::Pi()) {
	    notyetchanged_1 = false;
	    double multiplyStep = 0.1;
	    double oldDeltaAlpha = deltaAlpha;
	    deltaAlpha = multiplyStep * ((alphamax-alphamin)/Nsteps);
	    if(deltarpre<998. && deltaAlpha>oldDeltaAlpha) deltarpre *= (deltaAlpha/oldDeltaAlpha);
	    deltarhored = 0.01;
	    muvalratiomax = 0.01;
	    maxstepdist = 1.5;
            //xL = x0 - xmin;
	    //if(xL<=0) break;
	    //cosmin = xL/std::sqrt(xL*xL + yL*yL);
        }
        if(notyetchanged_2 && alpha>=0.3*TMath::Pi()) {
	    notyetchanged_2 = false;
	    double multiplyStep = 2.0;
	    double oldDeltaAlpha = deltaAlpha;
	    deltaAlpha = multiplyStep * ((alphamax-alphamin)/Nsteps);
	    if(deltarpre<998. && deltaAlpha>oldDeltaAlpha) deltarpre *= (deltaAlpha/oldDeltaAlpha);
	    deltarhored = 0.01;
	    muvalratiomax = 0.01;
	    maxstepdist = 2.0;
            //xL = x0 - xmin;
	    //if(xL<=0) break;
	    //cosmin = xL/std::sqrt(xL*xL + yL*yL);
        }
        if(notyetchanged_3 && alpha>0.5*TMath::Pi()) {
            notyetchanged_3 = false;
	    double multiplyStep = 4.0;
	    double oldDeltaAlpha = deltaAlpha;
	    deltaAlpha = multiplyStep * ((alphamax-alphamin)/Nsteps);
	    if(deltarpre<998. && deltaAlpha>oldDeltaAlpha) deltarpre *= (deltaAlpha/oldDeltaAlpha);
	    deltarpre *= (deltaAlpha/oldDeltaAlpha);
	    muvalratiomax = 0.05;
	    maxstepdist = 3.0;
            //xL = x0 - xmin;
	    //if(xL<=0) break;
	    //cosmin = xL/std::sqrt(xL*xL + yL*yL);
        }
    }

    applyStyleToGraph(h_);
    std::cout << "Finish interpolation." << std::endl;

    return h_;
}
TGraph *InterpolateDMXY(TGraph2D* h_Limit, double sqrt_gxgq=1.0, double xmin=1, double xmax=1000, double ymin=1, double ymax=200)
{
    std::cout << "Start interpolation." << std::endl;
    TGraph *h_ =  new TGraph();
    TRandom *r0 = new TRandom();

    int maxTries = 1000;
    int nnpts=0;
    int Nsteps = 100;
    int stepsWithoutHit = 0;
    double x(0.), y(0.),mu_val(0.);

    for(int j=0; j<Nsteps; j++) {
        x = j*(xmax-xmin)/Nsteps;
        for(int i=0; i<maxTries; i++) {
            y=r0->Uniform(ymin,ymax);
            mu_val = h_Limit->Interpolate(x,y);
            if(fabs(mu_val-sqrt_gxgq)/sqrt_gxgq<1e-3) {
                h_->SetPoint(nnpts,x,y);
                nnpts++;
                break;
            }
            if( i == maxTries - 1 ) stepsWithoutHit++;
        }
        //~ if( stepsWithoutHit > 10 ) break;
    }

    stepsWithoutHit = 0;
    //~ maxTries = 100;
    //~ Nsteps = 20;
    for(int j=0; j<Nsteps; j++) {
        y = j*(ymax-ymin)/Nsteps;
        for(int i=0; i<maxTries; i++) {
            x=r0->Uniform(xmin,xmax);
            mu_val = h_Limit->Interpolate(x,y);
            if(fabs(mu_val-sqrt_gxgq)/sqrt_gxgq<1e-3) {
                h_->SetPoint(nnpts,x,y);
                nnpts++;
                break;
            }
            if( i == maxTries - 1 ) stepsWithoutHit++;
        }
        //~ if( stepsWithoutHit > 10 ) break;
    }

    for(int i=0; i<50000; i++) {
        x=r0->Uniform(xmin,xmax);
        mu_val = h_Limit->Interpolate(x,0.1051);
        if(fabs(mu_val-sqrt_gxgq)/sqrt_gxgq<1e-4) {
            h_->SetPoint(nnpts,x,0.1051);
            nnpts++;
            break;
        }
        if( i == maxTries - 1 ) stepsWithoutHit++;
    }

    applyStyleToGraph(h_);
    std::cout << "Finish interpolation." << std::endl;
    h_->Sort(&TGraph::CompareY);
    return h_;
    return sortGraph(h_);
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
void setColorMap() {
    const Int_t     NRGBs         = 3;
    const Int_t     NCont         = 255;
    Double_t stops[NRGBs] = { 0.00, 0.2,1.0 };
    //~ Double_t red[NRGBs]   = { 23, 56,229  };
    //~ Double_t green[NRGBs] = { 77, 128,245  };
    //~ Double_t blue[NRGBs]  = { 40, 48,224 };
    Double_t red[NRGBs]   = { 23, 56,255  };
    Double_t green[NRGBs] = { 77, 128,255  };
    Double_t blue[NRGBs]  = { 40, 48,255 };
    for(int i = 0; i<NRGBs; i++) {
        red[i] = red[i]/255;
        green[i] = green[i]/255;
        blue[i] = blue[i]/255;
    }
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);
}
int getRelicColor() {
    return kAzure+3;
}
int getDiagonalColor() {
    return kAzure+3;
}
int getObservedColor() {
    return kRed-4;
}
double get_limit_from_file( TString file_path ) {

    double limit=-1;
    TFile myfile( file_path );
    if( myfile.IsOpen()==true && myfile.IsZombie()==false) {
        TTree *t = (TTree*)myfile.Get("limit");
        if( t != 0 ) {
            t->SetBranchAddress("limit",    & limit);
            t->GetEntry(0);
        }
        myfile.Close();
    }
    std::cout << file_path << " " << limit << std::endl;
    return limit;
}
