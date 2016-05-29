#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TString.h"
#include <math.h>
#include <iomanip>
using namespace std;

#include "Utils.h"

void drawLine(double x, double xx){

    TLine *ll;
    ll = new TLine(x, xx, x, 3.0e-4);
    ll->SetLineWidth(2);
    ll->SetLineColor(kGray+2);
    ll->SetLineStyle(3);
    ll->Draw();
}
void printTexLine( std::string preface, std::vector<double> vec ) {
    std::setw(3);
    std::setprecision(3);
    std::cout.precision(3);

    std::cout << preface;
    for( auto &v : vec ) {
        std::cout << " & ";
        std::cout << v;
    }
    std::cout << " \\\\" << std::endl;
}
void plotModelIndepXS(TString myfolder)
{

    double lumi = 2300;
    
    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->cd();
    
    gStyle->SetOptStat(0);
    bool is8TeV = true;
    bool savePlots = true;
    bool showObserved = true;
    gStyle->SetPadTopMargin(0.06);



    TCanvas *canv = new TCanvas("canv", "limits canvas",600., 600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);
    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    //t1->SetLogx(true);

    //TMultiGraph *mg = new TMultiGraph();

    std::vector<double> x_v;
    std::vector<double> exp_v;
    std::vector<double> p1s_v;
    std::vector<double> p2s_v;
    std::vector<double> m1s_v;
    std::vector<double> m2s_v;

    TGraphAsymmErrors * gr1s = new TGraphAsymmErrors();
    TGraphAsymmErrors * gr2s = new TGraphAsymmErrors();
    TGraphAsymmErrors * grExp = new TGraphAsymmErrors();

    int npoints = 0;
    for( int i=1; i<=13; i+=2) {
        TFile myfile( myfolder+"/"+Convert2TString( i )+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

        Double_t obs;
        Double_t m2s;
        Double_t m1s;
        Double_t exp;
        Double_t p1s;
        Double_t p2s;

        if ( myfile.IsOpen() && !myfile.IsZombie() ) {

            

            Double_t lim, limerr;
            TTree *t = (TTree*)myfile.Get("limit");
            t->SetBranchAddress("limit",    &lim);
            t->SetBranchAddress("limitErr", &limerr);

            // Expected  2.5%
            t->GetEntry(0);
            m2s  = lim;

            // Expected  16.0%
            t->GetEntry(1);
            m1s = lim;

            // Expected 50.0%
            t->GetEntry(2);
            exp = lim;

            // Expected 84.0%
            t->GetEntry(3);
            p1s = lim;

            // Expected 97.5%
            t->GetEntry(4);
            p2s = lim;

            // Observed limit
            t->GetEntry(5);
            obs = lim;

            // The limits are limits on signal strengths
            // We use a signal with 10 events, so we need to multiply by 10
            grExp->SetPoint( npoints, 70+10*i, 10 * exp / lumi );
            gr1s->SetPoint(  npoints, 70+10*i, 10 * exp / lumi );
            gr2s->SetPoint(  npoints, 70+10*i, 10 * exp / lumi );
            grExp->SetPointError(npoints,0,0,0,0);
            gr1s->SetPointError(npoints,0,0,10 * (exp-m1s) / lumi ,10 * (p1s-exp) / lumi );
            gr2s->SetPointError(npoints,0,0,10 * (exp-m2s) /  lumi ,10 * (p2s-exp) / lumi );

            //~ std::cout << (exp-m1s) / lumi << " " << exp / lumi << " " << (p1s-exp) / lumi << std::endl;
            npoints++;

            x_v.push_back( 70 + 10 * i );
            exp_v.push_back( 10 * exp );
            p1s_v.push_back( 10 * p1s );
            p2s_v.push_back( 10 * p2s );
            m1s_v.push_back( 10 * m1s );
            m2s_v.push_back( 10 * m2s );
        }


        myfile.Close();
    }


    gr1s->SetFillColor(kGreen);

    gr2s->SetFillColor(kYellow);

    grExp->SetLineWidth(2);
    grExp->SetLineStyle(2);
    grExp->SetMarkerStyle(24);



    TMultiGraph *mg = new TMultiGraph();

    mg->Add(gr2s,"3");
    mg->Add(gr1s,"3");
    //mg->Add(grObs,"P");
    mg->Add(grExp,"P");
//~ 
    mg->Draw("3A"); //RJ




    mg->SetMaximum(1e-1);
    mg->SetMinimum(1e-3);


    mg->GetYaxis()->SetTitle("#sigma^{BSM}_{vis} = #sigma A #varepsilon [pb]");
    mg->GetXaxis()->SetTitle("E_{T}^{miss} threshold [GeV]");
    mg->GetXaxis()->SetNdivisions(10 + 100*0 + 10000*0);

    float posx1 = 0.2+0.4;
    float posx2 = 0.48+0.45;
    float posy1 = 0.65;
    float posy2 = 0.92;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);

    //leg->AddEntry(grObs, "Observed", "PL");
    leg->AddEntry(grExp, "Expected", "P");
    leg->AddEntry(gr1s, "Expected #pm 1#sigma", "F");
    leg->AddEntry(gr2s, "Expected #pm 2#sigma", "F");
    leg->Draw();




    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);

    addText(0.7-0.02,0.995-0.02,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);

    //~ T = new TPaveText(0.50+0.01-0.025,0.97,0.68+0.01+0.025,0.82, "NDC");

    addText(0.1+0.06+0.05,0.3+0.2+0.05,0.835+0.02,0.898+0.02,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);
    //~ addText(0.20,0.45+0.21,0.83-0.05+0.02,0.87-0.02+0.02,"#it{pp} #rightarrow Z + E^{miss}_{T} #rightarrow #it{l^{+}l^{-}} + E^{miss}_{T}",kBlack);


    t1->RedrawAxis();
    if( savePlots ) {
        std::string plotName = "Zwimps";
        if(is8TeV) plotName += "_13TeV";
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS.png").c_str() );
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS.pdf").c_str() );
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS.eps").c_str() );
    }
    printTexLine( "\\ETm (\\GeVns{}) threshold",      x_v );
    printTexLine( "Exp. upper limit +2$\\sigma$",    p2s_v );
    printTexLine( "Exp. upper limit +1$\\sigma$",    p1s_v );
    printTexLine( "Exp. upper limit",               exp_v );
    printTexLine( "Exp. upper limit -1$\\sigma$",    m1s_v );
    printTexLine( "Exp. upper limit -2$\\sigma$",    m2s_v );




    
    //return;
}



