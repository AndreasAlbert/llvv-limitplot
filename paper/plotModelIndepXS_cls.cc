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
void plotModelIndepXS_cls(TString myfolder)
{

    double lumi = 2.300;

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

    std::vector<double> x_v;
    std::vector<double> obs_v;
    std::vector<double> exp_v;
    std::vector<double> p1s_v;
    std::vector<double> p2s_v;
    std::vector<double> m1s_v;
    std::vector<double> m2s_v;

    TGraphAsymmErrors * gr1s = new TGraphAsymmErrors();
    TGraphAsymmErrors * gr2s = new TGraphAsymmErrors();
    TGraphAsymmErrors * grExp = new TGraphAsymmErrors();
    TGraphAsymmErrors * grObs = new TGraphAsymmErrors();

    int npoints = 0;
    for( int i=1; i<=24; i+=1) {
        TString folder = myfolder+"/ModelIndep_"+Convert2TString( 70+10*i );
        Double_t obs = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.root" );
        Double_t m1s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.160.root" );
        Double_t m2s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.025.root" );
        Double_t exp = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.500.root" );
        Double_t p1s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.840.root" );
        Double_t p2s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.975.root" );

        // The limits are limits on signal strengths
        // We use a signal with 10 events, so we need to multiply by 10

        double threshold = 0;
        if(i<=13) {
            threshold = 70 + 10*i;
            if( i % 2 != 0 ) continue;
        }
        else{
            threshold = 250 + (i-14)*50;
            if( threshold != 300
            and threshold != 450
            and threshold != 600 ) {
            continue;
            }
        }


        grObs->SetPoint( npoints, threshold, 10 * obs / lumi );
        grExp->SetPoint( npoints, threshold, 10 * exp / lumi );
        gr1s->SetPoint(  npoints, threshold, 10 * exp / lumi );
        gr2s->SetPoint(  npoints, threshold, 10 * exp / lumi );

        grObs->SetPointError(npoints,0,0,0,0);
        grExp->SetPointError(npoints,0,0,0,0);

        double err_up1   = fabs(std::max(p1s,m1s) - exp);
        double err_down1 = fabs(std::min(p1s,m1s) - exp);
        double err_up2   = fabs(std::max(p2s,m2s) - exp);
        double err_down2 = fabs(std::min(p2s,m2s) - exp);
        gr1s->SetPointError(npoints,0,0,10 * err_down1 / lumi ,10 * err_up1 / lumi );
        gr2s->SetPointError(npoints,0,0,10 * err_down2 /  lumi ,10 * err_up2 / lumi );

        npoints++;

        x_v.push_back( threshold );
        obs_v.push_back( 10 * obs );
        exp_v.push_back( 10 * exp );
        p1s_v.push_back( 10 * p1s );
        p2s_v.push_back( 10 * p2s );
        m1s_v.push_back( 10 * m1s );
        m2s_v.push_back( 10 * m2s );
    }



    gr1s->SetFillColor(kGreen);

    gr2s->SetFillColor(kYellow);

    grExp->SetLineWidth(2);
    grExp->SetLineStyle(2);
    grExp->SetMarkerStyle(24);
    grExp->SetLineWidth(1);
    grExp->SetLineStyle(2);
    grExp->SetMarkerStyle(24);



    TMultiGraph *mg = new TMultiGraph();

    mg->Add(gr2s,"3");
    mg->Add(gr1s,"3");
    mg->Add(grObs,"P");
    mg->Add(grExp,"LP");
    mg->Draw("3A");




    mg->SetMaximum(1e2);
    mg->SetMinimum(1);


    mg->GetYaxis()->SetTitle("95% CL limit on #sigma^{BSM}_{vis} = #sigma A #varepsilon [fb]");
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

    leg->AddEntry(grObs, "Observed", "LP");
    leg->AddEntry(grExp, "Expected", "LP");
    leg->AddEntry(gr1s, "Expected #pm 1 s.d.", "F");
    leg->AddEntry(gr2s, "Expected #pm 2 s.d.", "F");
    leg->Draw();




    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);

    addText(0.7-0.04,0.995-0.04,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);


    addText(0.,0.45,0.94,1.0,"#bf{CMS}",kBlack);


    t1->RedrawAxis();
    if( savePlots ) {
        std::string plotName = "Zwimps";
        if(is8TeV) plotName += "_13TeV";
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS_cls.png").c_str() );
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS_cls.pdf").c_str() );
        canv->SaveAs( ("./"+plotName+"_ModelIndepXS_cls.eps").c_str() );
    }
    printTexLine( "\\ETm (\\GeVns{}) threshold",      x_v );
    printTexLine( "Exp. upper limit +2$\\sigma$",    p2s_v );
    printTexLine( "Exp. upper limit +1$\\sigma$",    p1s_v );
    printTexLine( "Exp. upper limit",               exp_v );
    printTexLine( "Exp. upper limit -1$\\sigma$",    m1s_v );
    printTexLine( "Exp. upper limit -2$\\sigma$",    m2s_v );
    printTexLine( "Obs. upper limit",    obs_v );




    //return;
}



