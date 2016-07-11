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

using namespace std;

#include "Utils.h"

double monoZ_X[17]={1.01,1.02,1.04,1.06,1.09,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.20001};

double monoJetX[5] = {1.5,1.6,1.7,1.8,1.9};
double monoJetY[5] = {10.00,4.91,2.91,2.01,1.60};

double LEPX[] = {1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,1.95001};
double LEPY[] = {9.19,6.82,5.18,4.03,3.19,2.58,2.11,1.75,1.48,1.26,1.08,0.94,0.82,0.73,0.64,0.58,0.52,0.47,0.43,0.43};

double monoZ_8TeV_X[] = {1.01,1.02,1.04,1.06,1.09,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.200001};
double monoZ_8TeV_Y[] = {1.65789e+136,3.86191e+74,1.21707e+39,5.83507e+26,6.1486e+17,1.3288e+16,3.2407e+07,33777.3,955.665,114.949,30.8855,10.6683,5.25326,2.7238,1.84864,0.942391,0.942391};




std::vector<int> du = {101,102,104,106,109,110,120,130,140,150,160,170,180,190,200,220};
std::vector<double> r_exp;
void plotUnpart(TString myfolder = "")
{
    TGraph *monoZ_obs = new TGraph();
    TGraph *monoZ_exp = new TGraph();
    TGraphAsymmErrors *monoZ_1s = new TGraphAsymmErrors();
    TGraphAsymmErrors *monoZ_2s = new TGraphAsymmErrors();
    int npoint=0;
    for ( auto & idu: du ) {
        TFile myfile( myfolder+"/Unpart_"+Convert2TString(idu)+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

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

            }
            r_exp.push_back(exp);

            double lambda_exp = 15/pow(exp,100.0/idu);
            double lambda_p1s = 15/pow(p1s,100.0/idu);
            double lambda_m1s = 15/pow(m1s,100.0/idu);
            double lambda_p2s = 15/pow(p2s,100.0/idu);
            double lambda_m2s = 15/pow(m2s,100.0/idu);
            double lambda_obs = 15/pow(obs,100.0/idu);
            std::cout << lambda_exp <<" " << lambda_p1s << " " <<lambda_m1s << std::endl;
            monoZ_exp->SetPoint(npoint,float(idu)/100,lambda_exp);
            monoZ_1s->SetPoint(npoint,float(idu)/100,lambda_exp);
            monoZ_2s->SetPoint(npoint,float(idu)/100,lambda_exp);

            monoZ_1s->SetPointError(npoint,0,0,lambda_exp-lambda_p1s, lambda_m1s-lambda_exp);
            monoZ_2s->SetPointError(npoint,0,0,lambda_exp-lambda_p2s, lambda_m2s-lambda_exp);

            monoZ_obs->SetPoint(npoint,float(idu)/100,lambda_obs);
            npoint++;
            myfile.Close();
    }
    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->cd();

    string parName = "Scaling dimension #it{d}_{U}";
    bool is8TeV = true;
    bool savePlots = true;
    bool showObserved = false;
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin (0.03);

    gStyle->SetHatchesLineWidth(20);


    TCanvas *canv = new TCanvas("canv", "limits canvas", 600.,600.);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);
    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    t1->SetLogx(false);

    TMultiGraph *mg = new TMultiGraph();


    TGraph *monoZ_exp2 = (TGraph*) monoZ_exp->Clone("monoZ_exp2");

    monoZ_exp->SetMarkerStyle(24);
    monoZ_exp->SetLineStyle(2);

    monoZ_1s->SetFillColor(kGreen);
    monoZ_2s->SetFillColor(kYellow);

    monoZ_obs->SetLineColor(kBlack);
    monoZ_obs->SetLineWidth(2);
    monoZ_obs->SetMarkerStyle(20);


    TGraph *monoJet = new TGraph(sizeof(monoJetX)/sizeof(*monoJetX),monoJetX,monoJetY);
    monoJet->SetLineWidth(-603);
    monoJet->SetFillStyle(3005);

    monoJet->SetFillColor(kViolet-5);
    monoJet->SetLineColor(kViolet-6);

    TGraph *monoZ_8TeV = new TGraph(sizeof(monoZ_8TeV_X)/sizeof(*monoZ_8TeV_X),monoZ_8TeV_X,monoZ_8TeV_Y);
    monoZ_8TeV->SetLineWidth(-603);
    monoZ_8TeV->SetFillStyle(3005);
    monoZ_8TeV->SetFillColor(kBlue-5);
    monoZ_8TeV->SetLineColor(kBlue-6);

    TGraph *lep = new TGraph(sizeof(LEPX)/sizeof(*LEPX),LEPX,LEPY);
    TGraph *lep2 = (TGraph*) lep->Clone("lep2");
    lep2->SetLineWidth(-803);
    lep2->SetFillStyle(3005);
    lep2->SetLineColor(kOrange+3);

    lep->SetLineWidth(-402);

    lep->SetFillColor(kOrange+7);
    lep->SetLineColor(kOrange+3);
    lep->SetLineWidth(-9003);

    mg->Add(monoZ_2s,"3");
    mg->Add(monoZ_1s,"3");
    mg->Add(monoZ_exp,"PL");

    mg->Add(monoZ_obs,"PL");
    mg->Add(monoJet);
    mg->Add(lep);
    mg->Add(lep2);
    mg->Add(monoZ_8TeV);
    mg->Draw("AC");
    mg->SetMinimum(0.3);
    mg->SetMaximum(130);
    mg->GetXaxis()->SetNdivisions(512);
    mg->GetXaxis()->SetTitle(parName.c_str());
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitle("#it{#Lambda}_{U} [TeV]");
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(1.35);
    mg->GetXaxis()->SetRangeUser(1.06,2.2-0.01);
    mg->GetXaxis()->SetRangeUser(1.06,2.25);

    addText(0.7-0.02,0.995-0.02,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);

    addText(0.64,0.96,0.835+0.02,0.898+0.02,"#bf{CMS} #it{Preliminary}",kBlack);
    addText(0.6,0.96,0.8,0.87,"#it{pp} #rightarrow Z#it{U} #rightarrow #it{l^{+}l^{-}U}",kBlack);
    addText(0.68,0.96,0.75,0.82,"spin = 0, #it{#lambda} = 1",kBlack);


    float posx1 = 0.65;
    float posx2 = 1;
    float posy1 = 0.56;
    float posy2 = 0.75;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->AddEntry(monoZ_obs, "Observed", "PL");
    leg->AddEntry(monoZ_exp, "Expected", "PCL");
    leg->AddEntry(monoZ_1s, "Expected #pm1#sigma", "F");
    leg->AddEntry(monoZ_2s, "Expected #pm2#sigma", "F");
    leg->AddEntry(monoJet, "CMS Mono-Jet (8TeV)", "FL");
    leg->AddEntry(monoZ_8TeV, "CMS Mono-Z (8 TeV)", "FL");
    leg->AddEntry(lep, "LEP reinterpretation", "FL");
    leg->Draw();


    t1->RedrawAxis();
    if( savePlots ) {
        string plotName = parName;
        string toReplace[] = {"_", "{", "}", "#", "^"};
        UInt_t nvars = sizeof(toReplace)/sizeof(string);
        for(UInt_t k=0; k<nvars; ++k) {
            int poschar = plotName.find(toReplace[k].c_str());
            while( poschar>-1 ) {
                plotName.replace(poschar, 1, "");
                poschar = plotName.find(toReplace[k].c_str());
            }
        }
        plotName = "Unpart";
        if(is8TeV) plotName += "_13TeV";
        canv->SaveAs( ("./"+plotName+"_LambdaU.png").c_str() );
        canv->SaveAs( ("./"+plotName+"_LambdaU.pdf").c_str() );
        canv->SaveAs( ("./"+plotName+"_LambdaU.eps").c_str() );
    }

}



