//
//  plot2D_DMV.cc
//
//
//  Created by RENJIE WANG on 11/6/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TString.h"

#include "Utils.h"

using namespace std;


void plotDM_EWK_1D(TString myfolder = "")
{


    TString tag = "";
    if(myfolder.Contains("/MEWK")) tag="MEWK";

    bool showObserved = true;

    vector<int> dm_masses;
    dm_masses.push_back(1);
    dm_masses.push_back(10);
    dm_masses.push_back(50);
    dm_masses.push_back(100);
    dm_masses.push_back(200);
    dm_masses.push_back(400);
    dm_masses.push_back(800);
    dm_masses.push_back(1300);

    vector<TString> K1;


    K1.push_back("1");


    vector<double> K1_num;

    K1_num.push_back(1);




    TGraph* h_exp = new TGraph();
    TGraphAsymmErrors* h_1s = new TGraphAsymmErrors();
    TGraphAsymmErrors* h_2s = new TGraphAsymmErrors();

    int npoint=0;

    for(int nmx=0; nmx<dm_masses.size(); nmx++) {

        TString str_mx = Convert2TString(dm_masses[nmx]);
        str_mx += "_";
        str_mx += "1";

        TFile myfile( myfolder+"/"+str_mx+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

        Double_t obs;
        Double_t m2s;
        Double_t m1s;
        Double_t exp;
        Double_t p1s;
        Double_t p2s;

        if ( myfile.IsOpen() && !myfile.IsZombie() ) {

            Double_t lim;
            TTree *t = (TTree*)myfile.Get("limit");
            t->SetBranchAddress("limit",    &lim);

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
        myfile.Close();


        //Lambda
        double Lambda_exp = 3000./pow(exp*1.0e6,1./6.);
        double Lambda_obs = 3000./pow(obs*1.0e6,1./6.);
        double Lambda_m1s = 3000./pow(m1s*1.0e6,1./6.);
        double Lambda_m2s = 3000./pow(m2s*1.0e6,1./6.);
        double Lambda_p1s = 3000./pow(p1s*1.0e6,1./6.);
        double Lambda_p2s = 3000./pow(p2s*1.0e6,1./6.);


        Lambda_p2s =  Lambda_p2s - Lambda_exp;
        Lambda_p1s =  Lambda_p1s - Lambda_exp;

        Lambda_m2s =  Lambda_exp - Lambda_m2s;
        Lambda_m1s =  Lambda_exp - Lambda_m1s;

        h_exp ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);
        h_1s    ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);
        h_2s    ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);

        h_1s    ->SetPointError(npoint,0,0,Lambda_p1s,Lambda_m1s);
        h_2s    ->SetPointError(npoint,0,0,Lambda_p2s,Lambda_m2s);

        npoint++;
    }
    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->cd();
    
    std::cout << "Create Canvas and pad." << std::endl;
    TCanvas *canv = new TCanvas("c", "c", 600,600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);

    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogx(false);
    t1->SetLogy(false);
    t1->SetLogz(false);

    h_exp->SetMarkerStyle(24);
    h_exp->SetLineStyle(2);

    h_1s->SetFillColor(kGreen);
    h_2s->SetFillColor(kYellow);

    h_2s->SetMinimum(150);
    h_2s->SetMaximum(900);
    

    h_2s->Draw("3A");
    h_1s->Draw("3,SAME");
    h_exp->Draw("LP,SAME");
    
    h_2s->GetXaxis()->SetTitle( "m_{#chi} [GeV]" );

    h_2s->GetXaxis()->SetRangeUser(-200,1350);
    h_2s->GetYaxis()->SetTitle( "90% CL limit on #Lambda [GeV]" );


    TLegend *leg = new TLegend(0.15,0.63,0.6,0.88);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);

    leg->AddEntry( h_exp, "Expected", "PL" );
    leg->AddEntry( h_1s, "Expected #pm 1#sigma", "F" );
    leg->AddEntry( h_2s, "Expected #pm 2#sigma", "F" );

    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.6,0.8,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);
    addText(0.17,0.37,0.2,0.3,"c_{1} = c_{2} = 1",kBlack);
    leg->Draw();
    canv->SaveAs("EWKDM_13TeV_Lambda.png");
    canv->SaveAs("EWKDM_13TeV_Lambda.pdf");
    
    delete canv;
}
