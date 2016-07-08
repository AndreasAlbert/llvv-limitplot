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


void plotDM_EWK_K1K2_obs(TString myfolder = "")
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

    K1.push_back("0.1");
    K1.push_back("0.2");
    K1.push_back("0.3");
    K1.push_back("0.5");
    K1.push_back("1");
    K1.push_back("2");
    K1.push_back("3");
    K1.push_back("5");
    K1.push_back("10");

    vector<double> K1_num;
    K1_num.push_back(0.1);
    K1_num.push_back(0.2);
    K1_num.push_back(0.3);
    K1_num.push_back(0.5);
    K1_num.push_back(1);
    K1_num.push_back(2);
    K1_num.push_back(3);
    K1_num.push_back(5);
    K1_num.push_back(10);



    TGraph2D* h_obs = new TGraph2D();
    TGraph2D* h_exp = new TGraph2D();
    TGraph2D* h_p1s = new TGraph2D();
    TGraph2D* h_m1s = new TGraph2D();

    int npoint=0;

    for(int nk1=0; nk1<K1.size(); nk1++) {

        for(int nmx=0; nmx<dm_masses.size(); nmx++) {

            TString str_mx = Convert2TString(dm_masses[nmx]);
            str_mx += "_";
            str_mx += K1[nk1];

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

            h_obs->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_obs);
            h_exp->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_exp);
            h_p1s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_p1s);
            h_m1s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_m1s);

            std::cout << npoint<< " " <<dm_masses[nmx]<<" " <<K1_num[nk1]<<" " <<exp << std::endl;
            npoint++;
        }

    }
    TGraph*   exclusion_exp =  InterpolateDM(tag,h_exp,300,800,1400,0.1,5);
    TGraph*   exclusion_p1s =  InterpolateDM(tag,h_p1s,300,800,1400,0.1,5);
    TGraph*   exclusion_m1s =  InterpolateDM(tag,h_m1s,300,800,1400,0.1,5);
    TGraph*   exclusion_obs =  InterpolateDM(tag,h_obs,300,800,1400,0.1,5);
    exclusion_exp -> SetMarkerColor( kBlack );
    exclusion_p1s -> SetMarkerColor( kGray );
    exclusion_m1s -> SetMarkerColor( kGray );
    exclusion_obs -> SetMarkerColor( kRed );
    exclusion_exp -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineColor( kGray );
    exclusion_m1s -> SetLineColor( kGray );
    exclusion_obs -> SetLineColor( kRed );


    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->cd();
    gStyle->SetPalette(55);
    std::cout << "Create Canvas and pad." << std::endl;
    TCanvas *canv = new TCanvas("c", "c", 700, 550);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);

    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogx(false);
    t1->SetLogy(true);
    t1->SetLogz(false);
    t1->SetRightMargin(0.2);


    TH2D* h2 = new TH2D("h2","",1000,0,1300,100,0.1,10);

    h_obs->SetHistogram(h2);
    h_obs->Draw("COLZ");

    h_obs->SetMaximum(6e2);
    h_obs->SetMinimum(2e2);

    h_obs->GetXaxis()->SetTitle("#it{m_{#chi}} [GeV]");
    h_obs->GetYaxis()->SetTitle("Coupling c_{1}/c_{2}");
    h_obs->GetZaxis()->SetTitle("90% CL observed limit on #Lambda [GeV]");
    h_obs->GetZaxis()->SetRangeUser(2e2,6e2);


    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.17,0.37,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Preliminary}}",kBlack);
    addText(0.17,0.37,0.2,0.3,"c_{2}=1",kBlack);
    canv->SaveAs("EWKDM_13TeV_k1k2_obs.png");
    canv->SaveAs("EWKDM_13TeV_k1k2_obs.pdf");

    delete canv;
}
