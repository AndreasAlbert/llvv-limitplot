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

void plotDM_EWK_K1K2_mu_obs(TString myfolder = "")
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
            double Lambda_exp = exp;
            double Lambda_obs = obs;
            double Lambda_m1s = m1s;
            double Lambda_m2s = m2s;
            double Lambda_p1s = p1s;
            double Lambda_p2s = p2s;

            h_obs->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_obs);
            h_exp->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_exp);
            h_p1s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_p1s);
            h_m1s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_m1s);

            std::cout << npoint<< " " <<dm_masses[nmx]<<" " <<K1_num[nk1]<<" " <<exp << std::endl;
            npoint++;
        }

    }

    TGraph*   exclusion_exp =  InterpolateDM(tag,h_exp,1,700,1400,0.1,5);
    TGraph*   exclusion_p1s =  InterpolateDM(tag,h_p1s,1,700,1400,0.1,5);
    TGraph*   exclusion_m1s =  InterpolateDM(tag,h_m1s,1,700,1400,0.1,5);
    TGraph*   exclusion_obs =  InterpolateDM(tag,h_obs,1,700,1400,0.1,5);
    exclusion_exp -> SetMarkerColor( kBlack );
    exclusion_p1s -> SetMarkerColor( kGray );
    exclusion_m1s -> SetMarkerColor( kGray );
    exclusion_obs -> SetMarkerColor( kBlue );
    exclusion_exp -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineColor( kBlack );
    exclusion_m1s -> SetLineColor( kBlack );
    exclusion_obs -> SetLineColor( kBlue );
    exclusion_exp -> SetLineWidth(6);
    exclusion_p1s -> SetLineWidth(6);
    exclusion_m1s -> SetLineWidth(6);
    exclusion_obs -> SetLineWidth(6);
    exclusion_p1s -> SetLineStyle( 4 );
    exclusion_m1s -> SetLineStyle( 4 );

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
    t1->SetLogz(true);
    t1->SetRightMargin(0.2);

    TH2D* h2 = new TH2D("h2","",1000,0,1300,100,0.1,10);

    h_obs->SetHistogram(h2);
    h_obs->Draw("COLZ");
    exclusion_exp->Draw("L SAMES");
    exclusion_p1s->Draw("L SAMES");
    exclusion_m1s->Draw("L SAMES");
    exclusion_obs->Draw("L SAMES");
    h_obs->SetMaximum(5);
    h_obs->SetMinimum(5e-3);

    h_obs->GetXaxis()->SetTitle("#it{m_{#chi}} [GeV]");
    h_obs->GetYaxis()->SetTitle("Coupling c_{1}/c_{2}");
    h_obs->GetZaxis()->SetTitle("90% CL observed limit on #sigma_{obs}/#sigma_{theo}");
    h_obs->GetZaxis()->SetRangeUser(5e-3,5);



    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.17,0.37,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Preliminary}}",kBlack);
    addText(0.17,0.47,0.15,0.35,"#Lambda=300 GeV, c_{2}=1",kBlack);



    float posx1 = 0.2;
    float posx2 = 0.4;
    float posy1 = 0.3;
    float posy2 = 0.45;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetLineColor(kBlack);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetHeader("#sigma/#sigma_{theo}=1:");
    leg->AddEntry(exclusion_obs, "Observed", "L");
    leg->AddEntry(exclusion_exp, "Expected", "L");
    leg->AddEntry(exclusion_p1s, "Expected - 1#sigma", "L");

    leg->Draw();



    canv->SaveAs("EWKDM_13TeV_k1k2_mu_obs.png");
    canv->SaveAs("EWKDM_13TeV_k1k2_mu_obs.pdf");

    delete canv;
}
