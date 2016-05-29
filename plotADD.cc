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

double XS[7] ={3.469e-02,3.323e-02,3.904e-02,4.828e-02,6.511e-02,9.103e-02,1.502e-01};

void plotADD(TString myfolder = "")
{


    TString tag = "";
    if(myfolder.Contains("/ADD")) tag="ADD";

    bool showObserved = true;

    vector<int> dm_masses;
    dm_masses.push_back(2);
    dm_masses.push_back(3);
    dm_masses.push_back(4);
    dm_masses.push_back(5);
    dm_masses.push_back(6);
    dm_masses.push_back(7);
    dm_masses.push_back(8);



    vector<double> obs_v;
    vector<double> m2s_v;
    vector<double> m1s_v;
    vector<double> exp_v;
    vector<double> p1s_v;
    vector<double> p2s_v;




    int npoint=0;
    for(int nmx=0; nmx<dm_masses.size(); nmx++) {

        TString str_mx = Convert2TString(dm_masses[nmx]);

        TFile myfile( myfolder+"/"+str_mx+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

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
        myfile.Close();



        double Lambda_exp = exp * XS[npoint];
        double Lambda_obs = obs * XS[npoint];
        double Lambda_m1s = m1s * XS[npoint];
        double Lambda_m2s = m2s * XS[npoint];
        double Lambda_p1s = p1s * XS[npoint];
        double Lambda_p2s = p2s * XS[npoint];

        Lambda_p2s =  Lambda_p2s - Lambda_exp;
        Lambda_p1s =  Lambda_p1s - Lambda_exp;

        Lambda_m2s =  Lambda_exp - Lambda_m2s;
        Lambda_m1s =  Lambda_exp - Lambda_m1s;

        m2s_v.push_back(Lambda_m2s);
        m1s_v.push_back(Lambda_m1s);
        exp_v.push_back(Lambda_exp);
        p1s_v.push_back(Lambda_p1s);
        p2s_v.push_back(Lambda_p2s);
        obs_v.push_back(Lambda_obs);


        npoint++;
    }

    
    TGraph * h_exp = new TGraph( );
    TGraphAsymmErrors * h_1s = new TGraphAsymmErrors( );
    TGraphAsymmErrors * h_2s = new TGraphAsymmErrors( );
    for( int i=0; i<dm_masses.size(); i++ ){
        h_exp->SetPoint(i,dm_masses[i],exp_v[i]);
        h_1s->SetPoint(i,dm_masses[i],exp_v[i]);
        h_2s->SetPoint(i,dm_masses[i],exp_v[i]);
        
        h_1s->SetPointError(i,0,0,m1s_v[i],p1s_v[i]);
        h_2s->SetPointError(i,0,0,m2s_v[i],p2s_v[i]);
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
    t1->SetLogy(true);
    t1->SetLogz(false);


    h_exp->SetMarkerStyle(24);
    h_exp->SetLineStyle(2);

    h_1s->SetFillColor(kGreen);
    h_2s->SetFillColor(kYellow);

    h_2s->SetMinimum(1e-3);
    h_2s->SetMaximum(10);
    

    h_2s->Draw("3A");
    h_1s->Draw("3,SAME");
    h_exp->Draw("LP,SAME");
    
    h_2s->GetXaxis()->SetTitle( "n_{D}" );
    h_2s->GetYaxis()->SetTitle( "95% CL limit on #sigma(pp#rightarrow ZG)#times BR [pb]" );


    TLegend *leg = new TLegend(0.15,0.6,0.6,0.85);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);

    leg->AddEntry( h_exp, "Expected", "PL" );
    leg->AddEntry( h_1s, "Expected #pm 1#sigma", "F" );
    leg->AddEntry( h_2s, "Expected #pm 2#sigma", "F" );

    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.6,0.85,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);
    leg->Draw();

    canv->SaveAs("ADD_13TeV_XS.png");
    canv->SaveAs("ADD_13TeV_XS.pdf");



}
