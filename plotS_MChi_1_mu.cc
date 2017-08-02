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

double XS[9] ={0.007327882113186,0.007079719965470,0.004326010725432,0.002530922638544,0.001845383347520,0.001122081502664,0.000421147237101,0.000079522094909,0.000021179814444};//,1.502e-01};

void plotS_MChi_1_mu(TString myfolder = "")
{


    TString tag = "";
    if(myfolder.Contains("/NLOS")) tag="NLOS";

    bool showObserved = true;

    vector<int> dm_masses;
    dm_masses.push_back(20);
    dm_masses.push_back(50);
    dm_masses.push_back(200);
    dm_masses.push_back(300);
    dm_masses.push_back(350);
    dm_masses.push_back(400);
    dm_masses.push_back(500);
    dm_masses.push_back(750);
    dm_masses.push_back(1000);    

    //vector<double> obs_v;
    vector<double> m2s_v;
    vector<double> m1s_v;
    vector<double> exp_v;
    vector<double> p1s_v;
    vector<double> p2s_v;
    vector<double> obs_v;
   // vector<double> theoX_v;


    int npoint=0;
    for(int nmx=0; nmx<dm_masses.size(); nmx++) {

        TString str_mv = Convert2TString(dm_masses[nmx]);

        TFile myfile( myfolder+"/1_"+str_mv+"/higgsCombineTest.Asymptotic.mH120.root" );

        Double_t obs;
        Double_t m2s;
        Double_t m1s;
        Double_t exp;
        Double_t p1s;
        Double_t p2s;
        //Double_t theoX;
        
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



        double Lambda_exp = exp ;
        double Lambda_obs = obs ;
        double Lambda_m1s = m1s ;
        double Lambda_m2s = m2s ;
        double Lambda_p1s = p1s ;
        double Lambda_p2s = p2s ;
       // double Lambda_theo = XS[npoint];
        
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
        //theoX_v.push_back(Lambda_theo);

        npoint++;
    }

    
    TGraph * h_exp = new TGraph( );
    TGraph * h_obs = new TGraph( );
    TGraph * h_theo = new TGraph( );
    TGraphAsymmErrors * h_1s = new TGraphAsymmErrors( );
    TGraphAsymmErrors * h_2s = new TGraphAsymmErrors( );
    for( int i=0; i<dm_masses.size(); i++ ){
        h_exp->SetPoint(i,dm_masses[i],exp_v[i]);
      //  h_theo->SetPoint(i,dm_masses[i],theoX_v[i]);
        h_1s->SetPoint(i,dm_masses[i],exp_v[i]);
        h_2s->SetPoint(i,dm_masses[i],exp_v[i]);
        h_obs->SetPoint(i,dm_masses[i],obs_v[i]);
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

    h_obs->SetLineColor(kBlack);
    h_obs->SetLineWidth(2);
    h_obs->SetMarkerStyle(20);

  /*  h_theo->SetLineColor(kRed);
    h_theo->SetLineWidth(2);
    h_theo->SetMarkerStyle(20);*/
    
    h_1s->SetFillColor(kGreen+1);
    h_2s->SetFillColor(kOrange);

    h_2s->SetMinimum(1.e-1);
    h_2s->SetMaximum(1.e04);
    

    h_2s->Draw("3A");
    h_1s->Draw("3,SAME");
    h_exp->Draw("LP,SAME");
    
    h_obs->Draw("LP,SAME");
   // h_theo->Draw("LP,SAME");
        
    h_2s->GetXaxis()->SetTitle( "#it{m}_{med} [GeV]" );
    h_2s->GetYaxis()->SetTitle( "95% CL limit on #sigma_{obs}/#sigma_{theo}" );


    TLegend *leg = new TLegend(0.60,0.72,0.90,0.92,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

    leg->AddEntry( h_exp, "Expected", "PL" );
    leg->AddEntry( h_obs, "Observed", "PL" );

    leg->AddEntry( h_1s, "Expected #pm 1 s.d", "F" );
    leg->AddEntry( h_2s, "Expected #pm 2 s.d", "F" );

    addText(0.7-0.15,0.995-0.15,0.94,0.996,"35.9 fb^{-1} (13 TeV)",kBlack);    
    leg->Draw();
    
    
   TPaveText *pt = new TPaveText(0.32,0.746,0.45,0.798,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetLineColor(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.03);
   TText *AText = pt->AddText("#splitline{Scalar mediator, g_{q} = 1.0}{Dirac DM, #it{m}_{DM} = 1 GeV, g_{DM} = 1.0 }");
   AText->SetTextAlign(22);
   pt->Draw();
   
       
   pt = new TPaveText(0.10,0.846,0.48,0.898,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetLineColor(0);
   pt->SetTextFont(42);
   AText = pt->AddText("#bf{CMS}");
   AText->SetTextAlign(22);
   pt->Draw();
    
    
    canv->SaveAs("mu_13TeV_NLOS_Mchi_1.png");
    canv->SaveAs("mu_13TeV_NLOS_Mchi_1.pdf");
}
