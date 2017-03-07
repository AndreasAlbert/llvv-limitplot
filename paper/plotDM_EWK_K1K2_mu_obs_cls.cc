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

void plotDM_EWK_K1K2_mu_obs_cls(TString myfolder = "")
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

    std::map<int,std::map<double,double>> xsec = {};
    xsec[800] = {{0.5,1.090941e-08},{1.0,1.1545e-08},{2.0,1.52245e-08},{3.0,2.21271e-08},{0.1,1.09813e-08},{0.2,1.09132e-08},{10.0,1.60578e-07},{0.3,1.0879e-08},{5.0,4.5566e-08}};
    xsec[1] = {{0.5,1.0535389183e-07},{1.0,1.121410811e-07},{2.0,1.476546082e-07},{3.0,2.117532859e-07},{0.1,1.0507171551e-07},{0.2,1.0466175528e-07},{10.0,1.47342024e-06},{0.3,1.0465149672e-07},{5.0,4.274439559e-07}};
    xsec[100] = {{0.5,9.2973e-08},{1.0,9.9009e-08},{2.0,1.30331e-07},{3.0,1.87483e-07},{0.1,9.2816e-08},{0.2,9.2463e-08},{10.0,1.30928e-06},{0.3,9.23959e-08},{5.0,3.7907e-07}};
    xsec[200] = {{0.5,7.34151e-08},{1.0,7.8078e-08},{2.0,1.027757e-07},{3.0,1.48132e-07},{0.1,7.3411e-08},{0.2,7.30941e-08},{10.0,1.0429e-06},{0.3,7.30237e-08},{5.0,3.0083e-07}};
    xsec[10] = {{0.5,1.0504691631e-07},{1.0,1.1196014804e-07},{2.0,1.4736422406e-07},{3.0,2.1126531429e-07},{0.1,1.04743093357e-07},{0.2,1.043530992e-07},{10.0,1.4665013603e-06},{0.3,1.0421880487e-07},{5.0,4.259205395e-07}};
    xsec[400] = {{0.5,4.07876e-08},{1.0,4.3298e-08},{2.0,5.70324e-08},{3.0,8.2485e-08},{0.1,4.09251e-08},{0.2,4.07138e-08},{10.0,5.8865e-07},{0.3,4.06172e-08},{5.0,1.68632e-07}};
    xsec[50] = {{0.5,1.012787e-07},{1.0,1.0785e-07},{2.0,1.41856e-07},{3.0,2.03561e-07},{0.1,1.01077e-07},{0.2,1.00691e-07},{10.0,1.41727e-06},{0.3,1.006385e-07},{5.0,4.1107e-07}};
    xsec[1300] = {{0.5,2.007769e-09},{1.0,2.1273e-09},{2.0,2.80375e-09},{3.0,4.06586e-09},{0.1,2.01819e-09},{0.2,2.00655e-09},{10.0,2.9253e-08},{0.3,2.00098e-09},{5.0,8.337e-09}};


    TGraph2D* h_obs = new TGraph2D();
    TGraph2D* h_exp = new TGraph2D();
    TGraph2D* h_p1s = new TGraph2D();
    TGraph2D* h_m1s = new TGraph2D();
    TGraph2D* h_p2s = new TGraph2D();
    TGraph2D* h_m2s = new TGraph2D();

    int npoint=0;

    for(int nk1=0; nk1<K1.size(); nk1++) {

        for(int nmx=0; nmx<dm_masses.size(); nmx++) {

            TString str_mx = Convert2TString(dm_masses[nmx]);
            str_mx += "_";
            str_mx += K1[nk1];

            TString folder = myfolder+"/MEWK_"+str_mx;
            double scale = 1/(xsec[dm_masses[nmx]][K1_num[nk1]] / 5e-9);
            Double_t obs = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.root" );
            Double_t m2s = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.025.root" );
            Double_t m1s = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.160.root" );
            Double_t exp = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.500.root" );
            Double_t p1s = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.840.root" );
            Double_t p2s = scale * get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.975.root" );


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
            h_p2s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_p2s);
            h_m1s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_m1s);
            h_m2s->SetPoint(npoint,dm_masses[nmx],K1_num[nk1],Lambda_m2s);

            std::cout << npoint<< " " <<dm_masses[nmx]<<" " <<K1_num[nk1]<<" " <<exp << std::endl;
            npoint++;
        }

    }

    TGraph*   exclusion_exp =  InterpolateDMXY(h_exp,1,900,1299,0.11,7);
    TGraph*   exclusion_p1s =  InterpolateDMXY(h_p1s,1,800,1299,0.11,7);
    TGraph*   exclusion_p2s =  InterpolateDMXY(h_p2s,1,800,1299,0.11,7);
    TGraph*   exclusion_m1s =  InterpolateDMXY(h_m1s,1,800,1299,0.11,7);
    TGraph*   exclusion_m2s =  InterpolateDMXY(h_m2s,1,800,1299,0.11,7);
    TGraph*   exclusion_obs =  InterpolateDMXY(h_obs,1,800,1299,0.11,7);


    TGraph * shade = new TGraph();
    int n_p1s = exclusion_p1s->GetN();
    for( int n=0; n < n_p1s; n++ ) {
        double x(0),y(0);
        exclusion_p1s->GetPoint(n,x,y);
        if(x==0 or y==0) continue;
        shade->SetPoint(n,x,y);
        std::cout << n << " " << x << " " << y<< std::endl;
    }
    int n_m1s = exclusion_m1s->GetN();
    for( int n=0; n < n_m1s; n++  ) {
        double x(0),y(0);
        exclusion_m1s->GetPoint(n_m1s-n-1,x,y);
        if(x==0 or y==0) continue;
        shade->SetPoint(n_p1s+n,x,y);
        std::cout << n_p1s+n+1 << " " << x << " " << y << std::endl;
    }
    std::cout << shade->GetN() << " " << exclusion_p1s->GetN() << " " << exclusion_m1s->GetN() << std::endl;
    shade->SetFillStyle(3244);
    shade->SetFillColor(kGray);
    exclusion_p1s->SetFillStyle(3244);
    exclusion_p1s->SetFillColor(kGray);

    exclusion_exp -> SetMarkerColor( kBlack );
    exclusion_p1s -> SetMarkerColor( kGray );
    exclusion_m1s -> SetMarkerColor( kGray );
    exclusion_p2s -> SetMarkerColor( kGray );
    exclusion_m2s -> SetMarkerColor( kGray );
    exclusion_obs -> SetMarkerColor( kBlue );
    exclusion_exp -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineColor( kGray+2 );
    exclusion_m1s -> SetLineColor( kGray+2 );
    exclusion_p2s -> SetLineColor( kGray+2 );
    exclusion_m2s -> SetLineColor( kGray+2 );
    exclusion_obs -> SetLineColor( kBlue );
    exclusion_exp -> SetLineWidth(4);
    exclusion_p1s -> SetLineWidth(4);
    exclusion_m1s -> SetLineWidth(4);
    exclusion_obs -> SetLineWidth(4);
    exclusion_p1s -> SetLineStyle( 2 );
    exclusion_m1s -> SetLineStyle( 2 );
    exclusion_p2s -> SetLineStyle( 2 );
    exclusion_m2s -> SetLineStyle( 2 );
    exclusion_exp -> SetLineStyle( 7 );
    exclusion_obs -> SetLineStyle( 1 );

    exclusion_m1s->RemovePoint(1);

    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->SetPadLeftMargin(0.08);
    TDR->SetPadRightMargin(0.01);
    TDR->cd();
    gStyle->SetPalette(55);
    std::cout << "Create Canvas and pad." << std::endl;
    TCanvas *canv = new TCanvas("c", "c", 700, 550);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);
    canv->SetGridy(0);
    // And change the colormap
    const Int_t     NRGBs         = 3;
    const Int_t     NCont         = 15;
    Double_t stops[NRGBs] = { 0.00, 0.5,1.0 };
    Double_t red[NRGBs]   = { 222 ,252,    254   };
    Double_t green[NRGBs] = { 45  ,146,    224  };
    Double_t blue[NRGBs]  = { 38  ,114,    210  };
    for(int i = 0; i<NRGBs; i++) {
        red[i] = red[i]/255;
        green[i] = green[i]/255;
        blue[i] = blue[i]/255;
    }
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);

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
    shade->Draw("F SAME");
    exclusion_exp->Draw("L SAMES");
    exclusion_p1s->Draw("L SAMES");
    exclusion_m1s->Draw("L SAMES");
    exclusion_obs->Draw("L SAMES");
    h_obs->SetMaximum(5);
    h_obs->SetMinimum(5e-3);

    h_obs->GetXaxis()->SetTitle("#it{m}_{#chi} [GeV]");
    h_obs->GetYaxis()->SetTitle("c_{1}    ");
    h_obs->GetYaxis()->SetTitleOffset(0.6);
    h_obs->GetZaxis()->SetTitle("95% CL observed limit on #sigma_{obs}/#sigma_{theo}");
    h_obs->GetZaxis()->SetRangeUser(5e-3,5);


    addText(0.7-0.15,0.995-0.15,0.95,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.,0.45,0.85,0.92,"#bf{CMS}",kWhite);

    addText(0.4,0.9,0.84,0.92,"#Lambda = 300 GeV",kBlack);



    float posx1 = 0.12;
    float posx2 = 0.52;
    float posy1 = 0.15;
    float posy2 = 0.4;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetLineColor(kBlack);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->SetHeader("#sigma_{obs}/#sigma_{theo} = 1:");
    leg->AddEntry(exclusion_obs, "Observed", "L");
    leg->AddEntry(exclusion_exp, "Expected", "L");
    leg->AddEntry(exclusion_p1s, "Expected #pm 1 s.d.   ", "FL");

    leg->Draw();

    canv->SaveAs("EWKDM_13TeV_k1k2_mu_obs_cls.png");
    canv->SaveAs("EWKDM_13TeV_k1k2_mu_obs_cls.pdf");

    delete canv;
}
