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


void plotDM_EWK_K1K2_obs_cls(TString myfolder = "")
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
    // And change the colormap
    const Int_t     NRGBs         = 3;
    const Int_t     NCont         = 15;
    Double_t stops[NRGBs] = { 0.00, 0.5,1.0 };
    Double_t red[NRGBs]   = { 254,252,222 };
    Double_t green[NRGBs] = { 224,146,45 };
    Double_t blue[NRGBs]  = { 210,114,38 };
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
    t1->SetLogz(false);
    t1->SetRightMargin(0.2);


    TH2D* h2 = new TH2D("h2","",1000,0,1300,100,0.1,10);

    h_obs->SetHistogram(h2);
    h_obs->Draw("COLZ");

    h_obs->SetMaximum(6e2);
    h_obs->SetMinimum(2e2);

    h_obs->GetXaxis()->SetTitle("#it{m}_{#chi} [GeV]");
    h_obs->GetYaxis()->SetTitle("c_{1}     ");
    h_obs->GetYaxis()->SetTitleOffset(0.6);
    h_obs->GetZaxis()->SetTitle("95% CL observed limit on #Lambda [GeV]");
    h_obs->GetZaxis()->SetRangeUser(2e2,6e2);


    addText(0.7-0.15,0.995-0.15,0.95,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.,0.45,0.85,0.92,"#bf{CMS}",kWhite);
    canv->SaveAs("EWKDM_13TeV_k1k2_obs_cls.png");
    canv->SaveAs("EWKDM_13TeV_k1k2_obs_cls.pdf");

    delete canv;
}
