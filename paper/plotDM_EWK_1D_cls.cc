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



void plotDM_EWK_1D_cls(TString myfolder = "")
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

    std::map<int,std::map<double,double>> xsec = {};
    xsec[800] = {{0.5,1.090941e-08},{1.0,1.1545e-08},{2.0,1.52245e-08},{3.0,2.21271e-08},{0.1,1.09813e-08},{0.2,1.09132e-08},{10.0,1.60578e-07},{0.3,1.0879e-08},{5.0,4.5566e-08}};
    xsec[1] = {{0.5,1.0535389183e-07},{1.0,1.121410811e-07},{2.0,1.476546082e-07},{3.0,2.117532859e-07},{0.1,1.0507171551e-07},{0.2,1.0466175528e-07},{10.0,1.47342024e-06},{0.3,1.0465149672e-07},{5.0,4.274439559e-07}};
    xsec[100] = {{0.5,9.2973e-08},{1.0,9.9009e-08},{2.0,1.30331e-07},{3.0,1.87483e-07},{0.1,9.2816e-08},{0.2,9.2463e-08},{10.0,1.30928e-06},{0.3,9.23959e-08},{5.0,3.7907e-07}};
    xsec[200] = {{0.5,7.34151e-08},{1.0,7.8078e-08},{2.0,1.027757e-07},{3.0,1.48132e-07},{0.1,7.3411e-08},{0.2,7.30941e-08},{10.0,1.0429e-06},{0.3,7.30237e-08},{5.0,3.0083e-07}};
    xsec[10] = {{0.5,1.0504691631e-07},{1.0,1.1196014804e-07},{2.0,1.4736422406e-07},{3.0,2.1126531429e-07},{0.1,1.04743093357e-07},{0.2,1.043530992e-07},{10.0,1.4665013603e-06},{0.3,1.0421880487e-07},{5.0,4.259205395e-07}};
    xsec[400] = {{0.5,4.07876e-08},{1.0,4.3298e-08},{2.0,5.70324e-08},{3.0,8.2485e-08},{0.1,4.09251e-08},{0.2,4.07138e-08},{10.0,5.8865e-07},{0.3,4.06172e-08},{5.0,1.68632e-07}};
    xsec[50] = {{0.5,1.012787e-07},{1.0,1.0785e-07},{2.0,1.41856e-07},{3.0,2.03561e-07},{0.1,1.01077e-07},{0.2,1.00691e-07},{10.0,1.41727e-06},{0.3,1.006385e-07},{5.0,4.1107e-07}};
    xsec[1300] = {{0.5,2.007769e-09},{1.0,2.1273e-09},{2.0,2.80375e-09},{3.0,4.06586e-09},{0.1,2.01819e-09},{0.2,2.00655e-09},{10.0,2.9253e-08},{0.3,2.00098e-09},{5.0,8.337e-09}};


    TGraph* h_exp = new TGraph();
    TGraph* h_obs = new TGraph();
    TGraphAsymmErrors* h_1s = new TGraphAsymmErrors();
    TGraphAsymmErrors* h_2s = new TGraphAsymmErrors();

    int npoint=0;

    for(int nmx=0; nmx<dm_masses.size(); nmx++) {

        TString str_mx = Convert2TString(dm_masses[nmx]);
        str_mx += "_";
        str_mx += "1";

        TString folder = myfolder+"/MEWK_"+str_mx;
        double scale = 1/( xsec[dm_masses[nmx]][1] / 5e-9);
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

        std::cout << "------- " << str_mx << std::endl;
        std::cout << Lambda_m2s << std::endl;
        std::cout << Lambda_m1s << std::endl;
        std::cout << Lambda_exp << " " <<Lambda_obs << std::endl;
        std::cout << Lambda_p1s << std::endl;
        std::cout << Lambda_p2s << std::endl;

        double err_up1   = fabs(std::max(Lambda_p1s,Lambda_m1s) - Lambda_exp);
        double err_down1 = fabs(std::min(Lambda_p1s,Lambda_m1s) - Lambda_exp);
        double err_up2   = fabs(std::max(Lambda_p2s,Lambda_m2s) - Lambda_exp);
        double err_down2 = fabs(std::min(Lambda_p2s,Lambda_m2s) - Lambda_exp);


        h_exp ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);
        h_1s    ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);
        h_2s    ->SetPoint(npoint,dm_masses[nmx],Lambda_exp);

        h_1s    ->SetPointError(npoint,0,0,err_down1,err_up1);
        h_2s    ->SetPointError(npoint,0,0,err_down2,err_up2);

        h_obs->SetPoint(npoint, dm_masses[nmx],Lambda_obs );
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

    h_obs->SetLineStyle(1);
    h_obs->SetMarkerStyle(20);

    h_1s->SetFillColor(kGreen+1);
    h_2s->SetFillColor(kOrange);

    h_2s->SetMinimum(150);
    h_2s->SetMaximum(900);


    h_2s->Draw("3A");
    h_1s->Draw("3,SAME");
    h_exp->Draw("LP,SAME");
    h_obs->Draw("LP,SAME");

    h_2s->GetXaxis()->SetTitle( "#it{m}_{#chi} [GeV]" );

    h_2s->GetXaxis()->SetRangeUser(-200,1350);
    h_2s->GetYaxis()->SetTitle( "95% CL limit on #Lambda [GeV]" );


    TLegend *leg = new TLegend(0.2,0.63,0.65,0.88);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);

    leg->AddEntry( h_obs, "Observed", "PL" );
    leg->AddEntry( h_exp, "Expected", "PL" );
    leg->AddEntry( h_1s, "Expected #pm 1 s.d.", "F" );
    leg->AddEntry( h_2s, "Expected #pm 2 s.d.", "F" );

    addText(0.63,0.995,0.95,1.0,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.08,0.45,0.88,.93,"#bf{CMS}",kBlack);

    addText(0.25,0.35,0.2,0.3,"c_{1} = 1",kBlack);
    leg->Draw();
    canv->SaveAs("EWKDM_13TeV_Lambda_cls.png");
    canv->SaveAs("EWKDM_13TeV_Lambda_cls.pdf");

    delete canv;
}
