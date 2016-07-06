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

// mass of proton: 938.27203 MeV/c^2
// mass of neutron: 939.56536 MeV/c^2
// reduced mass : 469.4591 MeV/c^2
const double mN=0.93827203;
double getSigmaSD(double obs_lam, double mx)
{
    // https://indico.cern.ch/event/489596/contributions/2009079/attachments/1227520/1797906/LHC_DM_WG_recommendation_document_V1.pdf
    // Equation 4.9
    double redMass = (mx*mN)/(mx+mN);
    double ret = 3.8e-29*pow(redMass,2)/pow(obs_lam,4);
    return ret; //cm^2
}



void plot_wimpxs_sd()
{
    string dbpars[] = {"1","10", "200", "500", "1000"};
    string parName = "WIMP mass #it{m}_{#chi} [GeV]";
    bool savePlots = true;
    bool showObserved = false;

    TCanvas *canv = new TCanvas("canv", "limits canvas", 600,600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);
    TStyle * TDR = createTdrStyle();
    TDR->cd();

    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    t1->SetLogx(true);

    
    TMultiGraph *mg = new TMultiGraph();

    TGraph *monoZ_13TeV_2p3fb_0 = new TGraph("interpolate_MA_observed.txt","%lg %lg");
    TGraph *monoZ_13TeV_2p3fb = new TGraph();
    int nd9=monoZ_13TeV_2p3fb_0->GetN(); //get ploted array dimention
    double m_med[3000];
    double m_dm[3000];
    for(Int_t i=0; i<nd9; i++) {
        monoZ_13TeV_2p3fb_0->GetPoint(i,m_med[i],m_dm[i]);
        double wimpxs = getSigmaSD(m_med[i],m_dm[i]);
        monoZ_13TeV_2p3fb->SetPoint(i,m_dm[i],wimpxs);
    }
    
    monoZ_13TeV_2p3fb->SetLineWidth(3);
    monoZ_13TeV_2p3fb->SetLineColor(kRed+1);
    monoZ_13TeV_2p3fb->SetMarkerColor(kRed+1);
    monoZ_13TeV_2p3fb->SetLineStyle(1);
    monoZ_13TeV_2p3fb->SetMarkerStyle(20);
    
    monoZ_13TeV_2p3fb->SetMarkerSize(0.8);


    TGraph *coupp = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/COUPP_2012_SD_flat_efficiency_model.dat","%lg %lg");
    coupp->SetLineWidth(2);
    coupp->SetLineColor(kPink+1);
    coupp->SetLineStyle(9);
    coupp->SetMarkerSize(0);


    TGraph *simple = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/SIMPLE_2011_SD.dat","%lg %lg");
    simple->SetLineWidth(2);
    simple->SetLineColor(kViolet-1);
    simple->SetLineStyle(3);
    simple->SetMarkerSize(0);


    TGraph *IceCube = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/IceCube_2011_SD_WW.dat","%lg %lg");
    IceCube->SetLineWidth(2);
    IceCube->SetLineColor(kViolet+3);
    IceCube->SetLineStyle(7);
    IceCube->SetMarkerSize(0);


    TGraph *IceCube_bb = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/IceCube_2011_SD_bb.dat","%lg %lg");
    IceCube_bb->SetLineWidth(2);
    IceCube_bb->SetLineColor(kViolet+3);
    IceCube_bb->SetLineStyle(3);
    IceCube_bb->SetMarkerSize(0);


    TGraph *picasso = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/PICASSO_2012_SD_Limits.dat","%lg %lg");
    picasso->SetLineWidth(2);
    picasso->SetLineColor(kCyan+1);
    picasso->SetLineStyle(6);
    picasso->SetMarkerSize(0);


    TGraph *pico2l = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/PICO-2L_SDp.dat","%lg %lg");
    pico2l->SetLineWidth(2);
    pico2l->SetLineColor(kPink+1);
    pico2l->SetLineStyle(9);
    pico2l->SetMarkerSize(0);


    TGraph *xenon100 = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SD/xenon100_90CL1301.6620v2.dat","%lg %lg");
    xenon100->SetLineWidth(2);
    xenon100->SetLineColor(kViolet-1);
    xenon100->SetLineStyle(3);
    xenon100->SetMarkerSize(0);


    mg->Add(pico2l);
    mg->Add(xenon100);
    mg->Add(IceCube);
 
    mg->Draw("ACP");
    monoZ_13TeV_2p3fb->Draw("P");
    mg->SetMinimum(1e-43);
    mg->SetMaximum(1e-35);
    mg->GetXaxis()->SetRangeUser(1,1100);
    mg->GetXaxis()->SetTitle(parName.c_str());

    mg->GetYaxis()->SetTitle("WIMP-nucleon cross section [cm^{2}]");

    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    char Buffer[1024];
    double iEcm_8  = 8.0;
    double iEcm_7  = 7.0;
    double iLumi_8 = 19712;
    double iLumi_7 = 5051;


    addText(0.21+0.55,0.37+0.55,0.22,0.12,"90% CL",kGray+2);
    addText(0.18,0.45,0.3,0.05,"Spin dependent",kGray+2);
    addText(0.70,0.83,0.63+0.05,0.56+0.05,"PICO-2L",kPink+1,18);
    addText(0.38-0.1,0.54-0.1,0.77,0.63,"XENON100",kViolet-1,-65);
    addText(0.75,0.95,0.53,0.46,"IceCube (WW)",kViolet+3,8);



    float posx1 = 0.2;
    float posx2 = 0.55;
    float posy1 = 0.25;
    float posy2 = 0.45;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);


    leg->AddEntry(monoZ_13TeV_2p3fb, "13 TeV, 2.3 fb^{-1}", "L");

    leg->Draw();

    //~ addText(0.76-0.6+0.02,0.96-0.6+0.02,0.835,0.898,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);
    addText(0.21,0.5,0.92,0.82,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);

    addText(0.77-0.05,0.97-0.05,0.755+0.1,0.818+0.1,"#splitline{#it{Axial vector}}{#it{coupling}}",kBlack);
    addText(0.67,0.90,0.345,0.408,"#it{Observed limits}",kBlack);

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
        plotName = "Zwimps_13TeV";
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SD.png").c_str() );
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SD.pdf").c_str() );
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SD.eps").c_str() );
    }

    //return;
}




