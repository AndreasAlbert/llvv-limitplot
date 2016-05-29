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

#include "Utils.h"


using namespace std;
const double mN=0.93827203;
double getSigmaSI(double obs_lam, double mx)
{
    // https://indico.cern.ch/event/489596/contributions/2009079/attachments/1227520/1797906/LHC_DM_WG_recommendation_document_V1.pdf
    // Equation 4.3
    double redMass = (mx*mN)/(mx+mN);
    return 1.1e-27 * pow(redMass,2)/pow(obs_lam,4); //cm^2
}
void plot_wimpxs_si()
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
    
    TGraph *monoZ_13TeV_2p3fb_0 = new TGraph("interpolate_MV.txt","%lg %lg");
    TGraph *monoZ_13TeV_2p3fb = new TGraph();
    int nd9=monoZ_13TeV_2p3fb_0->GetN(); //get ploted array dimention
    double ax[3000];
    double ay[3000];
    for(Int_t i=0; i<nd9; i++) {
        monoZ_13TeV_2p3fb_0->GetPoint(i,ax[i],ay[i]);
        double x2 = getSigmaSI(ay[i],ax[i]);
        monoZ_13TeV_2p3fb->SetPoint(i,ax[i],x2);
    }
        
    monoZ_13TeV_2p3fb->SetLineWidth(3);
    monoZ_13TeV_2p3fb->SetLineColor(kRed+1);
    monoZ_13TeV_2p3fb->SetMarkerColor(kRed+1);
    monoZ_13TeV_2p3fb->SetLineStyle(1);
    monoZ_13TeV_2p3fb->SetMarkerStyle(20);
    monoZ_13TeV_2p3fb->SetMarkerSize(0.8);
    
    TGraph *Lux2013 = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SI/LUX_2013_90CL.dat","%lg %lg");
    Lux2013->SetLineWidth(2);
    Lux2013->SetLineColor(kPink+7);
    //Lux2013->SetFillStyle();
    //Lux2013->SetFillColor(kPink+3);
    //Lux2013->SetLineWidth(3002);
    Lux2013->SetLineStyle(9);
    Lux2013->SetMarkerSize(0);


    TGraph *cdmslite = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SI/CDMSlite.dat","%lg %lg");
    cdmslite->SetLineWidth(2);
    cdmslite->SetLineColor(kViolet-1);
    cdmslite->SetLineStyle(2);
    cdmslite->SetMarkerSize(0);


    TGraph *superCDMS2014 = new TGraph("~/code/llvv/limitplots/wimpxs/wimpPlotter_VectorM/wimpPlotter_SI/superCDMS_2014.dat","%lg %lg");
    superCDMS2014->SetLineWidth(2);
    superCDMS2014->SetLineColor(kViolet+2);
    superCDMS2014->SetLineStyle(8);
    superCDMS2014->SetMarkerSize(0);

    mg->Add(monoZ_13TeV_2p3fb,"P");
    mg->Add(Lux2013,"CP");
    mg->Add(cdmslite,"CP");

    mg->Draw("A");
    mg->SetMinimum(1e-46);
    mg->SetMaximum(1e-29);
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
    addText(0.18,0.45,0.3,0.05,"Spin Independent",kGray+2);
    addText(0.21,0.5,0.92,0.82,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);
    addText(0.77-0.05,0.97-0.05,0.755+0.1,0.818+0.1,"#splitline{#it{Vector}}{#it{coupling}}",kBlack);
    addText(0.67,0.90,0.345,0.408,"#it{Expected limits}",kBlack);

    addText(0.75,0.92,0.29,0.24,"LUX 2013",kPink+7,10);
    addText(0.41,0.55,0.35,0.5,"CDMSlite",kViolet-1);
    //#################################################
    float posx1 = 0.2;
    float posx2 = 0.55;
    float posy1 = 0.6;
    float posy2 = 0.7;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);


    leg->AddEntry(monoZ_13TeV_2p3fb, "13 TeV, 2.3 fb^{-1}", "L");

    leg->Draw();


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
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SI.png").c_str() );
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SI.pdf").c_str() );
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SI.eps").c_str() );
    }

    //return;
}




