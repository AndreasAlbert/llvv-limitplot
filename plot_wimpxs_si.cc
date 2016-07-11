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

TGraph *superCDMS() {
    int i0 = -1;
    double *lX = new double[1000];
    double *lY = new double[1000];
    i0++; lX[i0] = 3.5946953342351033; lY[i0]= 9.848667679721256e-40;
    i0++; lX[i0] = 3.75095871431554;   lY[i0] = 4.7369269651270765e-40;
    i0++; lX[i0] = 4.055296786109914;  lY[i0] = 1.731452556791703e-40;
    i0++; lX[i0] = 4.859185478655973;  lY[i0] = 3.335637195211155e-41;
    i0++; lX[i0] = 5.443071002530319;  lY[i0] = 1.4979478319670555e-41;
    i0++; lX[i0] = 6.568479751568511;  lY[i0] = 4.257347230887382e-42;
    i0++; lX[i0] = 7.759712545805536;  lY[i0] = 1.4200965776277166e-42;
    i0++; lX[i0] = 10.195931085528889; lY[i0] = 5.073387251235342e-43;
    i0++; lX[i0] = 13.636666866311272; lY[i0] = 2.930133564061426e-43;
    i0++; lX[i0] = 18.04552501288121;  lY[i0] = 2.0321090962591862e-43;
    i0++; lX[i0] = 23.377104302070727; lY[i0] = 1.941241066769181e-43;
    i0++; lX[i0] = 29.857390063416926; lY[i0] = 2.0321090962591862e-43;
    TGraph *lLimit = new TGraph(i0,lX,lY);
    lLimit->SetLineWidth(2.);
    lLimit->SetLineColor(kBlue+2);
    return lLimit;
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

    TGraph *monoZ_13TeV_2p3fb_0 = new TGraph("interpolate_MV_observed.txt","%lg %lg");
    monoZ_13TeV_2p3fb_0 = sortGraph(monoZ_13TeV_2p3fb_0);
    TGraph *monoZ_13TeV_2p3fb = new TGraph();
    int nd9=monoZ_13TeV_2p3fb_0->GetN(); //get ploted array dimention
    double m_med[3000];
    double m_dm[3000];
    for(Int_t i=0; i<nd9; i++) {
        monoZ_13TeV_2p3fb_0->GetPoint(i,m_med[i],m_dm[i]);
        double wimpxs = getSigmaSI(m_med[i],m_dm[i]);
        monoZ_13TeV_2p3fb->SetPoint(i,m_dm[i],wimpxs);
    }

    monoZ_13TeV_2p3fb->SetLineWidth(3);
    monoZ_13TeV_2p3fb->SetLineColor(kRed+1);
    monoZ_13TeV_2p3fb->SetMarkerColor(kRed+1);
    monoZ_13TeV_2p3fb->SetLineStyle(1);
    monoZ_13TeV_2p3fb->SetMarkerStyle(20);
    monoZ_13TeV_2p3fb->SetMarkerSize(0.8);


        TString path_to_dat =  "/.automount/home/home__home1/institut_3a/albert/code/llvv/limitplots/dat/";

    TGraph *Lux2013 = new TGraph(path_to_dat + "LUX_2013_90CL.dat","%lg %lg");
    Lux2013->SetLineWidth(2);
    Lux2013->SetLineColor(kPink+7);
    //Lux2013->SetFillStyle();
    //Lux2013->SetFillColor(kPink+3);
    //Lux2013->SetLineWidth(3002);
    Lux2013->SetLineStyle(1);
    Lux2013->SetMarkerSize(0);


    TGraph * supercdms = superCDMS();

    mg->Add(monoZ_13TeV_2p3fb,"L");
    mg->Add(Lux2013,"C");
    mg->Add(supercdms,"C");

    mg->Draw("A");
    mg->SetMinimum(1e-47);
    mg->SetMaximum(1e-28);
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
    addText(0.21,0.5,0.92,0.82,"#splitline{#bf{CMS}}{#it{Preliminary}}",kBlack);
    addText(0.77-0.05,0.97-0.05,0.755+0.1,0.818+0.1,"#splitline{#it{Vector}}{#it{coupling}}",kBlack);
    addText(0.67,0.90,0.345,0.408,"#it{Observed limits}",kBlack);

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
    leg->AddEntry(supercdms, "SuperCDMS", "L");
    leg->AddEntry(Lux2013, "LUX", "L");

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




