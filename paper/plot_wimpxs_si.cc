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
double getSigmaSI(double mmed, double mx)
{
    // https://arxiv.org/pdf/1603.04156v1.pdf
    // Equation 4.3
    double redMass = (mx*mN)/(mx+mN);
    return 6.9e-29 * pow(redMass,2)/pow(mmed,4); //cm^2
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
    string parName = "#it{m}_{#chi} [GeV]";
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

    TGraph *monoZ_13TeV_2p3fb_0 = new TGraph("interpolate_MV_observed_0p25.txt","%lg %lg");
    //monoZ_13TeV_2p3fb_0 = sortGraph(monoZ_13TeV_2p3fb_0);
    TGraph *monoZ_13TeV_2p3fb = new TGraph();
    int nd9=monoZ_13TeV_2p3fb_0->GetN(); //get ploted array dimention
    double m_med[10000];
    double m_dm[10000];
    for(Int_t i=0; i<nd9; i++) {
        monoZ_13TeV_2p3fb_0->GetPoint(i,m_med[i],m_dm[i]);
        double wimpxs = getSigmaSI(m_med[i],m_dm[i]);
        monoZ_13TeV_2p3fb->SetPoint(i,m_dm[i],wimpxs);
    }

    dumpGraphToFile( monoZ_13TeV_2p3fb, "monoZll_observed_90CL_SI.txt" );

    monoZ_13TeV_2p3fb->SetLineWidth(305);
    monoZ_13TeV_2p3fb->SetFillStyle(3005);
    monoZ_13TeV_2p3fb->SetLineColor(kRed+1);
    monoZ_13TeV_2p3fb->SetFillColor(kRed+1);
    monoZ_13TeV_2p3fb->SetMarkerColor(kRed+1);
    monoZ_13TeV_2p3fb->SetLineStyle(1);
    monoZ_13TeV_2p3fb->SetMarkerStyle(20);
    monoZ_13TeV_2p3fb->SetMarkerSize(0.8);


    TString path_to_dat =  "/user/albert/cernbox/share/dmlimits/";

    TGraph *lux2016 = new TGraph(path_to_dat + "lux2016.txt","%lg %lg");
    lux2016->SetLineWidth(2);
    lux2016->SetLineColor(kPink+7);
    lux2016->SetLineStyle(1);
    lux2016->SetMarkerSize(0);

    TGraph *PandaXII = new TGraph(path_to_dat + "pandax.txt","%lg %lg");
    PandaXII->SetLineWidth(2);
    PandaXII->SetLineColor(kGreen-3);
    PandaXII->SetLineStyle(1);
    PandaXII->SetMarkerSize(0);

    TGraph *CresstII = new TGraph(path_to_dat + "cresstii.txt","%lg %lg");
    CresstII->SetLineWidth(2);
    CresstII->SetLineColor(kOrange-3);
    CresstII->SetLineStyle(1);
    CresstII->SetMarkerSize(0);

    TGraph *supercdms = new TGraph(path_to_dat + "cdmslite2015.txt","%lg %lg");
    supercdms->SetLineWidth(2.);
    supercdms->SetLineColor(kBlue+2);
    supercdms->SetLineStyle(1);
    supercdms->SetMarkerSize(0);

    mg->Add(monoZ_13TeV_2p3fb,"L");
    mg->Add(lux2016,"C");
    mg->Add(supercdms,"C");
    mg->Add(PandaXII,"C");
    mg->Add(CresstII,"C");

    mg->Draw("A");
    mg->SetMinimum(1e-50);
    mg->SetMaximum(1e-35);
    mg->GetXaxis()->SetLimits(m_dm[0], 1e3);
    mg->GetXaxis()->SetTitle(parName.c_str());
    mg->GetXaxis()->SetTitleOffset(1.00);

    mg->GetYaxis()->SetTitle("DM-nucleon cross section [cm^{2}]");
    mg->GetYaxis()->SetTitleOffset(1.30);



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


    //~ addText(0.21+0.55,0.37+0.55,0.22,0.12,"90% CL",kBlack);
    addText(0.18,0.95,0.24,0.14,"Spin-independent          90% CL",kBlack);

    addText(0.20, 0.32, 0.94, 0.878, "#bf{CMS}", kBlack);

    //#################################################
    float posx1 = 0.6;
    float posx2 = 0.91;
    float posy1 = 0.45;
    float posy2 = 0.65;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);


    leg->AddEntry(supercdms, "CDMSLite 2015", "L");
    leg->AddEntry(lux2016, "LUX 2016", "L");
    leg->AddEntry(PandaXII, "PandaX-II", "L");
    leg->AddEntry(CresstII, "CRESST-II", "L");

    leg->Draw();

    posx1 = 0.17;
    posx2 = 0.91;
    posy1 = 0.2;
    posy2 = 0.37;

    TLegend *leg2 = new TLegend(posx1, posy1, posx2, posy2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetLineColor(0);
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);

    leg2->AddEntry(monoZ_13TeV_2p3fb, "#splitline{Vector mediator, Dirac fermion #chi}{#it{g}_{q} = 0.25, #it{g}_{#chi} = 1 (13 TeV, 2.3 fb^{-1})}", "L");

    leg2->Draw();

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
        canv->SaveAs( ("./"+plotName+"_DM_NucleonXS_SI.root").c_str() );
    }

    //return;
}




