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
#include "TGraph2D.h"
#include "TGraphSmooth.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"

#include "Utils.h"

using namespace std;


void plot2D_DMA_obs_1p00_cls(TString myfolder = "", TString tag="", bool savePlots=true)
{


    // if(myfolder.Contains("/MV")) tag="MV";
    // if(myfolder.Contains("/MA")) tag="MA";

    tag="MA";

    bool showObserved = true;
    double sigma_theo = 0.15;
    vector<int> dm_masses = {1, 2, 5, 7, 10, 12, 17, 22, 25, 27, 32, 37, 42, 47, 50, 52, 57, 60, 62, 67, 70, 72, 77, 80, 82, 87, 90, 92, 97, 100, 102, 107, 110, 112, 117, 122, 125, 127, 132, 135, 137, 142, 147, 150, 152, 157, 162, 167, 172, 175, 177, 182, 187, 192, 197, 200, 202, 207, 212, 217, 222, 227, 232, 237, 242, 247, 250, 252, 257, 262, 267, 272, 277, 282, 287, 292, 297, 300, 400, 500, 600, 700, 800, 900, 1000};


    /// Map:    DM Mass -> Mediator masses
    std::map<int, std::vector<int>> MassMap;

    MassMap[    1 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[    2 ] = {14};
    MassMap[    5 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[    7 ] = {5};
    MassMap[   10 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   12 ] = {34};
    MassMap[   17 ] = {25};
    MassMap[   22 ] = {54};
    MassMap[   25 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   27 ] = {45};
    MassMap[   32 ] = {74};
    MassMap[   37 ] = {65};
    MassMap[   42 ] = {94};
    MassMap[   47 ] = {85};
    MassMap[   50 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   52 ] = {114};
    MassMap[   57 ] = {105};
    MassMap[   60 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   62 ] = {134};
    MassMap[   67 ] = {125};
    MassMap[   70 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   72 ] = {154};
    MassMap[   77 ] = {145};
    MassMap[   80 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   82 ] = {174};
    MassMap[   87 ] = {165};
    MassMap[   90 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[   92 ] = {194};
    MassMap[   97 ] = {185};
    MassMap[  100 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  102 ] = {214};
    MassMap[  107 ] = {205};
    MassMap[  110 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  112 ] = {234};
    MassMap[  117 ] = {225};
    MassMap[  122 ] = {254};
    MassMap[  125 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  127 ] = {245};
    MassMap[  132 ] = {274};
    MassMap[  135 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  137 ] = {265};
    MassMap[  142 ] = {294};
    MassMap[  147 ] = {285};
    MassMap[  150 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  152 ] = {314};
    MassMap[  157 ] = {305};
    MassMap[  162 ] = {334};
    MassMap[  167 ] = {325};
    MassMap[  172 ] = {354};
    MassMap[  175 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  177 ] = {345};
    MassMap[  182 ] = {374};
    MassMap[  187 ] = {365};
    MassMap[  192 ] = {394};
    MassMap[  197 ] = {385};
    MassMap[  200 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  202 ] = {414};
    MassMap[  207 ] = {405};
    MassMap[  212 ] = {434};
    MassMap[  217 ] = {425};
    MassMap[  222 ] = {454};
    MassMap[  227 ] = {445};
    MassMap[  232 ] = {474};
    MassMap[  237 ] = {465};
    MassMap[  242 ] = {494};
    MassMap[  247 ] = {485};
    MassMap[  250 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  252 ] = {514};
    MassMap[  257 ] = {505};
    MassMap[  262 ] = {534};
    MassMap[  267 ] = {525};
    MassMap[  272 ] = {554};
    MassMap[  277 ] = {545};
    MassMap[  282 ] = {574};
    MassMap[  287 ] = {565};
    MassMap[  292 ] = {594};
    MassMap[  297 ] = {585};
    MassMap[  300 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  400 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  500 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  600 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  700 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  800 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[  900 ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 1000 ] = {    20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};



    TGraph2D* h_Limit_obs      = new TGraph2D(); h_Limit_obs->SetName("Observed");
    TGraph2D* h_Limit_exp      = new TGraph2D(); h_Limit_exp->SetName("Expected");
    TGraph2D* h_Limit_p1s      = new TGraph2D(); h_Limit_p1s->SetName("ExpectedPlusSigma");
    TGraph2D* h_Limit_m1s      = new TGraph2D(); h_Limit_m1s->SetName("ExpectedMinusSigma");

    TGraph2D* h_Limit_theo_m1s = new TGraph2D(); h_Limit_theo_m1s->SetName("ObservedMinusSigma");
    TGraph2D* h_Limit_theo_p1s = new TGraph2D(); h_Limit_theo_p1s->SetName("ObservedPlusSigma");

    int npoint=0;
    for(size_t nmx=0; nmx<dm_masses.size(); nmx++) {

      if(dm_masses[nmx]>300.) continue;

        std::vector<int> mv_masses = MassMap[ dm_masses[nmx] ];

        int sizemv = mv_masses.size();

        for(int nmv=0; nmv<sizemv; nmv++) {

	  if(mv_masses[nmv]>1150.) continue;
            TString str_mx = Convert2TString(dm_masses[nmx]);
            TString str_mv = Convert2TString(mv_masses[nmv]);
            TString folder = myfolder+"/"+tag+"_"+str_mx+"_"+str_mv+"_1p00";
            Double_t obs = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.root" );
            Double_t m1s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.160.root" );
            Double_t exp = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.500.root" );
            Double_t p1s = get_limit_from_file( folder + "/higgsCombineZwimps01jets.HybridNew.mH120.quant0.840.root" );

            if( obs < 0 ) continue;
            if(obs<0.1) obs = 0.1;

            h_Limit_obs->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],obs);
            h_Limit_exp->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],exp);
            h_Limit_p1s->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],p1s);
            h_Limit_m1s->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],m1s);

            // For the theory uncertainty, vary the limit on mu by the uncertainty
            h_Limit_theo_p1s->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],(1+sigma_theo)*obs);
            h_Limit_theo_m1s->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],(1-sigma_theo)*obs);
            npoint++;
        }


    }
    std::cout << "Reading Finished." << std::endl;

    // Format:
    //   InterpolateDM(TString tag, TGraph2D* h_Limit, double sqrt_gxgq=1.0, double xmin=1, double xmax=1000, double ymin=1, double ymax=200, double xctr=100, double yctr=1, int nsteps=100)
    TGraph* exclusion_exp          = InterpolateDM(tag, h_Limit_exp     , 1,  1, 800, 0, 250, 150, 1, 100);
    TGraph* exclusion_p1s          = InterpolateDM(tag, h_Limit_p1s     , 1,  1, 800, 0, 250, 150, 1, 100);
    TGraph* exclusion_m1s          = InterpolateDM(tag, h_Limit_m1s     , 1,  1, 800, 0, 250, 150, 1, 100);
    TGraph* exclusion_obs          = InterpolateDM(tag, h_Limit_obs     , 1,  1, 800, 0, 250, 150, 1, 100);
    TGraph* exclusion_obs_theo_p1s = InterpolateDM(tag, h_Limit_theo_p1s, 1,  1, 800, 0, 250, 150, 1, 100);
    TGraph* exclusion_obs_theo_m1s = InterpolateDM(tag, h_Limit_theo_m1s, 1,  1, 800, 0, 250, 150, 1, 100);

    TGraph* diagonal = new TGraph();
    diagonal->SetPoint(0,0,0);
    diagonal->SetPoint(1,2e3,1e3);
    diagonal->SetLineWidth(2);
    diagonal->SetLineStyle(9);
    diagonal->SetLineColor( getDiagonalColor() );

    dumpGraphToFile( exclusion_exp, "interpolate_"+tag+"_expected_1p00.txt" );
    dumpGraphToFile( exclusion_obs, "interpolate_"+tag+"_observed_1p00.txt" );

    exclusion_exp -> SetMarkerColor( kBlack );
    exclusion_p1s -> SetMarkerColor( kGray );
    exclusion_m1s -> SetMarkerColor( kGray );
    exclusion_obs -> SetMarkerColor( getObservedColor() );
    exclusion_exp -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineStyle( 3 );
    exclusion_m1s -> SetLineColor( kBlack );
    exclusion_m1s -> SetLineStyle( 3 );
    exclusion_obs -> SetLineColor( getObservedColor() );

    exclusion_obs_theo_p1s -> SetLineStyle( 7 );
    exclusion_obs_theo_m1s -> SetLineStyle( 7 );
    exclusion_obs_theo_p1s -> SetLineColor(getObservedColor());
    exclusion_obs_theo_m1s -> SetLineColor(getObservedColor());

    exclusion_exp -> SetLineWidth( 2 );
    exclusion_p1s -> SetLineWidth( 2 );
    exclusion_m1s -> SetLineWidth( 2 );
    exclusion_obs -> SetLineWidth( 2 );
    exclusion_obs_theo_p1s -> SetLineWidth( 2 );
    exclusion_obs_theo_m1s -> SetLineWidth( 2 );
    // Use TDR as basis
    TStyle * TDR = createTdrStyle();
    TDR->cd();
    // And change the colormap
    setColorMap();

    // Relic Density
    TGraph * h_relic = new TGraph("relic_a_1.txt");
    h_relic->SetLineColor(getRelicColor());
    h_relic->SetLineStyle(1);
    h_relic->SetLineWidth(2);
    h_relic->SetFillStyle(3005);
    h_relic->SetFillColor(getRelicColor());
    TGraph * h_relic2 = new TGraph("list_scan_A_g25_MD_xxd_A_gq1.0_gdm1.0_contour_10.txt");
    h_relic2->SetLineColor(getRelicColor());
    h_relic2->SetLineStyle(1);
    h_relic2->SetLineWidth(2);
    h_relic2->SetFillStyle(3005);
    h_relic2->SetFillColor(getRelicColor());

    ///////////////////////////////////////////
    //////////// TRANSPOSED PLOT
    std::cout << "Create Canvas and pad." << std::endl;
    TCanvas *canv = new TCanvas("c2", "c2", 700, 600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);

    TPad *t1 = new TPad("t2","t2", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogz(true);
    //~ t1->SetLogy(true);
    //~ t1->SetLogx(true);
    t1->SetRightMargin(0.2);


    TH2D* h3 = new TH2D("h3","",100,1,1e3,100,1,2.5e2);

    h_Limit_obs->GetYaxis()->SetLimits(1,5e+2);
    h_Limit_obs->SetHistogram(h3);

    h_Limit_obs->Draw("COLZ");
    diagonal->Draw("L SAMES");
    //~ h_relic->Draw("L SAMES");
    h_relic2->Draw("L SAMES");
    exclusion_exp->Draw("L SAMES");
    exclusion_p1s->Draw("L SAMES");
    exclusion_m1s->Draw("L SAMES");
    exclusion_obs->Draw("L SAMES");
    exclusion_obs_theo_p1s->Draw("L SAMES");
    exclusion_obs_theo_m1s->Draw("L SAMES");

    h_Limit_obs->SetMaximum(50);
    h_Limit_obs->SetMinimum(0.5);

    h_Limit_obs->GetYaxis()->SetTitle("#it{m}_{#chi} [GeV]");
    h_Limit_obs->GetXaxis()->SetTitle("#it{M}_{med} [GeV]");
    h_Limit_obs->GetZaxis()->SetTitle("95% CL observed limit on #sigma_{obs}/#sigma_{theo}");
    h_Limit_obs->GetZaxis()->SetTitleOffset(1.15);
    h_Limit_obs->GetZaxis()->SetRangeUser(0.1,10);

    addText(0.55, 0.84, 0.946, 0.996, "2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.16, 0.27, 0.88, 0.93, "#bf{CMS}", kBlack);

    double dy = 0.03;
    double dx = -0.085;
    if( tag == "MV" ) addText(0.50+dx,0.75+dx,0.60+dy,0.68+dy,"Vector mediator",kBlack, 0, 0.04);
    if( tag == "MA" ) addText(0.50+dx,0.75+dx,0.60+dy,0.68+dy,"Axial-vector mediator",kBlack, 0, 0.04);
    addText(0.50+dx,0.75+dx,0.52+dy,0.60+dy,"Dirac fermion #chi",kBlack, 0, 0.04);
    addText(0.50+dx,0.75+dx,0.44+dy,0.52+dy,"#it{g}_{#chi} = 1, #it{g}_{q} = 1",kBlack, 0, 0.04);
    addText(-0.12,0.9,0.6+dy,1.0+dy,"#it{M}_{med} = 2#it{m}_{#chi}",getDiagonalColor(), 65, 0.04);
    addText(0.42,0.65,0.35,0.45,"#Omega_{c} #times h^{2} = 0.12",getRelicColor(), 45, 0.04);

    float posx1 = 0.48;
    float posx2 = 0.78;
    float posy1 = 0.72;
    float posy2 = 0.92;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetLineColor(kBlack);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->SetHeader("#sigma/#sigma_{theo}=1:");
    leg->AddEntry(exclusion_obs, "Observed", "L");
    leg->AddEntry(exclusion_obs_theo_p1s, "Theo. uncertainty", "L");
    leg->AddEntry(exclusion_exp, "Expected", "L");
    leg->AddEntry(exclusion_p1s, "Expected #pm 1 s.d.", "L");


    leg->Draw();
    //~ leg->Draw();
    if(savePlots){
        canv->SaveAs("Limit_Mx_"+tag+"_obs_1p00_cls.png");
        canv->SaveAs("Limit_Mx_"+tag+"_obs_1p00_cls.pdf");
    }

    delete canv;

}
