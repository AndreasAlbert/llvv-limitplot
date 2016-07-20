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


void plot2D_DMV_obs_1p00(TString myfolder = "", TString tag="")
{


    if(myfolder.Contains("/MV")) tag="MV";
    if(myfolder.Contains("/MA")) tag="MA";

    bool showObserved = true;
    double sigma_theo = 0.05;
    //~ vector<int> dm_masses = {1,5,10,20,25,50,40,50,60,70,80,90,100,150,500,1000};
    vector<int> dm_masses = {1, 5,10,25, 50, 60, 70, 80, 90, 100, 110, 125, 135, 150, 175, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000};


    /// Map:    DM Mass -> Mediator masses
    std::map<int, std::vector<int>> MassMap;
    MassMap[ 1    ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 5    ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 10   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 25   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 50   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 60   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 70   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 80   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 90   ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 100  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 110  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 125  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 135  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 150  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 175  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 200  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 250  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 300  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 400  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 500  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 600  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 700  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 800  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 900  ] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};
    MassMap[ 1000 ] = {    20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000};

    TGraph2D* h_Limit_obs = new TGraph2D(); // Observed
    TGraph2D* h_Limit_exp = new TGraph2D(); // Expected
    TGraph2D* h_Limit_p1s = new TGraph2D(); // Expected - 1 sigma
    TGraph2D* h_Limit_m1s = new TGraph2D(); // Expected + 1 sigma

    TGraph2D* h_Limit_theo_m1s = new TGraph2D(); // Observed + 1 sigma_theo
    TGraph2D* h_Limit_theo_p1s = new TGraph2D(); // Observed - 1 sigma_theo

    int npoint=0;
    for(size_t nmx=0; nmx<dm_masses.size(); nmx++) {

        std::vector<int> mv_masses = MassMap[ dm_masses[nmx] ];

        int sizemv = mv_masses.size();

        for(int nmv=0; nmv<sizemv; nmv++) {

            TString str_mx = Convert2TString(dm_masses[nmx]);
            TString str_mv = Convert2TString(mv_masses[nmv]);
            std::cout << myfolder+"/"+tag+"_"+str_mx+"_"+str_mv+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" << std::endl;
            TFile myfile( myfolder+"/"+tag+"_"+str_mx+"_"+str_mv+"_gQ1p00/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

            Double_t obs(0.);
            Double_t m2s(0.);
            Double_t m1s(0.);
            Double_t exp(0.);
            Double_t p1s(0.);
            Double_t p2s(0.);

            if ( myfile.IsOpen() && !myfile.IsZombie() ) {

                Double_t lim, limerr;
                TTree *t = (TTree*)myfile.Get("limit");
                if( t == 0 ) continue;
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


            double sqrt_gxgq = exp;
            cout << "mx: " << dm_masses[nmx] << " mv: " << mv_masses[nmv] << " exp: " << sqrt_gxgq  << " p1s: " << p1s << " m1s: " << m1s << " obs: " << obs <<  endl;


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
    TGraph* exclusion_exp          = InterpolateDM(tag, h_Limit_exp     , 1, 10, 700, 0, 200);
    TGraph* exclusion_p1s          = InterpolateDM(tag, h_Limit_p1s     , 1, 10, 700, 0, 200);
    TGraph* exclusion_m1s          = InterpolateDM(tag, h_Limit_m1s     , 1, 10, 700, 0, 200);
    TGraph* exclusion_obs          = InterpolateDM(tag, h_Limit_obs     , 1, 10, 700, 0, 200);
    TGraph* exclusion_obs_theo_p1s = InterpolateDM(tag, h_Limit_theo_p1s, 1,  1, 700, 0, 200);
    TGraph* exclusion_obs_theo_m1s = InterpolateDM(tag, h_Limit_theo_m1s, 1,  1, 700, 0, 200);

    TGraph* diagonal = new TGraph();
    diagonal->SetPoint(0,0,0);
    diagonal->SetPoint(1,2e3,1e3);
    diagonal->SetLineWidth(2);
    diagonal->SetLineColor( kGray );
    
    dumpGraphToFile( exclusion_exp, "interpolate_"+tag+"_expected_1p00.txt" );
    dumpGraphToFile( exclusion_obs, "interpolate_"+tag+"_observed_1p00.txt" );

    exclusion_exp -> SetMarkerColor( kBlack );
    exclusion_p1s -> SetMarkerColor( kGray );
    exclusion_m1s -> SetMarkerColor( kGray );
    exclusion_obs -> SetMarkerColor( kRed );
    exclusion_exp -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineColor( kBlack );
    exclusion_p1s -> SetLineStyle( 2 );
    exclusion_m1s -> SetLineColor( kBlack );
    exclusion_m1s -> SetLineStyle( 2 );
    exclusion_obs -> SetLineColor( kRed );

    exclusion_obs_theo_p1s -> SetLineStyle( 2 );
    exclusion_obs_theo_m1s -> SetLineStyle( 2 );
    exclusion_obs_theo_p1s -> SetLineColor(kRed);
    exclusion_obs_theo_m1s -> SetLineColor(kRed);

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
    const Int_t     NRGBs         = 5;
    const Int_t     NCont         = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.50, 0.50, 1.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.50, 1.00, 1.00, 0.60, 0.50 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.50, 0.40, 0.50 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);



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


    TH2D* h3 = new TH2D("h3","",100,1,1e3,100,1,2e2);
    
    h_Limit_obs->GetYaxis()->SetLimits(1,5e+2);
    h_Limit_obs->SetHistogram(h3);
    
    h_Limit_obs->Draw("COLZ");
    exclusion_exp->Draw("L SAMES");
    exclusion_p1s->Draw("L SAMES");
    exclusion_m1s->Draw("L SAMES");
    exclusion_obs->Draw("L SAMES");
    exclusion_obs_theo_p1s->Draw("L SAMES");
    exclusion_obs_theo_m1s->Draw("L SAMES");
    
    //~ diagonal->Draw("L SAMES");

    h_Limit_obs->SetMaximum(50);
    h_Limit_obs->SetMinimum(0.5);

    h_Limit_obs->GetYaxis()->SetTitle("#it{m_{#chi}} [GeV]");
    h_Limit_obs->GetXaxis()->SetTitle("#it{M_{med}} [GeV]");
    h_Limit_obs->GetZaxis()->SetTitle("95% CL observed limit on #sigma_{obs}/#sigma_{theo}");
    h_Limit_obs->GetZaxis()->SetRangeUser(0.01,50);

    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.17,0.37,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Preliminary}}",kBlack);

    if( tag == "MV" ) addText(0.18,0.42,0.835-0.1,0.898-0.1,"Vector, g_{q}=1.0",kBlack);
    if( tag == "MA" ) addText(0.18,0.42,0.835-0.1,0.898-0.1,"Axial Vector, g_{q}=1.0",kBlack);;

    float posx1 = 0.5;
    float posx2 = 0.75;
    float posy1 = 0.75;
    float posy2 = 0.9;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetLineColor(kBlack);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetHeader("#sigma/#sigma_{theo}=1:");
    leg->AddEntry(exclusion_obs, "Observed", "L");
    leg->AddEntry(exclusion_obs_theo_p1s, "Theory Uncertainty", "L");
    leg->AddEntry(exclusion_exp, "Expected", "L");
    leg->AddEntry(exclusion_p1s, "Expected #pm 1#sigma", "L");
    
    
    leg->Draw();
    //~ leg->Draw();
    canv->SaveAs("Limit_Mx_"+tag+"_obs_1p00.png");
    canv->SaveAs("Limit_Mx_"+tag+"_obs_1p00.pdf");

    TFile * f = new TFile("out.root","RECREATE");
    h3->Write();
    h_Limit_obs-> Write();
    exclusion_obs -> Write();
    f->Write();
    f->Close();
    delete canv;

}
