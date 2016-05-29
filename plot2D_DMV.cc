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

TGraph *InterpolateDM(TString tag, TGraph2D* h_Limit, double sqrt_gxgq=1.0, double MinMV=1, double MaxMV=2000, double MaxDM=1000)
{
    std::cout << "Start interpolation." << std::endl;
    TGraph *h_ =  new TGraph();
    TRandom *r0 = new TRandom();

    ofstream myfile;
    myfile.open ("interpolate_"+tag+".txt");

    int maxTries = 50000;
    int nnpts=0;
    int Nsteps = 1000;
    int stepsWithoutHit = 0;
    
    for(int j=0; j<Nsteps; j++) {
        double DMmass = j*(MaxDM-0)/Nsteps;
        for(int i=0; i<maxTries; i++) {
            double Lambda=r0->Uniform(MinMV,MaxMV);
            double mu_val = h_Limit->Interpolate(DMmass,Lambda);
            if(fabs(mu_val-sqrt_gxgq)<1e-3) {
                h_->SetPoint(nnpts,DMmass,Lambda);
                nnpts++;
                myfile << DMmass << " " << Lambda << endl;
                break;
            }
            if( i == maxTries - 1 ) stepsWithoutHit++;
        }
        if( stepsWithoutHit > 10 ) break;
    }
    myfile.close();
    std::cout << "Finish interpolation." << std::endl;
    return h_;
}


void plot2D_DMV(TString myfolder = "")
{


    TString tag = "";
    if(myfolder.Contains("/MV")) tag="MV";
    if(myfolder.Contains("/MA")) tag="MA";

    bool showObserved = true;
    
    vector<int> dm_masses = {1,10,30,40,50,60,70,80,90,100,150,500,1000};


    /// Map:    DM Mass -> Mediator masses
    std::map<int, std::vector<int>> MassMap;

    if( tag == "MV" ) {
        MassMap[ 1   ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 5   ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 10  ] = {10,20,100,5000};
        MassMap[ 20  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 30  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 40  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 50  ] = {10,50,95,200,300,5000};
        MassMap[ 60  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 70  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 80  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 90  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 100 ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 150 ] = {10,200,295,500,5000};
        MassMap[ 500 ] = {995,2000,5000};
        MassMap[ 1000] = {10,1000,1995,5000};
    } else if( tag == "MA" ) {
        MassMap[ 1   ] = {10,20,200,300,500,2000,5000};
        MassMap[ 5   ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 10  ] = {10,20,100,5000};
        MassMap[ 20  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 30  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 40  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 50  ] = {10,50,95,200,300};
        MassMap[ 60  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 70  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 80  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 90  ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 100 ] = {10,20,50,100,200,300,500,1000,2000,5000};
        MassMap[ 150 ] = {10,200,295,500,5000};
        MassMap[ 500 ] = {10,995,2000,5000};
        MassMap[ 1000] = {10,1000,1995,5000};
    }

    TGraph2D* h_Limit = new TGraph2D();
    TGraph2D* h_Limit_p1s = new TGraph2D();
    TGraph2D* h_Limit_m1s = new TGraph2D();

    TGraph2D* h_Limit_t = new TGraph2D();

    int npoint=0;
    for(int nmx=0; nmx<dm_masses.size(); nmx++) {

        std::vector<int> mv_masses = MassMap[ dm_masses[nmx] ];

        int sizemv = mv_masses.size();

        for(int nmv=0; nmv<sizemv; nmv++) {

            TString str_mx = Convert2TString(dm_masses[nmx]);
            TString str_mv = Convert2TString(mv_masses[nmv]);

            TFile myfile( myfolder+"/"+str_mx+"_"+str_mv+"/higgsCombineZwimps01jets.Asymptotic.mH1.root" );

            Double_t obs;
            Double_t m2s;
            Double_t m1s;
            Double_t exp;
            Double_t p1s;
            Double_t p2s;

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


            double sqrt_gxgq = exp;
            cout << "mx: " << dm_masses[nmx] << " mv: " << mv_masses[nmv] << " sqrt_gxgq: " << sqrt_gxgq << endl;


            h_Limit->SetPoint(npoint,dm_masses[nmx],mv_masses[nmv],sqrt_gxgq);
            h_Limit_t->SetPoint(npoint,mv_masses[nmv],dm_masses[nmx],sqrt_gxgq);
            
            h_Limit_p1s->SetPoint(npoint,dm_masses[nmx],mv_masses[nmv],p1s);
            h_Limit_m1s->SetPoint(npoint,dm_masses[nmx],mv_masses[nmv],m1s);
            npoint++;
        }


    }
    std::cout << "Reading Finished." << std::endl;
    TGraph*   h_g1 =  InterpolateDM(tag,h_Limit,1.0,1.,500,200);
    h_g1->SetLineWidth(3);
    h_g1->SetMarkerStyle(20);
    h_g1->SetLineStyle(4);
    h_g1->SetMarkerSize(0.5);
    h_g1->SetMarkerColor(2);
    h_g1->SetLineColor(2);
    
    TGraph*   h_g1_t = (TGraph*)  h_g1->Clone();
    Double_t x, y;
    for( int i = 1; i < h_g1->GetN(); i++ ) {
        h_g1->GetPoint(i,x,y);
        h_g1_t->SetPoint(i,y,x);
    }
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
    //////////// MMED VS MCHI PLOT
    std::cout << "Create Canvas and pad." << std::endl;
   
    TCanvas *canv = new TCanvas("c", "c", 700,600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);

    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogx(true);
    t1->SetLogy(true);
    t1->SetLogz(true);
    t1->SetRightMargin(0.2);

    
    TH2D* h2 = new TH2D("h2","",1000,1,1e3,1000,1,1e5);


    float posx1 = 0.48;
    float posx2 = 0.78;
    float posy1 = 0.75;
    float posy2 = 0.90;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->AddEntry(h_g1,"#sigma_{obs}/#sigma_{theo}=1","PL");
    
    h_Limit->GetYaxis()->SetLimits(1,1e+5);
    h_Limit->SetHistogram(h2);
    h_Limit->Draw("COLZ");
    h_g1->Draw("P sames");
    

    h_Limit->SetMaximum(50);
    h_Limit->SetMinimum(0.5);

    h_Limit->GetXaxis()->SetTitle("#it{m_{#chi}} [GeV]");
    h_Limit->GetYaxis()->SetTitle("#it{M_{med}} [GeV]");

    h_Limit->GetZaxis()->SetTitle("90% CL expected limit on #sigma_{obs}/#sigma_{theo}");
    h_Limit->GetZaxis()->SetRangeUser(0.01,50);


    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.17,0.37,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);

    if( tag == "MV" ) addText(0.13,0.37,0.835-0.1,0.898-0.1,"Vector",kBlack);
    if( tag == "MA" ) addText(0.13,0.37,0.835-0.1,0.898-0.1,"Axial Vector",kBlack);

    canv->cd();
    leg->Draw();
    canv->SaveAs("Limit_Mx_"+tag+".png");
    canv->SaveAs("Limit_Mx_"+tag+".pdf");
    canv->Clear();


    ///////////////////////////////////////////
    //////////// TRANSPOSED PLOT
    std::cout << "Create Canvas and pad." << std::endl;
    canv = new TCanvas("c2", "c2", 700, 600);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);

    t1 = new TPad("t2","t2", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogx(true);
    t1->SetLogy(true);
    t1->SetLogz(true);
    t1->SetRightMargin(0.2);


    TH2D* h3 = new TH2D("h3","",1000,1,1e5,1000,1,1e3);
    
    h_Limit_t->GetYaxis()->SetLimits(1,1e+3);
    h_Limit_t->SetHistogram(h3);
    
    h_Limit_t->Draw("COLZ");
    h_g1->Draw("P sames");

    h_Limit_t->SetMaximum(50);
    h_Limit_t->SetMinimum(0.5);

    h_Limit_t->GetYaxis()->SetTitle("#it{m_{#chi}} [GeV]");
    h_Limit_t->GetXaxis()->SetTitle("#it{M_{med}} [GeV]");
    h_Limit_t->GetZaxis()->SetTitle("90% CL expected limit on #sigma_{obs}/#sigma_{theo}");
    h_Limit_t->GetZaxis()->SetRangeUser(0.01,50);

    
    h_Limit_t ->Draw("COLZ");
    h_g1_t    ->Draw("P sames");

    addText(0.7-0.15,0.995-0.15,0.94,0.996,"2.3 fb^{-1} (13 TeV)",kBlack);
    addText(0.17,0.37,0.835+0.01,0.898+0.01,"#splitline{#bf{CMS}}{#it{Work in Progress}}",kBlack);

    if( tag == "MV" ) addText(0.13,0.37,0.835-0.1,0.898-0.1,"Vector",kBlack);
    if( tag == "MA" ) addText(0.13,0.37,0.835-0.1,0.898-0.1,"Axial Vector",kBlack);;
    

    leg->Draw();
    canv->SaveAs("Limit_Mx_"+tag+"_t.png");
    canv->SaveAs("Limit_Mx_"+tag+"_t.pdf");


    TFile * outfile = new TFile( "Limit_Mx_"+tag+".root","RECREATE" );
    outfile->cd();
    h_Limit->Write();
    h_g1->Write();
    h2->Write();
    h_Limit_t->Write();
    h_g1_t->Write();
    h3->Write();
    canv->Write();
    outfile->Close();
    
    delete canv;

}
