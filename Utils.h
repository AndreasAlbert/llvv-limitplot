#include <iostream>

#include <sstream>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAxis.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>



void addText(double x1, double x2, double y1, double y2, TString TEXT, Color_t color, Float_t angle = 0)
{
    TPaveText* T = new TPaveText(x1,y1,x2,y2, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    T->SetTextColor(color);
    TText *text = T->AddText(TEXT);
    text->SetTextAngle(angle);
    text->SetTextAlign(22);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);
};

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (double number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

TString Convert2TString (int number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (float number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (double number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TStyle * createTdrStyle() {
    TStyle * tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    //// For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(600); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);

    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    //tdrStyle->SetPadGridX(1);
    //tdrStyle->SetPadGridY(1);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);

    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);

    // For the histo:
    // tdrStyle->SetHistFillColor(1);
    // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    // tdrStyle->SetLegoInnerR(Float_t rad = 0->5);
    // tdrStyle->SetNumberContours(Int_t number = 20);

    tdrStyle->SetEndErrorSize(2);
    //tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(.5);

    tdrStyle->SetMarkerStyle(20);

    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);

    //For the date:
    tdrStyle->SetOptDate(0);
    // tdrStyle->SetDateX(Float_t x = 0->01);
    // tdrStyle->SetDateY(Float_t y = 0->01);

    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(43);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);


    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.16);
    tdrStyle->SetPadRightMargin(0.05);

    // For the Global title:
    tdrStyle->SetTitleFont(43);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);

    // For the axis titles:
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");

    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);

    // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);

    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);


    tdrStyle->cd();
    return tdrStyle;
}


