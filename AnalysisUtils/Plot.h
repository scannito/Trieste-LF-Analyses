#pragma once

#include "Riostream.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include <array>
#include <string>

#include "AnalysisUtils/Parameters.h"

inline void SetHistoStyle(TH1* h1, Color_t color)
{
    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(color);
    h1->SetMarkerSize(1.5);
    h1->SetLineColor(color);
    h1->SetLineWidth(2);
    h1->SetFillStyle(3001);
    h1->SetFillColor(color);
    h1->GetYaxis()->SetRangeUser(0.0, 1.2 * h1->GetMaximum());
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleOffset(1.0);
    h1->GetYaxis()->SetLabelSize(0.045);
}

inline void SetGraphStyle(TGraphAsymmErrors* graph, Style_t markstyle, Color_t markcolor, Size_t marksize, Style_t linestyle, Color_t linecolor, Width_t linewidth, Style_t fillstyle, Color_t fillcolor, Float_t alpha, TMultiGraph* mg) 
{													
    graph->SetMarkerStyle(markstyle);
    graph->SetMarkerColor(markcolor);
    graph->SetMarkerSize(marksize);
    graph->SetLineStyle(linestyle);
    graph->SetLineColor(linecolor);
    graph->SetLineWidth(linewidth);
    graph->SetFillStyle(fillstyle);
    graph->SetFillColorAlpha(fillcolor,alpha);
    mg->Add(graph);
}

inline std::array<TCanvas*, nbin_deltay> PlotHistograms(TH1D* h1Yield[nbin_deltay][nbin_mult], TH1D* h1YieldMB[nbin_deltay], std::string outPath, std::string name) 
{
    std::array<TCanvas*, nbin_deltay> cYield;
    TPad* topPad[nbin_deltay]; 
    TPad* bottomPad[nbin_deltay];

    TH1D* h1YieldRatio[nbin_deltay][nbin_mult];

    TLegend* leg1 [nbin_deltay];
    TLegend* leg2 [nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        std::string cName = "c" + name;
        cYield[i] = new TCanvas(Form(cName.c_str(), i), Form(cName.c_str(), i), 800, 800);
        cYield[i]->cd();
        //gPad->SetMargin(0.16,0.01,0.13,0.06);
        gStyle->SetOptStat(0);

        topPad[i] = new TPad("topPad", "Top Pad", 0, bottomPadHeight, 1, 1);
        topPad[i]->SetBottomMargin(0);
        topPad[i]->SetLogy();
        topPad[i]->Draw();

        bottomPad[i] = new TPad("bottomPad", "Bottom Pad", 0, 0, 1, bottomPadHeight);
        bottomPad[i]->SetTopMargin(0);
        bottomPad[i]->SetBottomMargin(0.3);
        bottomPad[i]->SetLogy();
        bottomPad[i]->Draw();

        topPad[i]->cd();
        
        leg1[i] = new TLegend(0.5, 0.82, 0.8, 0.85);
        leg1[i]->SetHeader("#bf{This work}");
        leg1[i]->SetTextSize(0.05);
        leg1[i]->SetLineWidth(0);

        leg2[i] = new TLegend(0.5, 0.62, 0.8, 0.82);
        leg2[i]->SetHeader(Form("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < %f", deltay_axis[i]));
        leg2[i]->SetTextSize(0.035);
        leg2[i]->SetLineWidth(0);
        leg2[i]->SetNColumns(2);
        
        for (int j = 0; j < nbin_mult; j++) {
            h1Yield[i][j]->SetMarkerStyle(20);
            h1Yield[i][j]->SetMarkerColor(Colors[j]);
            h1Yield[i][j]->SetMarkerSize(1.5);
            h1Yield[i][j]->SetLineColor(Colors[j]);
            h1Yield[i][j]->SetLineWidth(2);
            h1Yield[i][j]->SetFillStyle(3001);
            h1Yield[i][j]->SetFillColor(Colors[j]);
            h1Yield[i][j]->GetXaxis()->SetLabelOffset(0.5);
            h1Yield[i][j]->GetYaxis()->SetTitleSize(0.045);
            h1Yield[i][j]->GetYaxis()->SetTitleOffset(1.0);
            h1Yield[i][j]->GetYaxis()->SetLabelSize(0.045);
            h1Yield[i][j]->GetYaxis()->SetRangeUser(1e-3, 1e1);
            h1Yield[i][j]->GetYaxis()->SetRangeUser(0.4e-5, 1.3e-1);

            if (j == 0) h1Yield[i][j]->Draw();
            else h1Yield[i][j]->Draw("same");

            leg2[i]->AddEntry(h1Yield[i][j], Form("%i-%i %%", (int)mult_axis[j], (int)mult_axis[j+1]), "p");
        }

        leg1[i]->Draw("same");
        leg2[i]->Draw("same");

        bottomPad[i]->cd();

        for (int j = 0; j < nbin_mult; j++) {
            h1YieldRatio[i][j] = (TH1D*)h1Yield[i][j]->Clone(Form("h1YieldRatio%i_%i", i, j));
            h1YieldRatio[i][j]->Divide(h1YieldMB[i]);
            h1YieldRatio[i][j]->SetTitle("; #it{p}_{T} (GeV/#it{c}); Ratio to 0-100 %");
            h1YieldRatio[i][j]->SetMarkerStyle(20);
            h1YieldRatio[i][j]->SetMarkerColor(Colors[j]);
            h1YieldRatio[i][j]->SetMarkerSize(1.5);
            h1YieldRatio[i][j]->SetLineColor(Colors[j]);
            h1YieldRatio[i][j]->SetLineWidth(2);
            h1YieldRatio[i][j]->SetFillStyle(3001);
            h1YieldRatio[i][j]->SetFillColor(Colors[j]);
            h1YieldRatio[i][j]->GetXaxis()->SetLabelOffset(0.03);
            h1YieldRatio[i][j]->GetXaxis()->SetNdivisions(515);
            h1YieldRatio[i][j]->GetXaxis()->SetLabelSize(0.045 / scaleFactor);
            h1YieldRatio[i][j]->GetXaxis()->SetTitleSize(0.045 / scaleFactor);
            h1YieldRatio[i][j]->GetXaxis()->SetTitleOffset(1.2);
            h1YieldRatio[i][j]->GetYaxis()->SetTitleSize(0.045 / scaleFactor - 0.01);
            h1YieldRatio[i][j]->GetYaxis()->SetTitleOffset(0.45);
            h1YieldRatio[i][j]->GetYaxis()->SetLabelSize(0.045 / scaleFactor);
            h1YieldRatio[i][j]->GetYaxis()->SetRangeUser(0.5e-1, 1.5e1);

            if (j == 0) h1YieldRatio[i][j]->Draw();
            else h1YieldRatio[i][j]->Draw("same");
        }

        std::string outName = outPath + name + ".root";
        cYield[i]->SaveAs(Form(outName.c_str(), i));
        outName = outPath + name + ".pdf";
        cYield[i]->SaveAs(Form(outName.c_str(), i));
    }

    return cYield;
}