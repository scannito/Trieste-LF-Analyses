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

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/Plot.h"

#include "LFInvMassFitter.h"

int main(int argc, char* argv[])
{
    TFile* outFile = new TFile(mOutFileName.c_str(), "RECREATE");
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhiAssocYield[i][j]->Write();
        }
    }
    outFile->Close();

    std::array<TCanvas*, nbin_deltay> cPhiK0SYield = PlotHistograms(h1PhiK0SYield, h1PhiK0SYieldMB, mOutPath, "rawSpectrumK0SDY%i2D");
    std::array<TCanvas*, nbin_deltay> cPhiPiTPCYield = PlotHistograms(h1PhiPiTPCYield, h1PhiPiTPCYieldMB, mOutPath, "rawSpectrumPiTPCDY%i2D");
    std::array<TCanvas*, nbin_deltay> cPhiPiTOFYield = PlotHistograms(h1PhiPiTOFYield, h1PhiPiTOFYieldMB, mOutPath, "rawSpectrumPiTOFDY%i2D");

    Double_t PhiK0SYieldpTint[nbin_deltay][nbin_mult] = {0}, errPhiK0SYieldpTint[nbin_deltay][nbin_mult] = {0};

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pT::K0S; k++) {
                PhiK0SYieldpTint[i][j] += PhiK0SYieldpTdiff[i][j][k] * (pTK0S_axis[k+1] - pTK0S_axis[k]);
                errPhiK0SYieldpTint[i][j] += TMath::Power(errPhiK0SYieldpTdiff[i][j][k] * (pTK0S_axis[k+1] - pTK0S_axis[k]), 2);
            }
            errPhiK0SYieldpTint[i][j] = TMath::Sqrt(errPhiK0SYieldpTint[i][j]);
        }
    }

    TMultiGraph* mgK0S = new TMultiGraph();
    TGraphAsymmErrors* K0S[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        K0S[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiK0SYieldpTint[i], errmult.data(), errmult.data(), errPhiK0SYieldpTint[i], errPhiK0SYieldpTint[i]);
        if (i == 0) {
            PlotFeatures(K0S[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgK0S);
        } else if (i == 1) {
            PlotFeatures(K0S[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgK0S);
        } else if (i == 2) {
            PlotFeatures(K0S[i], 33, kRed+1, 2, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgK0S);
        }
    }

    TCanvas* cK0S = new TCanvas("cK0S", "cK0S", 1000, 800);
    cK0S->cd();
    //gPad->SetLogy();
    gStyle->SetOptStat(0);
    //gPad->SetMargin(0.16,0.01,0.13,0.06)
    mgK0S->SetTitle("; #LTdN_{ch}/d#eta#GT_{|#eta|<0.5} ; 1/N_{ev,#phi} dN_{K^{0}_{S}}/d#it{y}");
    mgK0S->Draw("AP");

    TLegend* legK0S1 = new TLegend(0.15, 0.82, 0.35, 0.85);
    legK0S1->SetHeader("#bf{This work}");
    legK0S1->SetTextSize(0.05);
    legK0S1->SetLineWidth(0);
    legK0S1->Draw("same");

    TLegend* legK0S2 = new TLegend(0.15, 0.62, 0.35, 0.82);
    legK0S2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legK0S2->AddEntry(K0S[0], "Inclusive", "p");
    legK0S2->AddEntry(K0S[1], "|#it{#Deltay}| < 0.5", "p");
    legK0S2->SetTextSize(0.05);
    legK0S2->SetLineWidth(0);
    legK0S2->Draw("same");

    /*string outNameRawK0SYield = outPath + "rawK0SYield.root";
    cK0S->SaveAs(outNameRawK0SYield.c_str());
    outNameRawK0SYield = outPath + "rawK0SYield.pdf";
    cK0S->SaveAs(outNameRawK0SYield.c_str());*/

    //********************************************************************************************

    Double_t PhiPiYieldpTint[nbin_deltay][nbin_mult] = {0}, errPhiPiYieldpTint[nbin_deltay][nbin_mult] = {0};

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pT::Pi; k++) {
                PhiPiYieldpTint[i][j] += PhiPiYieldpTdiff[i][j][k] * (pTPi_axis[k+1] - pTPi_axis[k]);
                errPhiPiYieldpTint[i][j] += TMath::Power(errPhiPiYieldpTdiff[i][j][k] * (pTPi_axis[k+1] - pTPi_axis[k]), 2);
            }
            errPhiPiYieldpTint[i][j] = TMath::Sqrt(errPhiPiYieldpTint[i][j]);
        }
    }

    TMultiGraph* mgPi = new TMultiGraph();
    TGraphAsymmErrors* Pi[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        Pi[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiPiYieldpTint[i], errmult.data(), errmult.data(), errPhiPiYieldpTint[i], errPhiPiYieldpTint[i]);
        if (i == 0) {
            PlotFeatures(Pi[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgPi);
        } else if (i == 1) {
            PlotFeatures(Pi[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgPi);
        } else if (i == 2) {
            PlotFeatures(Pi[i], 33, kRed+1, 2, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgPi);
        }
    }

    TCanvas* cPi = new TCanvas("cPi", "cPi", 1000, 800);
    cPi->cd();
    //gPad->SetLogy();
    gStyle->SetOptStat(0);
    mgPi->SetTitle("; #LTdN_{ch}/d#eta#GT_{|#eta|<0.5} ; #frac{1}{N_{ev,#phi}} (#pi^{+}+#pi^{#minus})");
    mgPi->Draw("AP");

    TLegend* legPi1 = new TLegend(0.15, 0.82, 0.35, 0.85);
    legPi1->SetHeader("#bf{This work}");
    legPi1->SetTextSize(0.05);
    legPi1->SetLineWidth(0);
    legPi1->Draw("same");

    TLegend* legPi2 = new TLegend(0.15, 0.62, 0.35, 0.82);
    legPi2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legPi2->AddEntry(Pi[0], "Inclusive", "p");
    legPi2->AddEntry(Pi[1], "|#it{#Deltay}| < 0.5", "p");
    legPi2->SetTextSize(0.05);
    legPi2->SetLineWidth(0);
    legPi2->Draw("same");

    /*string outNameRawPiYield = outPath + "rawPiYield.root";
    cPi->SaveAs(outNameRawPiYield.c_str());
    outNameRawPiYield = outPath + "rawPiYield.pdf";
    cPi->SaveAs(outNameRawPiYield.c_str());*/
}
