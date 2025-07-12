#include "Riostream.h"
#include "TApplication.h"
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

int main(int argc, char* argv[])
{
    TApplication app("PhiK0SAnalysis", &argc, argv);

    std::cout << "Starting Phi K0S Analysis..." << std::endl;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.json>" << std::endl;
        return 1;
    }

    std::array<std::array<TH1*, nbin_mult>, nbin_deltay> h1PhiK0SYield;
    std::array<std::array<TH1*, nbin_mult>, nbin_deltay> h1PhiPiTPCYield;
    std::array<std::array<TH1*, nbin_mult>, nbin_deltay> h1PhiPiTOFYield;

    TFile* inFileK0S = TFile::Open("../AnalysisSte/data/PhiK0SYields.root", "READ");
    TFile* inFilePiTPC = TFile::Open("../AnalysisSte/data/PhiPiTPCYields.root", "READ");
    TFile* inFilePiTOF = TFile::Open("../AnalysisSte/data/PhiPiTOFYields.root", "READ");

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhiK0SYield[i][j] = (TH1*)inFileK0S->Get(Form("h1PhiK0SYield_%i_%i", i, j));
            h1PhiK0SYield[i][j]->SetDirectory(0);
            h1PhiPiTPCYield[i][j] = (TH1*)inFilePiTPC->Get(Form("h1PhiPiTPCYield_%i_%i", i, j));
            h1PhiPiTPCYield[i][j]->SetDirectory(0);
            h1PhiPiTOFYield[i][j] = (TH1*)inFilePiTOF->Get(Form("h1PhiPiTOFYield_%i_%i", i, j));
            h1PhiPiTOFYield[i][j]->SetDirectory(0);
        }
    }

    inFileK0S->Close();
    inFilePiTPC->Close();
    inFilePiTOF->Close();
    
    /*std::array<TCanvas*, nbin_deltay> cPhiK0SYield = PlotHistograms(h1PhiK0SYield, h1PhiK0SYieldMB, mOutPath, "rawSpectrumK0SDY%i2D");
    std::array<TCanvas*, nbin_deltay> cPhiPiTPCYield = PlotHistograms(h1PhiPiTPCYield, h1PhiPiTPCYieldMB, mOutPath, "rawSpectrumPiTPCDY%i2D");
    std::array<TCanvas*, nbin_deltay> cPhiPiTOFYield = PlotHistograms(h1PhiPiTOFYield, h1PhiPiTOFYieldMB, mOutPath, "rawSpectrumPiTOFDY%i2D");*/

    std::array<std::array<Double_t, nbin_mult>, nbin_deltay> PhiK0SYieldpTint{}, errPhiK0SYieldpTint{};
    std::array<std::array<Double_t, nbin_mult>, nbin_deltay> PhiPiYieldpTint{}, errPhiPiYieldpTint{};
    std::array<std::array<Double_t, nbin_mult>, nbin_deltay> RatioK0SPiYieldpTint{}, errRatioK0SPiYieldpTint{};

    TH1D* h1PhiPiYield[nbin_deltay][nbin_mult];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhiPiYield[i][j] = new TH1D(Form("h1PhiPiYield_%i_%i", i, j), "; #it{p}_{T} (GeV/#it{c}); 1/N_{ev,#phi} d^{2}N_{#pi}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]", nbin_pT::Pi, pT_axis::Pi.data());
            h1PhiPiYield[i][j]->SetDirectory(0);

            for (int k = 0; k < nbin_pT::Pi; k++) {
                if (k < 4) {
                    if (k == 0) {
                        h1PhiPiYield[i][j]->SetBinContent(k+1, 0.0);
                        h1PhiPiYield[i][j]->SetBinError(k+1, 0.0);
                    } else {
                        h1PhiPiYield[i][j]->SetBinContent(k+1, h1PhiPiTPCYield[i][j]->GetBinContent(k+1));
                        h1PhiPiYield[i][j]->SetBinError(k+1, h1PhiPiTPCYield[i][j]->GetBinError(k+1));
                    }
                } else if (k == 4) {
                    Double_t wTPC = 1.0/ TMath::Power(h1PhiPiTPCYield[i][j]->GetBinError(k+1), 2);
                    Double_t wTOF = 1.0/ TMath::Power(h1PhiPiTOFYield[i][j]->GetBinError(k+1), 2);
                    h1PhiPiYield[i][j]->SetBinContent(k+1, (h1PhiPiTPCYield[i][j]->GetBinContent(k+1) * wTPC + h1PhiPiTOFYield[i][j]->GetBinContent(k+1) * wTOF) / (wTPC + wTOF));
                    h1PhiPiYield[i][j]->SetBinError(k+1, 1.0 / TMath::Sqrt(wTPC + wTOF));
                } else {
                    h1PhiPiYield[i][j]->SetBinContent(k+1, h1PhiPiTOFYield[i][j]->GetBinContent(k+1));
                    h1PhiPiYield[i][j]->SetBinError(k+1, h1PhiPiTOFYield[i][j]->GetBinError(k+1));
                };
            }
        }
    }

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pT::K0S; k++) {
                PhiK0SYieldpTint[i][j] += h1PhiK0SYield[i][j]->GetBinContent(k+1) * (pT_axis::K0S[k+1] - pT_axis::K0S[k]);
                errPhiK0SYieldpTint[i][j] += TMath::Power(h1PhiK0SYield[i][j]->GetBinError(k+1) * (pT_axis::K0S[k+1] - pT_axis::K0S[k]), 2);
            }
            errPhiK0SYieldpTint[i][j] = TMath::Sqrt(errPhiK0SYieldpTint[i][j]);

            for (int k = 0; k < nbin_pT::Pi; k++) {
                PhiPiYieldpTint[i][j] += h1PhiPiYield[i][j]->GetBinContent(k+1) * (pT_axis::Pi[k+1] - pT_axis::Pi[k]);
                errPhiPiYieldpTint[i][j] += TMath::Power(h1PhiPiYield[i][j]->GetBinError(k+1) * (pT_axis::Pi[k+1] - pT_axis::Pi[k]), 2);
            }
            errPhiPiYieldpTint[i][j] = TMath::Sqrt(errPhiPiYieldpTint[i][j]);

            RatioK0SPiYieldpTint[i][j] = 2 * PhiK0SYieldpTint[i][j] / PhiPiYieldpTint[i][j];
            errRatioK0SPiYieldpTint[i][j] = RatioK0SPiYieldpTint[i][j] * TMath::Sqrt(TMath::Power(errPhiK0SYieldpTint[i][j] / PhiK0SYieldpTint[i][j], 2) + TMath::Power(errPhiPiYieldpTint[i][j] / PhiPiYieldpTint[i][j], 2));
        }
    }

    TMultiGraph* mgK0S = new TMultiGraph();
    TGraphAsymmErrors* K0S[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        K0S[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiK0SYieldpTint[i].data(), errmult.data(), errmult.data(), errPhiK0SYieldpTint[i].data(), errPhiK0SYieldpTint[i].data());
        if (i == 0) {
            SetGraphStyle(K0S[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgK0S);
        } else if (i == 1) {
            SetGraphStyle(K0S[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgK0S);
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

    TMultiGraph* mgPi = new TMultiGraph();
    TGraphAsymmErrors* Pi[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        Pi[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiPiYieldpTint[i].data(), errmult.data(), errmult.data(), errPhiPiYieldpTint[i].data(), errPhiPiYieldpTint[i].data());
        if (i == 0) {
            SetGraphStyle(Pi[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgPi);
        } else if (i == 1) {
            SetGraphStyle(Pi[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgPi);
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

    TMultiGraph* mgRatio = new TMultiGraph();
    TGraphAsymmErrors* RatioK0SPi[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        RatioK0SPi[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), RatioK0SPiYieldpTint[i].data(), errmult.data(), errmult.data(), errRatioK0SPiYieldpTint[i].data(), errRatioK0SPiYieldpTint[i].data());
        if (i == 0) {
            SetGraphStyle(RatioK0SPi[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgRatio);
        } else if (i == 1) {
            SetGraphStyle(RatioK0SPi[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgRatio);
        }
    }

    TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1000, 800);
    cRatio->cd();
    //gPad->SetLogy();
    gStyle->SetOptStat(0);
    mgRatio->SetTitle("; #LTdN_{ch}/d#eta#GT_{|#eta|<0.5} ; #frac{2}{N_{ev,#phi}} K^{0}_{S} / (#pi^{+}+#pi^{#minus})");
    mgRatio->Draw("AP");

    TLegend* legRatio1 = new TLegend(0.15, 0.82, 0.35, 0.85);
    legRatio1->SetHeader("#bf{This work}");
    legRatio1->SetTextSize(0.05);
    legRatio1->SetLineWidth(0);
    legRatio1->Draw("same");    

    TLegend* legRatio2 = new TLegend(0.15, 0.62, 0.35, 0.82);   
    legRatio2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legRatio2->AddEntry(RatioK0SPi[0], "Inclusive", "p");
    legRatio2->AddEntry(RatioK0SPi[1], "|#it{#Deltay}| < 0.5", "p");
    legRatio2->SetTextSize(0.05);
    legRatio2->SetLineWidth(0);
    legRatio2->Draw("same");    

    std::cout << "Phi K0S Analysis completed successfully!" << std::endl;

    app.Run();
}
