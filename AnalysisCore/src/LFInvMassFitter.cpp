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

#include <utility>
#include <vector>
#include <string>

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/Plot.h"

#include "LFInvMassFitter.h"

using namespace RooFit;

ClassImp(LFInvMassFitter);

LFInvMassFitter::LFInvMassFitter(const std::string& assoc, const std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay>& Histo2D, 
                                 /*const std::array<std::vector<TH2*>, nbin_deltay>& HistoMultInt2D,*/
                                 const std::string& OutPath, const std::string& OutFileName, int mode) : 
                                 TNamed(), mAssocParticleType(StringToAssoc(assoc)), mSetHisto2D(Histo2D), /*mSetHistoMultInt2D(HistoMultInt2D),*/ 
                                 /*mOutPath(OutPath), mOutFileName(OutFileName),*/ mMode(mode)
{
    std::cout << "LFInvMassFitter initialized for associated particle = " << mAssocParticleType << std::endl;
    if (mAssocParticleType == Unknown) {
        std::cerr << "Error: Unknown associated particle type. Please check the input." << std::endl;
        return;
    }
}

LFInvMassFitter::LFInvMassFitter(const std::string& assoc, const char* filename, int nbin_pT, const std::string& histoname) : 
                                 TNamed(), mAssocParticleType(StringToAssoc(assoc))
{
    std::cout << "LFInvMassFitter initialized for associated particle = " << mAssocParticleType << std::endl;
    HistogramAcquisition(filename, nbin_pT, histoname);
    if (mAssocParticleType == Unknown) {
        std::cerr << "Error: Unknown associated particle type. Please check the input." << std::endl;
        return;
    }
}

LFInvMassFitter::~LFInvMassFitter()
{
    for (auto& histo2DArray : mSetHisto2D) {
        for (auto& histo2D : histo2DArray) {
            for (auto& histo : histo2D) {
                if (histo) 
                    delete histo;
            }
        }
    }
    /*for (auto& histoMultInt : mSetHistoMultInt2D) {
        for (auto& histo : histoMultInt) {
            if (histo) 
                delete histo;
        }
    }*/
}

LFInvMassFitter::AssocParticleType LFInvMassFitter::StringToAssoc(const std::string& assoc) 
{
    std::map<std::string, AssocParticleType> map = {
        {"K0S", K0S},
        {"Pi", Pi}
    };

    std::string key(assoc);
    //std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    auto it = map.find(key);
    return it != map.end() ? it->second : AssocParticleType::Unknown;
}

std::string LFInvMassFitter::AssocToSymbol(AssocParticleType assoc) 
{
    switch (assoc) {
        case K0S: return "K^{0}_{S}";
        case Pi: return "#pi";
        default: return "Unknown";
    }
}

void LFInvMassFitter::HistogramAcquisition(const char* filename, int nbin_pT, const std::string& histoname)
{
    TFile* file = TFile::Open(filename);
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::cout << "Acquiring histograms from file: " << filename << std::endl;

    for (int i = 1; i <= nbin_deltay; ++i) {
        for (int j = 0; j < nbin_mult; ++j) {
            //std::vector<TH2*> histos;
            for (int k = 0; k < nbin_pT; ++k) {
                std::string hName = histoname + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                TH2* htemp = (TH2*)file->Get(hName.data());
                htemp->SetDirectory(0);
                mSetHisto2D[i-1][j].push_back(htemp);
            }
        }
    }

    file->Close();
}

std::pair<Double_t, Double_t> LFInvMassFitter::GetPhiPurityAndError(TH1* h1PhiInvMass, std::string nameCanvas, Int_t isDataOrReco, 
                                                                    Int_t isK0SOrPi, std::vector<Int_t> indices, bool printCanvas)
{
    h1PhiInvMass->SetTitle("; #it{M}(K^{+}K^{#minus}) (GeV/#it{c}^{2}); Counts");
    Double_t binsize = h1PhiInvMass->GetXaxis()->GetBinWidth(1);

    TF1* fitVoigtPolyPur = new TF1("fitVoigtPolyPur", VoigtPoly, 0.987, 1.06, 7);
    if (isDataOrReco == 2) {
        if (isK0SOrPi == 1) {
            fitVoigtPolyPur->SetParameter(0, 10);
            if (indices.size() == 2 && (indices[0] == 0 && (indices[1] == 0 || indices[1] == 5))) fitVoigtPolyPur->SetParameter(0, 100);
            if (indices.size() == 2 && (indices[0] == 2 && indices[1] == 6)) fitVoigtPolyPur->SetParameter(0, 50);
        }
        else if (isK0SOrPi == 2) fitVoigtPolyPur->SetParameter(0, 1.75);
    }
    fitVoigtPolyPur->FixParameter(1, 0.00426);
    fitVoigtPolyPur->SetParameter(2, 1.019);
    fitVoigtPolyPur->SetParameter(3, 0.001);
    fitVoigtPolyPur->SetNpx(400);
    fitVoigtPolyPur->SetLineColor(kRed);

    TFitResultPtr fitResultVoigtPolyPur = h1PhiInvMass->Fit("fitVoigtPolyPur", "RS0");
    const Double_t* paramsVoigtPolyPur = fitResultVoigtPolyPur->GetParams();

    TF1* fitVoigtPolybisPur = new TF1("fitVoigtPolybisPur", VoigtPoly, 0.987, 1.06, 7);
    fitVoigtPolybisPur->SetParameter(0, paramsVoigtPolyPur[0]);
    fitVoigtPolybisPur->SetParameter(1, paramsVoigtPolyPur[1]);
    fitVoigtPolybisPur->SetParameter(2, paramsVoigtPolyPur[2]);
    fitVoigtPolybisPur->SetParameter(3, paramsVoigtPolyPur[3]);
    fitVoigtPolybisPur->FixParameter(4, 0);
    fitVoigtPolybisPur->FixParameter(5, 0);
    fitVoigtPolybisPur->FixParameter(6, 0);

    TF1* fitVoigtPur = new TF1("fitVoigtPur", Voigt, 0.987, 1.05, 4);
    fitVoigtPur->SetParameter(0, paramsVoigtPolyPur[0]);
    fitVoigtPur->SetParameter(1, paramsVoigtPolyPur[1]);
    fitVoigtPur->SetParameter(2, paramsVoigtPolyPur[2]);
    fitVoigtPur->SetParameter(3, paramsVoigtPolyPur[3]);
    fitVoigtPur->SetNpx(400);
    fitVoigtPur->SetLineColor(kBlue);

    if (printCanvas) {
        TCanvas* fitPhiPur = new TCanvas(nameCanvas.c_str(), nameCanvas.c_str(), 800, 800);
        fitPhiPur->cd();
        gStyle->SetOptStat(0);
        h1PhiInvMass->Draw();
        fitVoigtPolyPur->Draw("same");
        fitVoigtPur->Draw("same");
        
        TLegend* legPhi1Pur = new TLegend(0.55, 0.4, 0.75, 0.43);
        legPhi1Pur->SetHeader("#bf{This work}");
        legPhi1Pur->SetTextSize(0.05);
        legPhi1Pur->SetLineWidth(0);
        legPhi1Pur->Draw("same");

        TLegend* legPhi2Pur = new TLegend(0.55, 0.2, 0.75, 0.4);
        legPhi2Pur->SetHeader("#phi #rightarrow K^{+}K^{#minus}");
        legPhi2Pur->AddEntry(fitVoigtPolyPur, "Voigt + bkg", "l");
        legPhi2Pur->AddEntry(fitVoigtPur, "Voigt", "l");
        legPhi2Pur->SetTextSize(0.05);
        legPhi2Pur->SetLineWidth(0);
        legPhi2Pur->Draw("same");

        fitPhiPur->cd();
        
        TLine* line1 = new TLine(lowmPhiPur, h1PhiInvMass->GetMinimum(), lowmPhiPur, h1PhiInvMass->GetMaximum());
        TLine* line2 = new TLine(upmPhiPur, h1PhiInvMass->GetMinimum(), upmPhiPur, h1PhiInvMass->GetMaximum());

        line1->SetLineColor(kBlack);
        line1->SetLineStyle(kDashed);
        line1->SetLineWidth(2);
        line1->Draw("same");

        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDashed);
        line2->SetLineWidth(2);
        line2->Draw("same");
    }

    Double_t integralVoigtPolyPur = fitVoigtPolyPur->Integral(lowmPhiPur, upmPhiPur) / binsize;
    Double_t integralVoigt1Pur = fitVoigtPur->Integral(lowmPhiPur, upmPhiPur) / binsize;

    TMatrixDSym covMatrixVoigtPolyPur = fitResultVoigtPolyPur->GetCovarianceMatrix();
    Double_t errintegralVoigtPolyPur = fitVoigtPolyPur->IntegralError(lowmPhiPur, upmPhiPur, paramsVoigtPolyPur, covMatrixVoigtPolyPur.GetMatrixArray()) / binsize;
    Double_t errintegralVoigt1Pur = fitVoigtPolybisPur->IntegralError(lowmPhiPur, upmPhiPur, paramsVoigtPolyPur, covMatrixVoigtPolyPur.GetMatrixArray()) / binsize;

    Double_t purityVoigt = integralVoigt1Pur / integralVoigtPolyPur;
    Double_t errpurityVoigt = purityVoigt * TMath::Sqrt(TMath::Power(errintegralVoigt1Pur / integralVoigt1Pur, 2) + TMath::Power(errintegralVoigtPolyPur / integralVoigtPolyPur, 2));

    return std::make_pair(purityVoigt, errpurityVoigt);
}

std::pair<Double_t, Double_t> LFInvMassFitter::FitPhiK0S(TH2* h2PhiK0SInvMass, std::vector<Int_t> indices, TFile* file,
                                                const std::vector<Double_t>& params, const std::vector<Double_t>& lowLimits, const std::vector<Double_t>& upLimits)
{
    // Definisci le variabili x e y
    RooRealVar x("x", "x", h2PhiK0SInvMass->GetXaxis()->GetXmin(), h2PhiK0SInvMass->GetXaxis()->GetXmax());
    RooRealVar y("y", "y", h2PhiK0SInvMass->GetYaxis()->GetXmin(), h2PhiK0SInvMass->GetYaxis()->GetXmax());

    // Converte l'istogramma 2D in un RooDataHist
    RooDataHist data("data", "data", RooArgList(x, y), h2PhiK0SInvMass);

    // Definisci i parametri per la Double Sided Crystal Ball e pol1 per l'asse x
    RooRealVar alpha1CB("alpha1CB", "alpha1CB", params.at(0), lowLimits.at(0), upLimits.at(0));
    RooRealVar alpha2CB("alpha2CB", "alpha2CB", params.at(1), lowLimits.at(1), upLimits.at(1));
    RooRealVar n1CB("n1CB", "n1CB", params.at(2), lowLimits.at(2), upLimits.at(2));
    RooRealVar n2CB("n2CB", "n2CB", params.at(3), lowLimits.at(3), upLimits.at(3));
    RooRealVar meanCB("meanCB", "meanCB", params.at(4), lowLimits.at(4), upLimits.at(4));
    RooRealVar sigmaCB("sigmaCB", "sigmaCB", params.at(5), lowLimits.at(5), upLimits.at(5));
    RooCrystalBall dsCrystalBall("dsCrystalBall", "DoubleSidedCrystalBall1", x, meanCB, sigmaCB, alpha1CB, n1CB, alpha2CB, n2CB);

    RooRealVar a1("a1", "a1", params.at(6), lowLimits.at(6), upLimits.at(6));
    RooPolynomial bkgDSCB("bkgDSCB", "bkgDSCB", x, RooArgList(a1));

    RooRealVar meanV("meanV", "meanV", 1.02, 0.987, 1.2);
    RooRealVar width("width", "width", 0.00426);
    width.setConstant(true);
    RooRealVar sigmaV("sigmaV", "sigmaV", 0.001, 0.0001, 0.01);
    RooVoigtian voigt("voigt", "Voigtian", y, meanV, width, sigmaV);

    RooRealVar b0("b0", "b0", 5, 0, 10);
    RooRealVar b1("b1", "b1", 7, 0, 50);
    RooRealVar b2("b2", "b2", 2, 0, 50);
    RooGenericPdf bkgVoigt("bkgVoigt", "bkgVoigt", "b0 + b1*y + b2*sqrt(y-0.987)", RooArgList(y, b0, b1, b2));
    
    RooProdPdf sigsig("sigsig", "sigsig", RooArgList(dsCrystalBall, voigt));
    RooProdPdf sigbkg("sigbkg", "sigbkg", RooArgList(dsCrystalBall, bkgVoigt));
    RooProdPdf bkgsig("bkgsig", "bkgsig", RooArgList(bkgDSCB, voigt));
    RooProdPdf bkgbkg("bkgbkg", "bkgbkg", RooArgList(bkgDSCB, bkgVoigt));

    RooRealVar nsigsig("nsigsig", "nsigsig", 1000, 0, 1000000);
    RooRealVar nsigbkg("nsigbkg", "nsigbkg", 50000, 0, 2000000);
    RooRealVar nbkgsig("nbkgsig", "nbkgsig", 5000, 0, 250000);
    RooRealVar nbkgbkg("nbkgbkg", "nbkgbkg", 10000, 0, 500000);

    RooAddPdf model("model", "model", RooArgList(sigsig, sigbkg, bkgsig, bkgbkg), RooArgList(nsigsig, nsigbkg, nbkgsig, nbkgbkg));

    /*RooAddPdf* model;
    if (indices.size() == 3 && indices[0] == 2 && (indices[1] == 5 || indices[1] == 6) && indices[2 == 6]) model = new RooAddPdf("model", "model", RooArgList(dsCrystalBall), RooArgList(nsig));
    else model = new RooAddPdf("model", "model", RooArgList(dsCrystalBall, bkgDSCB), RooArgList(nsig, nbkg));*/

    // Fitta il modello ai dati
    RooFitResult* result = model.fitTo(data, Optimize(1), Extended(1), Save(1));

    meanCB.Print();
    sigmaCB.Print();
    alpha1CB.Print();
    n1CB.Print();
    alpha2CB.Print();
    n2CB.Print();
    a1.Print();

    meanV.Print();
    width.Print();
    sigmaV.Print();
    b0.Print();
    b1.Print();
    b2.Print();

    nsigsig.Print();
    nsigbkg.Print();
    nbkgsig.Print();
    nbkgbkg.Print();

    Int_t lowedge = h2PhiK0SInvMass->GetXaxis()->FindFixBin(meanCB.getVal() - 3 * sigmaCB.getVal());
    Int_t upedge = h2PhiK0SInvMass->GetXaxis()->FindFixBin(meanCB.getVal() + 3 * sigmaCB.getVal());

    Double_t lowfitK0S = h2PhiK0SInvMass->GetXaxis()->GetBinLowEdge(lowedge);
    Double_t upfitK0S = h2PhiK0SInvMass->GetXaxis()->GetBinLowEdge(upedge +1);

    TCanvas* cPhiK0SInvMass;
    if (indices.size() == 2) cPhiK0SInvMass = new TCanvas(Form("cPhiK0SInvMass_%d_%d", indices[0], indices[1]), Form("cPhiK0SInvMass_%d_%d", indices[0], indices[1]), 800, 800);
    else if (indices.size() == 3) cPhiK0SInvMass = new TCanvas(Form("cPhiK0SInvMass_%d_%d_%d", indices[0], indices[1], indices[2]), Form("cPhiK0SInvMass_%d_%d_%d", indices[0], indices[1], indices[2]), 800, 800);
    cPhiK0SInvMass->cd();
    gPad->SetMargin(0.16,0.03,0.13,0.06);
    gStyle->SetOptStat(0);

    /*RooPlot* frame = x.frame();
    if (indices.size() == 2) {
        frame->SetName(Form("frame_%d_%d", indices[0], indices[1]));
        frame->SetTitle(Form("Mult int, %f - %f (GeV/#it{c}); #it{M}(#pi^{+} + #pi^{#minus}) (GeV/#it{c}^{2}); Counts", pTK0S_axis[indices[1]], pTK0S_axis[indices[1]+1]));
    }
    else if (indices.size() == 3) {
        frame->SetName(Form("frame_%d_%d_%d", indices[0], indices[1], indices[2]));
        frame->SetTitle(Form("%f - %f %%, %f - %f (GeV/#it{c}); #it{M}(#pi^{+} + #pi^{#minus}) (GeV/#it{c}^{2}); Counts", mult_axis[indices[1]], mult_axis[indices[1]+1], pTK0S_axis[indices[2]], pTK0S_axis[indices[2]+1]));
    }
    data.plotOn(frame);
    model.plotOn(frame);
    //model.plotOn(frame, Components(sigsig), LineColor(kRed), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    //model.plotOn(frame, Components(sigbkg), LineColor(kGreen), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    //model.plotOn(frame, Components(bkgsig), LineColor(kBlue), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    //model.plotOn(frame, Components(bkgbkg), LineColor(kMagenta), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    frame->Draw();*/

    h2PhiK0SInvMass->Draw("LEGO2");
    h2PhiK0SInvMass->SetTitle("; #it{M}(#pi^{+}#pi^{#minus}) (GeV/#it{c}^{2}); #it{M}(K^{+}K^{#minus}) (GeV/#it{c}^{2}); Counts");
    h2PhiK0SInvMass->GetXaxis()->SetTitleOffset(1.7);
    h2PhiK0SInvMass->GetXaxis()->SetLabelSize(0.03);
    h2PhiK0SInvMass->GetYaxis()->SetTitleOffset(2.3);
    h2PhiK0SInvMass->GetYaxis()->SetLabelSize(0.03);
    h2PhiK0SInvMass->GetZaxis()->SetTitleOffset(1.7);
    h2PhiK0SInvMass->GetZaxis()->SetLabelSize(0.03);

    TH2F* hist2D_fit = (TH2F*)model.createHistogram("hist2D_fit", x, Binning(h2PhiK0SInvMass->GetNbinsX()), YVar(y, Binning(h2PhiK0SInvMass->GetNbinsY())));
    hist2D_fit->SetLineColor(kRed);
    hist2D_fit->GetZaxis()->SetTitleOffset(1.5);
    hist2D_fit->Draw("SURF SAME");

    file->cd();
    cPhiK0SInvMass->Write();
    delete cPhiK0SInvMass;

    // Calcola l'integrale della funzione prodotto nel range specificato
    x.setRange("signal", lowfitK0S, upfitK0S);
    RooAbsReal* integralsigsig = sigsig.createIntegral(RooArgSet(x, y), NormSet(x, y), Range("signal"));

    RooProduct sigyield("sigyield", "sigyield", RooArgList(nsigsig, *integralsigsig));
    Double_t PhiK0SYieldpTdiff = sigyield.getVal();
    Double_t errPhiK0SYieldpTdiff = sigyield.getPropagatedError(*result, RooArgSet(x, y));

    return std::make_pair(PhiK0SYieldpTdiff, errPhiK0SYieldpTdiff);
}

std::pair<Double_t, Double_t> LFInvMassFitter::FitPhiPi(TH2* h2PhiPiInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file,
                                               const std::vector<Double_t>& params, const std::vector<Double_t>& lowLimits, const std::vector<Double_t>& upLimits)
{
    // Definisci le variabili x e y
    RooRealVar x("x", "x", h2PhiPiInvMass->GetXaxis()->GetXmin(), h2PhiPiInvMass->GetXaxis()->GetXmax());
    RooRealVar y("y", "y", h2PhiPiInvMass->GetYaxis()->GetXmin(), h2PhiPiInvMass->GetYaxis()->GetXmax());

    // Converte l'istogramma 2D in un RooDataHist
    RooDataHist data("data", "data", RooArgList(x, y), h2PhiPiInvMass);

    // Definisci i parametri per la Double Sided Crystal Ball e pol1 per l'asse x
    RooRealVar alpha1CB("alpha1CB", "alpha1CB", params.at(0), lowLimits.at(0), upLimits.at(0));
    RooRealVar alpha2CB("alpha2CB", "alpha2CB", params.at(1), lowLimits.at(1), upLimits.at(1));
    RooRealVar n1CB("n1CB", "n1CB", params.at(2), lowLimits.at(2), upLimits.at(2));
    RooRealVar n2CB("n2CB", "n2CB", params.at(3), lowLimits.at(3), upLimits.at(3));
    RooRealVar meanCB("meanCB", "meanCB", params.at(4), lowLimits.at(4), upLimits.at(4));
    RooRealVar sigmaCB("sigmaCB", "sigmaCB", params.at(5), lowLimits.at(5), upLimits.at(5));
    RooCrystalBall dsCrystalBall("dsCrystalBall", "DoubleSidedCrystalBall1", x, meanCB, sigmaCB, alpha1CB, n1CB, alpha2CB, n2CB);

    RooRealVar meanG1("meanG1", "meanG1", -7., -8., -6.);
    RooRealVar sigmaG1("sigmaG1", "sigmaG1", 0.2, 0.1, 1.);
    RooGaussian gauss1("gauss1", "gauss1", x, meanG1, sigmaG1);

    RooRealVar meanG2("meanG2", "meanG2", params.at(6), lowLimits.at(6), upLimits.at(6));
    RooRealVar sigmaG2("sigmaG2", "sigmaG2", params.at(7), lowLimits.at(7), upLimits.at(7));
    RooGaussian gauss2("gauss2", "gauss2", x, meanG2, sigmaG2);

    RooRealVar meanG3("meanG3", "meanG3", 2., 0., 7.);
    RooRealVar sigmaG3("sigmaG3", "sigmaG3", 1., 0.001, 5.);
    RooGaussian gauss3("gauss3", "gauss3", x, meanG3, sigmaG3);

    RooRealVar meanG4("meanG4", "meanG4", 4., 2., 6.);
    RooRealVar sigmaG4("sigmaG4", "sigmaG4", 0.2, 0.1, 2.);
    RooGaussian gauss4("gauss4", "gauss4", x, meanG4, sigmaG4);

    RooRealVar meanV("meanV", "meanV", 1.02, 0.987, 1.2);
    RooRealVar width("width", "width", 0.00426);
    width.setConstant(true);
    RooRealVar sigmaV("sigmaV", "sigmaV", 0.001, 0.0001, 0.01);
    RooVoigtian voigt("voigt", "Voigtian", y, meanV, width, sigmaV);

    RooRealVar b0("b0", "b0", 5, 0, 10);
    RooRealVar b1("b1", "b1", 7, 0, 50);
    RooRealVar b2("b2", "b2", 2, 0, 50);
    RooGenericPdf bkgVoigt("bkgVoigt", "bkgVoigt", "b0 + b1*y + b2*sqrt(y-0.987)", RooArgList(y, b0, b1, b2));

    RooProdPdf sigsig("sigsig", "sigsig", RooArgList(dsCrystalBall, voigt));
    RooProdPdf sigbkg("sigbkg", "sigbkg", RooArgList(dsCrystalBall, bkgVoigt));
    RooProdPdf missig1sig("missig1sig", "missig1sig", RooArgList(gauss1, voigt));
    RooProdPdf missig1bkg("missig1bkg", "missig1bkg", RooArgList(gauss1, bkgVoigt));
    RooProdPdf missig2sig("missig2sig", "missig2sig", RooArgList(gauss2, voigt));
    RooProdPdf missig2bkg("missig2bkg", "missig2bkg", RooArgList(gauss2, bkgVoigt));
    RooProdPdf missig3sig("missig3sig", "missig3sig", RooArgList(gauss3, voigt));
    RooProdPdf missig3bkg("missig3bkg", "missig3bkg", RooArgList(gauss3, bkgVoigt));
    RooProdPdf missig4sig("missig4sig", "missig4sig", RooArgList(gauss4, voigt));
    RooProdPdf missig4bkg("missig4bkg", "missig4bkg", RooArgList(gauss4, bkgVoigt));

    RooRealVar nsigsig("nsigsig", "nsigsig", 100000, 0, 90000000);
    RooRealVar nsigbkg("nsigbkg", "nsigbkg", 50000, 0, 50000000);
    RooRealVar nmissig1sig("nmissig1sig", "nmissig1sig", 100000, 0, 90000000);
    RooRealVar nmissig1bkg("nmissig1bkg", "nmissig1bkg", 50000, 0, 50000000);
    RooRealVar nmissig2sig("nmissig2sig", "nmissig2sig", 10000, 0, 900000);
    RooRealVar nmissig2bkg("nmissig2bkg", "nmissig2bkg", 5000, 0, 500000);
    RooRealVar nmissig3sig("nmissig3sig", "nmissig3sig", 10000, 0, 900000);
    RooRealVar nmissig3bkg("nmissig3bkg", "nmissig3bkg", 5000, 0, 500000);
    RooRealVar nmissig4sig("nmissig4sig", "nmissig4sig", 10000, 0, 900000);
    RooRealVar nmissig4bkg("nmissig4bkg", "nmissig4bkg", 5000, 0, 500000);

    //RooAddPdf model("model", "model", RooArgList(sig1, sig2), RooArgList(nsig1, nsig2));
    RooAddPdf* model;//("model", "model", RooArgList(dsCrystalBall, gauss), RooArgList(nsig1, nsig2));
    if (isDataOrMcReco == 0) {
        if (isTPCOrTOF == 0) {
            if (indices.size() == 2) {
                if (indices[1] < 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig3sig, sigbkg, missig3bkg), RooArgList(nsigsig, nmissig3sig, nsigbkg, nmissig3bkg));
            } else if (indices.size() == 3) {
                if (indices[2] < 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig3sig, sigbkg, missig3bkg), RooArgList(nsigsig, nmissig3sig, nsigbkg, nmissig3bkg));
            }
        } else if (isTPCOrTOF == 1) {
            if (indices.size() == 2){
                if (indices[1] < 6) {
                    if (indices[1] == 1) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig1sig, sigbkg, missig1bkg), RooArgList(nsigsig, nmissig1sig, nsigbkg, nmissig1bkg));
                    else model = new RooAddPdf("model", "model", RooArgList(sigsig, sigbkg), RooArgList(nsigsig, nsigbkg));
                }
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
            } else if (indices.size() == 3) {
                if (indices[2] < 6) {
                    if (indices[2] == 1) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig1sig, sigbkg, missig1bkg), RooArgList(nsigsig, nmissig1sig, nsigbkg, nmissig1bkg));
                    else model = new RooAddPdf("model", "model", RooArgList(sigsig, sigbkg), RooArgList(nsigsig, nsigbkg));
                }
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                if (indices[0] == 2 && indices[1] == 9 && (indices[2] == 4 || indices[2] == 5)) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
            }
        }
    } else if (isDataOrMcReco == 1) {
        if (isTPCOrTOF == 0) {
            if (indices.size() == 2) {
                if (indices[1] < 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig3sig, sigbkg, missig3bkg), RooArgList(nsigsig, nmissig3sig, nsigbkg, nmissig3bkg));
            } else if (indices.size() == 3) {
                if (indices[2] < 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig3sig, sigbkg, missig3bkg), RooArgList(nsigsig, nmissig3sig, nsigbkg, nmissig3bkg));
            }
        } else if (isTPCOrTOF == 1) {
            if (indices.size() == 2){
                if (indices[1] < 3) {
                    if (indices[1] == 1) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig1sig, sigbkg, missig1bkg), RooArgList(nsigsig, nmissig1sig, nsigbkg, nmissig1bkg));
                    else model = new RooAddPdf("model", "model", RooArgList(sigsig, sigbkg), RooArgList(nsigsig, nsigbkg));
                }
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                if (indices[0] == 2 && indices[1] == 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, missig4sig, sigbkg, missig2bkg, missig4bkg), RooArgList(nsigsig, nmissig2sig, nmissig4sig, nsigbkg, nmissig2bkg, nmissig4bkg));
            } else if (indices.size() == 3) {
                if (indices[2] < 3) {
                    if (indices[2] == 1) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig1sig, sigbkg, missig1bkg), RooArgList(nsigsig, nmissig1sig, nsigbkg, nmissig1bkg));
                    else model = new RooAddPdf("model", "model", RooArgList(sigsig, sigbkg), RooArgList(nsigsig, nsigbkg));
                }
                else model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, sigbkg, missig2bkg), RooArgList(nsigsig, nmissig2sig, nsigbkg, nmissig2bkg));
                if (indices[0] == 0 && (indices[1] == 8 || indices[1] == 9) && indices[2] == 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, missig4sig, sigbkg, missig2bkg, missig4bkg), RooArgList(nsigsig, nmissig2sig, nmissig4sig, nsigbkg, nmissig2bkg, nmissig4bkg));
                if (indices[0] == 1 && (indices[1] == 8 || indices[1] == 9) && indices[2] == 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, missig4sig, sigbkg, missig2bkg, missig4bkg), RooArgList(nsigsig, nmissig2sig, nmissig4sig, nsigbkg, nmissig2bkg, nmissig4bkg));
                if (indices[0] == 2 && (indices[1] == 5 || indices[1] == 6 ||
                    indices[1] == 7 || indices[1] == 8 || indices[1] == 9) && indices[2] == 6) model = new RooAddPdf("model", "model", RooArgList(sigsig, missig2sig, missig4sig, sigbkg, missig2bkg, missig4bkg), RooArgList(nsigsig, nmissig2sig, nmissig4sig, nsigbkg, nmissig2bkg, nmissig4bkg));
            }
        }
    }

    // Fitta il modello ai dati
    RooFitResult* result = model->fitTo(data, Optimize(1), Extended(1), Save(1));

    meanCB.Print();
    sigmaCB.Print();
    alpha1CB.Print();
    n1CB.Print();
    alpha2CB.Print();
    n2CB.Print();

    Int_t lowedge = h2PhiPiInvMass->GetXaxis()->FindFixBin(meanCB.getVal() - 3 * sigmaCB.getVal());
    Int_t upedge = h2PhiPiInvMass->GetXaxis()->FindFixBin(meanCB.getVal() + 3 * sigmaCB.getVal());

    Double_t lowfitPi = h2PhiPiInvMass->GetXaxis()->GetBinLowEdge(lowedge);
    Double_t upfitPi = h2PhiPiInvMass->GetXaxis()->GetBinLowEdge(upedge +1);

    TCanvas* cPhiPiInvMass;
    if (indices.size() == 2) cPhiPiInvMass = new TCanvas(Form("cPhiPiInvMass_%d_%d", indices[0], indices[1]), Form("cPhiPiInvMass_%d_%d", indices[0], indices[1]), 800, 800);
    else if (indices.size() == 3) cPhiPiInvMass = new TCanvas(Form("cPhiPiInvMass_%d_%d_%d", indices[0], indices[1], indices[2]), Form("cPhiPiInvMass_%d_%d_%d", indices[0], indices[1], indices[2]), 800, 800);
    cPhiPiInvMass->cd();
    gPad->SetMargin(0.16,0.03,0.13,0.06);
    gStyle->SetOptStat(0);

    /*RooPlot* frame = x.frame();
    frame->SetTitle("; #it{M}(#pi^{#pm} + #pi^{#mp}) (GeV/#it{c}^{2}); Counts");
    if (indices.size() == 2) {
        frame->SetName(Form("frame_%d_%d", indices[0], indices[1]));
        frame->SetTitle(Form("Mult int, %f - %f (GeV/#it{c}); n#sigma(#pi^{#pm}); Counts", pTPi_axis[indices[1]], pTPi_axis[indices[1]+1]));
    }
    else if (indices.size() == 3) {
        frame->SetName(Form("frame_%d_%d_%d", indices[0], indices[1], indices[2]));
        frame->SetTitle(Form("%f - %f %%, %f - %f (GeV/#it{c}); n#sigma(#pi^{#pm}); Counts", mult_axis[indices[1]], mult_axis[indices[1]+1], pTPi_axis[indices[2]], pTPi_axis[indices[2]+1]));
    }
    data.plotOn(frame);
    model->plotOn(frame);
    model->plotOn(frame, Components(dsCrystalBall), LineColor(kRed), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components(gauss1), LineColor(kGreen), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components(gauss2), LineColor(kGreen), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components(gauss3), LineColor(kGreen), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components(gauss4), LineColor(kOrange), LineStyle(kSolid), Normalization(1.0, RooAbsReal::RelativeExpected));
    frame->Draw();*/      

    h2PhiPiInvMass->Draw("LEGO2");
    h2PhiPiInvMass->SetTitle("; n#sigma(#pi^{#pm}); #it{M}(K^{+}K^{#minus}) (GeV/#it{c}^{2}); Counts");
    h2PhiPiInvMass->GetXaxis()->SetTitleOffset(1.7);
    h2PhiPiInvMass->GetXaxis()->SetLabelSize(0.03);
    h2PhiPiInvMass->GetYaxis()->SetTitleOffset(2.3);
    h2PhiPiInvMass->GetYaxis()->SetLabelSize(0.03);
    h2PhiPiInvMass->GetZaxis()->SetTitleOffset(1.7);
    h2PhiPiInvMass->GetZaxis()->SetLabelSize(0.03);

    TH2F* hist2D_fit = (TH2F*)model->createHistogram("hist2D_fit", x, Binning(h2PhiPiInvMass->GetNbinsX()), YVar(y, Binning(h2PhiPiInvMass->GetNbinsY())));
    hist2D_fit->SetLineColor(kRed);
    hist2D_fit->GetZaxis()->SetTitleOffset(1.5);
    hist2D_fit->Draw("SURF SAME");

    file->cd();
    cPhiPiInvMass->Write();
    delete cPhiPiInvMass;

    // Calcola l'integrale della funzione prodotto nel range specificato
    x.setRange("signal", lowfitPi, upfitPi);
    RooAbsReal* integralsigsig = sigsig.createIntegral(RooArgSet(x, y), NormSet(x, y), Range("signal"));

    RooProduct sigyield("sigyield", "sigyield", RooArgList(nsigsig, *integralsigsig));
    Double_t PhiPiYieldpTdiff = sigyield.getVal();
    Double_t errPhiPiYieldpTdiff = sigyield.getPropagatedError(*result, RooArgSet(x, y));

    return std::make_pair(PhiPiYieldpTdiff, errPhiPiYieldpTdiff);
}

std::pair<Double_t, Double_t> LFInvMassFitter::FitPhiAssoc(TH2* h2PhiAssocInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file/*,
                                                           const std::vector<Double_t>& params, const std::vector<Double_t>& lowLimits, const std::vector<Double_t>& upLimits*/)
{
    switch (mAssocParticleType) {
        case K0S: // K0S
            return FitPhiK0S(h2PhiAssocInvMass, indices, file/*, params, lowLimits, upLimits*/);
        case Pi: // Pi
            return FitPhiPi(h2PhiAssocInvMass, indices, isTPCOrTOF, isDataOrMcReco, file/*, params, lowLimits, upLimits*/);
        default:
            throw std::invalid_argument("Invalid association type specified.");
    }
}

void LFInvMassFitter::ExportYields(const char* filename, Int_t nbin_pT, const std::vector<Double_t>& pT_axis, const std::string& hSetName, Int_t isTPCOrTOF) 
{
    TFile* outputFile = TFile::Open(filename, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error opening output file: " << filename << std::endl;
        return;
    }
    
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(j);
            TH1* h1PhiAssocYield = new TH1D(hName.data(), hName.data(), nbin_pT, pT_axis.data());       
            h1PhiAssocYield->SetTitle(Form("; #it{p}_{T} (GeV/#it{c}); 1/N_{ev,#phi} d^{2}N_{%s}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]", AssocToSymbol(mAssocParticleType).data()));

            for (int k = 0; k < nbin_pT; k++) {
                std::cout << "Processing bin: " << i << ", " << j << ", " << k << std::endl;
                auto [PhiAssocYieldpTdiff, errPhiAssocYieldpTdiff] = FitPhiAssoc(mSetHisto2D[i][j][k], {i, j, k}, isTPCOrTOF, mMode, outputFile);
                PhiAssocYieldpTdiff /= deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pT_axis[k+1] - pT_axis[k]) /*/ mNEvents*/;
                errPhiAssocYieldpTdiff /= deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pT_axis[k+1] - pT_axis[k]) /*/ mNEvents*/;

                h1PhiAssocYield->SetBinContent(k + 1, PhiAssocYieldpTdiff);
                h1PhiAssocYield->SetBinError(k + 1, errPhiAssocYieldpTdiff);
            }

            outputFile->cd();
            h1PhiAssocYield->Write();
            delete h1PhiAssocYield;
        }
    }

    outputFile->Close();
}

void LFInvMassFitter::CheckValidMembers()
{
    int id = 0;
    for (auto& histo2DArray : mSetHisto2D) {
        for (auto& histo2D : histo2DArray) {
            for (auto& histo : histo2D) {
                if (!histo)
                    std::cerr << "Error: One of the TH2 histograms in mSetHisto2D is not initialized." << std::endl;
                else 
                    std::cout << id << ": TH2 histogram is valid: " << histo << std::endl;
                id++;
            }
        }
    }
}