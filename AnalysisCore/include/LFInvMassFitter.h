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

#include <utility>
#include <vector>
#include <string>

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooProduct.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooProdPdf.h"
#include "RooCBShape.h"
#include "RooCrystalBall.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"

class LFInvMassFitter : public TNamed 
{
    public:
        LFInvMassFitter();
        LFInvMassFitter(const char* name, const char* title);
        ~LFInvMassFitter();

        std::pair<Double_t, Double_t> GetPhiPurityAndError(TH1F* h1PhiInvMass, std::string nameCanvas, Int_t isDataOrReco, Int_t isK0SOrPi, std::vector<Int_t> indices, bool printCanvas = false);
        std::pair<Double_t, Double_t> FitPhiK0S(TH2F* h2PhiK0SInvMass, std::vector<Int_t> indices, TFile* file,
                                                const std::vector<Double_t>& params = {1., 1., 5., 5., 0.49, 0.003, -1.}, 
                                                const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., 0.48, 0.001, -1.7},
                                                const std::vector<Double_t>& upLimits = {2., 2., 10., 10., 0.5, 0.01, 10});
        std::pair<Double_t, Double_t> FitPhiPi(TH2F* h2PhiPiInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file,
                                               const std::vector<Double_t>& params = {1., 1., 10., 10., 0., 1., 7., 1.}, 
                                               const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., -1., 0.001, 3., 0.001}, 
                                               const std::vector<Double_t>& upLimits = {5., 5., 10., 10., 1.5, 2.5, 10., 5.});

        void fitDataAndMCClosure2D(int mode);

    private:
        TH1F* mHisto;

    ClassDef(LFInvMassFitter, 1);
};

