#pragma once

#include "Riostream.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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
#include <array>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
//#include<type_traits>

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

#include "AnalysisUtils/Parameters.h"

//template <typename H, typename = std::enable_if_t<std::is_base_of_v<TH1, H>>>
class LFInvMassFitter : public TNamed 
{
public:
    enum AssocParticleType { K0S, Pi, Unknown };

    LFInvMassFitter() = default;
    LFInvMassFitter(const std::string& assoc, const std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay>& Histo2D, 
                    /*const std::array<std::vector<TH2*>, nbin_deltay>& HistoMultInt2D,*/
                    const std::string& OutPath, const std::string& OutFileName, int mode = 0);
    LFInvMassFitter(const std::string& assoc, const char* filename, int nbin_pT, const std::string& histoname);
    ~LFInvMassFitter();

    std::pair<Double_t, Double_t> GetPhiPurityAndError(TH1* h1PhiInvMass, std::string nameCanvas, Int_t isDataOrReco, Int_t isK0SOrPi, std::vector<Int_t> indices, bool printCanvas = false);
    
    std::pair<Double_t, Double_t> FitPhiAssoc(TH2* h2PhiAssocInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file/*,
                                                const std::vector<Double_t>& params, const std::vector<Double_t>& lowLimits, const std::vector<Double_t>& upLimits*/);

    void ExportYields(const char* filename, Int_t nbin_pT, const std::vector<Double_t>& pT_axis, const std::string& hSetName, Int_t isTPCOrTOF);

    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> mSetHisto2D{};

    void CheckValidMembers();

private:
    AssocParticleType mAssocParticleType{K0S}; // Default association particle type

    //H* mHisto;
    //std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> mSetHisto2D{};
    //std::array<std::vector<TH2*>, nbin_deltay> mSetHistoMultInt2D{};

    //std::string mOutPath;
    //std::string mOutFileName;

    int mNEvents{0}; // Number of processed events
    int mMode{0}; // 0: Data, 1: Closure test

    AssocParticleType StringToAssoc(const std::string& assoc);
    std::string AssocToSymbol(AssocParticleType assoc);

    void HistogramAcquisition(const char* filename, int nbin_pT, const std::string& histoname);

    std::pair<Double_t, Double_t> FitPhiK0S(TH2* h2PhiK0SInvMass, std::vector<Int_t> indices, TFile* file, Double_t nsigma = 6.0,
                                            const std::vector<Double_t>& params = {1., 1., 5., 5., 0.49, 0.003, -1.}, 
                                            const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., 0.48, 0.001, -1.7},
                                            const std::vector<Double_t>& upLimits = {2., 2., 10., 10., 0.5, 0.01, 10});
    std::pair<Double_t, Double_t> FitPhiPi(TH2* h2PhiPiInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file, Double_t nsigma = 4.0,
                                            const std::vector<Double_t>& params = {1., 1., 10., 10., 0., 1., 7., 1.}, 
                                            const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., -1., 0.001, 3., 0.001}, 
                                            const std::vector<Double_t>& upLimits = {5., 5., 10., 10., 1.5, 2.5, 10., 5.});

    ClassDef(LFInvMassFitter, 1);
};

