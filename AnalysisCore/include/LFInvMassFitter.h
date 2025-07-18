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
#include <unordered_map>
#include <algorithm>
//#include<type_traits>

#include "RooWorkspace.h"

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/AxesUtils.h"
#include "AnalysisUtils/ParticleTypes.h"

struct Yield
{
    double  value;
    double  error;
};

//template <typename H, typename = std::enable_if_t<std::is_base_of_v<TH1, H>>>
class LFInvMassFitter : public TNamed 
{
public:
    LFInvMassFitter() = default;
    LFInvMassFitter(const std::string& assoc, const std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay>& Histo2D, 
                    /*const std::array<std::vector<TH2*>, nbin_deltay>& HistoMultInt2D,*/
                    const std::string& OutPath, const std::string& OutFileName, int mode = 0);
    LFInvMassFitter(const std::string& assoc, const char* filename, const std::vector<std::string>& requiredKeys,
                    const std::vector<AxisCut>& slicing, int nbin_pT, const std::string& histoname);
    ~LFInvMassFitter();

    std::pair<Double_t, Double_t> GetPhiPurityAndError(TH1* h1PhiInvMass, std::string nameCanvas, Int_t isDataOrReco, Int_t isK0SOrPi, std::vector<Int_t> indices, bool printCanvas = false);
    
    std::pair<Double_t, Double_t> FitPhiAssoc(TH2* h2PhiAssocInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file/*,
                                                const std::vector<Double_t>& params, const std::vector<Double_t>& lowLimits, const std::vector<Double_t>& upLimits*/);

    Yield DoFit(std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, Double_t nSigma, TFile* file);
    Yield DoFit2(std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, Double_t nSigma);

    void ExportYields(Int_t nbin_pT, const std::vector<Double_t>& pT_axis, const std::string& hSetName, Int_t isTPCOrTOF);
    void ExportYields2(Int_t nbin_pT, const std::vector<Double_t>& pT_axis, const std::string& hSetName, Int_t isTPCOrTOF);

    //void CheckValidMembers();

    RooAbsPdf* CreateBackgroundFitFunction() const;
    RooAbsPdf* CreateSignalFitFunction();

    void DoFit();

private:

    ParticleType mParticleType = ParticleType::Unknown; // Default association particle type

    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> mSetHisto2D{};
    std::unordered_map<std::vector<int>, TH2*, VectorHash> mSetHisto{}; // Using unordered_map for faster access

    RooWorkspace* mWorkspace{nullptr}; // Pointer to the RooWorkspace for fitting

    std::string mOutputFileName;

    int mNEvents{0}; // Number of processed events
    int mMode{0}; // 0: Data, 1: Closure test

    //void HistogramAcquisition(const char* filename, int nbin_pT, const std::string& histoname);
    void HistogramAcquisition(const std::map<std::string, std::string>& meta, const std::vector<AxisCut>& slicing);

    std::pair<Double_t, Double_t> FitPhiK0S(TH2* h2PhiK0SInvMass, std::vector<Int_t> indices, TFile* file, Double_t nsigma = 6.0,
                                            const std::vector<Double_t>& params = {1., 1., 5., 5., 0.49, 0.003, -1.}, 
                                            const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., 0.48, 0.001, -1.7},
                                            const std::vector<Double_t>& upLimits = {2., 2., 10., 10., 0.5, 0.01, 10});
    std::pair<Double_t, Double_t> FitPhiPi(TH2* h2PhiPiInvMass, std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, TFile* file, Double_t nsigma = 4.0,
                                            const std::vector<Double_t>& params = {1., 1., 10., 10., 0., 1., 7., 1.}, 
                                            const std::vector<Double_t>& lowLimits = {1., 1., 1., 1., -1., 0.001, 3., 0.001}, 
                                            const std::vector<Double_t>& upLimits = {5., 5., 10., 10., 1.5, 2.5, 10., 5.});

    Yield FitPhiK0S2(std::vector<Int_t> indices, Double_t nSigma, TFile* file);
    Yield FitPhiPi2(std::vector<Int_t> indices, Int_t isTPCOrTOF, Int_t isDataOrMcReco, Double_t nSigma, TFile* file);

    void FillWorkspace(std::vector<Int_t> indices) const;

    ClassDef(LFInvMassFitter, 1);
};

