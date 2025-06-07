#pragma once

#include "Riostream.h"
#include "TNamed.h"
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
#include <vector>
#include <string>

class ObjectHandler
{
public:
    ObjectHandler();
    ObjectHandler(const char* filename, const std::vector<std::string>& requiredKeys);
    ~ObjectHandler();

    TH2* GetHisto2D();
    TH2* GetHistoMultInt2D();

    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> GetSetHisto2D(int nbin_pT, const std::string& hSetName);
    std::array<std::vector<TH2*>, nbin_deltay> GetSetHistoMultInt2D(int nbin_pT, const std::string& hSetName);

    std::string GetOutPath() const { return mOutPath; }
    std::string GetOutFileName() const { return mOutFileName; }

private:
    std::vector<std::string> fRequiredKeys;

    THnSparseF* mTHnSparse;
    int mNEvents;
    std::string mOutPath;
    std::string mOutFileName;

    void JSONParser(const char* filename);

    TH2* Project2D(Int_t axistocut, Int_t binlow, Int_t binup, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", std::string hname = "");
    TH2* Project2D(Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", std::string hname = "");
    TH1* Project1D(Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistocut3, Int_t binlow3, Int_t binup3, Int_t axistoproj, Option_t* option = "", std::string hname = "");
};