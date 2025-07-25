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
#include <memory>
#include <map>

#include "AnalysisUtils/AxesUtils.h"

class THnSparseProjector 
{
public:
    THnSparseProjector() = default;
    THnSparseProjector(const char* filename, const std::vector<std::string>& requiredKeys);
    ~THnSparseProjector() = default;

    void ExportProjections(int nbin_pT, const std::vector<AxisCut>& slicing, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject);

private:
    std::unique_ptr<THnSparse> mTHnSparse;

    std::string mOutputFileName;

    void THnSparseAcquisition(const std::map<std::string, std::string>& meta);

    TH2* Project2D(const char* hname, const std::vector<AxisCut>& axiscut, const std::pair<Int_t, Int_t>& axixtoproject, Option_t* option = "");
    TH1* Project1D(const char* hname, const std::vector<AxisCut>& axiscut, Int_t axixtoproject, Option_t* option = "");
};