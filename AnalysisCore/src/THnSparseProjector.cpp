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
#include <memory>
#include <map>
#include <set>
#include <filesystem>

#include <rapidjson/document.h>

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/Plot.h"

#include "JSONReader.h"
#include "THnSparseProjector.h"

THnSparseProjector::THnSparseProjector(const char* filename, const std::vector<std::string>& requiredKeys)
{ 
    JSONReader jsonReader(filename, requiredKeys);
    THnSparseAcquisition(jsonReader.GetMeta());
    std::cout << "THnSparseProjector initialized with file: " << filename << std::endl; 
}

void THnSparseProjector::ExportProjections(int nbin_pT, const std::vector<AxisCut>& slicing, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
{
    static std::set<std::string> initializedFiles;

    // Delete the file if it exists, but only once per session
    if (initializedFiles.find(mOutputFileName) == initializedFiles.end()) {
        if (std::filesystem::exists(mOutputFileName)) {
            std::filesystem::remove(mOutputFileName);
            std::cout << "File " << mOutputFileName << " deleted (first use in this session)." << std::endl;
        }
        initializedFiles.insert(mOutputFileName);
    }

    std::unique_ptr<TFile> outputFile = std::unique_ptr<TFile>(TFile::Open(mOutputFileName.c_str(), "UPDATE"));
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error opening output file: " << mOutputFileName << std::endl;
        return;
    }

    auto cutSets = ExpandAxisCuts(slicing);

    for (const auto& combo : cutSets) {
        std::cout << "Writing histogram for bins:";
        std::string hName = hSetName;
        for (const auto& cut : combo) {
            std::cout << " axis " << cut.axis << " bin " << cut.binLow << ",";
            hName += "_" + std::to_string(cut.binLow);
        }
        std::cout << std::endl;

        TH2* h2 = Project2D(hName.c_str(), combo, axixtoproject);

        if (h2) {
            h2->Write();
            delete h2;
        } else {
            std::cerr << "Warning: failed to project histogram " << hName << "\n";
        }
    }

    /*for (int i = 0; i <= nbin_deltay; ++i) {
        for (int j = 0; j < nbin_mult; ++j) {
            for (int k = 0; k < nbin_pT; ++k) {
                std::cout << "Writing histogram for delta y bin: " << i << ", mult bin: " << j << ", pT bin: " << k << std::endl;
                std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                TH2* h2 = Project2D(hName.c_str(), {{0, i+1, i+1}, {1, j+1, j+1}, {2, k+1, k+1}}, axixtoproject);
                h2->Write();
                delete h2;
            }
        }
    }*/
}

void THnSparseProjector::THnSparseAcquisition(const std::map<std::string, std::string>& meta)
{
    std::unique_ptr<TFile> inputFile = std::unique_ptr<TFile>(TFile::Open(meta.at("inputFile").c_str()));
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << meta.at("inputFile") << std::endl;
        return;
    }

    THnSparse* tempTHnSparse = (THnSparse*)inputFile->Get(meta.at("objectPath").c_str());
    if (!tempTHnSparse) {
        std::cerr << "Error retrieving TObject: " << meta.at("objectPath") << std::endl;
        return;
    }
    THnSparse* cloneTHnSparse = (THnSparse*)tempTHnSparse->Clone("CloneTHnSparse");
    mTHnSparse.reset(cloneTHnSparse);

    mOutputFileName = meta.at("outputFile");
}

TH2* THnSparseProjector::Project2D(const char* hname, const std::vector<AxisCut>& axiscut, const std::pair<Int_t, Int_t>& axixtoproject, Option_t* option) 
{ 
    for (const auto& cut : axiscut)
        mTHnSparse->GetAxis(cut.axis)->SetRange(cut.binLow, cut.binHigh);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject.first, axixtoproject.second, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}

TH1* THnSparseProjector::Project1D(const char* hname, const std::vector<AxisCut>& axiscut, Int_t axixtoproject, Option_t* option) 
{ 
    for (const auto& cut : axiscut)
        mTHnSparse->GetAxis(cut.axis)->SetRange(cut.binLow, cut.binHigh);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}
