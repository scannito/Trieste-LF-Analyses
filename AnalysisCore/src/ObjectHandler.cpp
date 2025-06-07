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

#include <rapidjson/document.h>

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/Plot.h"

#include "AnalysisCore/include/ObjectHandler.h"

ObjectHandler::ObjectHandler() = default;
ObjectHandler::ObjectHandler(const char* filename, const std::vector<std::string>& requiredKeys) : fRequiredKeys(requiredKeys) { JSONParser(filename); }
ObjectHandler::~ObjectHandler()
{
    if (mTHnSparse)
        delete mTHnSparse;

}

TH2* ObjectHandler::GetHisto2D()
{
    if (!mTHnSparse) {
        std::cerr << "THnSparse is not initialized." << std::endl;
        return nullptr;
    }
    return Project2D(0, 0, 0, 1, 2);
}

TH2* ObjectHandler::GetHistoMultInt2D()
{
    if (!mTHnSparse) {
        std::cerr << "THnSparse is not initialized." << std::endl;
        return nullptr;
    }
    return Project2D(0, 0, 0, 1, 2);
}

std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> ObjectHandler::GetSetHisto2D(int nbin_pT, const std::string& hSetName)
{
    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> setHisto2D;
    for (int i = 0; i < nbin_deltay; ++i) {
        for (int j = 0; j < nbin_mult; ++j) {
            std::vector<TH2*> histos;
            for (int k = 0; k < nbin_pT; ++k) {
                std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                histos.push_back(Project2D(0, j + 1, j + 1, 1, k + 1, k + 1, 3, 2, "", hName));
            }
            setHisto2D[i][j] = histos;
        }
    }
    return setHisto2D;
}

std::array<std::vector<TH2*>, nbin_deltay> ObjectHandler::GetSetHistoMultInt2D(int nbin_pT, const std::string& hSetName)
{
    std::array<std::vector<TH2*>, nbin_deltay> setHistoMultInt2D;
    for (int i = 0; i < nbin_deltay; ++i) {
        std::vector<TH2*> histos;
        for (int k = 0; k < nbin_pT; ++k) {
            std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(k);
            histos.push_back(Project2D(0, 1, nbin_mult, 1, k+1, k+1, 3, 2, "", hName));
        }
        setHistoMultInt2D[i] = histos;
    }
    return setHistoMultInt2D;
}

void ObjectHandler::JSONParser(const char* filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return;
    }

    std::string jsonStr((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Parse JSON
    rapidjson::Document document;
    document.Parse(jsonStr.data());

    if (document.HasParseError()) {
        std::cerr << "Error parsing JSON" << std::endl;
        return;
    }

    // Convert JSON to std::map
    std::map<std::string, std::string> meta;
    for (auto itr = document.MemberBegin(); itr != document.MemberEnd(); ++itr) {
        meta[itr->name.GetString()] = itr->value.GetString();
    }

    bool allKeysPresent = true;
    for (const auto& key : fRequiredKeys) {
        if (meta.find(key) == meta.end()) {
            std::cerr << "Missing required key: " << key << std::endl;
            allKeysPresent = false;
        }
    }

    if (!allKeysPresent) {
        std::cerr << "Please provide all required keys in the JSON file" << std::endl;
        return;
    }

    TFile* file1 = TFile::Open(meta["inputFile"].data());
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error opening input file: " << meta["inputFile"] << std::endl;
        return;
    }
    TDirectoryFile* phik0shortanalysis = (TDirectoryFile*)file1->Get(meta["analysisDir"].data());
    if (!phik0shortanalysis) {
        std::cerr << "Error retrieving analysis directory: " << meta["analysisDir"] << std::endl;
        return;
    }

    TDirectoryFile* eventHist = (TDirectoryFile*)phik0shortanalysis->Get(meta["eventHistDir"].data());
    if (!eventHist) {
        std::cerr << "Error retrieving event histogram directory: " << meta["eventHistDir"] << std::endl;
        return;
    }

    TH1F* hEventSelection = (TH1F*)eventHist->Get(meta["eventHistName"].data());
    hEventSelection->SetDirectory(0);
    Int_t binNumber = hEventSelection->GetXaxis()->FindBin(meta["binEventHistName"].data());
    Double_t nEventsPhi = hEventSelection->GetBinContent(binNumber);

    TDirectoryFile* PhiAssocHist = (TDirectoryFile*)phik0shortanalysis->Get(meta["PhiAssocDir"].data());
    if (!PhiAssocHist) {
        std::cerr << "Error retrieving Phi Assoc histogram directory: " << meta["PhiAssocDir"] << std::endl;
        return;
    }

    THnSparseF* hnPhiAssoc = (THnSparseF*)PhiAssocHist->Get(meta["PhiAssocInvMassHistName"].data());
    if (!hnPhiAssoc) {
        std::cerr << "Error retrieving THnSparse: " << meta["PhiAssocInvMassHistName"] << std::endl;
        return;
    }

    std::string outPath = meta["outputPath"];
    std::string outFileName = meta["outputFile"];

    mTHnSparse = (THnSparseF*)hnPhiAssoc->Clone("hnPhiAssocClone");
    mNEvents = nEventsPhi;
    mOutPath = outPath;
    mOutFileName = outFileName;

    delete hEventSelection;
    delete hnPhiAssoc;
    delete PhiAssocHist;
    delete eventHist;
    delete phik0shortanalysis;
    file1->Close();
}

TH2* ObjectHandler::Project2D(Int_t axistocut, Int_t binlow, Int_t binup, Int_t axistoproj1, Int_t axistoproj2, Option_t* option, std::string hname) 
{ 
    mTHnSparse->GetAxis(axistocut)->SetRange(binlow, binup);
    TH2* h2 = (TH2*)mTHnSparse->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.data());
    h2->SetDirectory(0);
    return h2;
}

TH2* ObjectHandler::Project2D(Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistoproj1, Int_t axistoproj2, Option_t* option, std::string hname) 
{ 
    mTHnSparse->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    mTHnSparse->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    TH2* h2 = (TH2*)mTHnSparse->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.data());
    h2->SetDirectory(0);
    return h2;
}
TH1* ObjectHandler::Project1D(Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistocut3, Int_t binlow3, Int_t binup3, Int_t axistoproj, Option_t* option, std::string hname) 
{ 
    mTHnSparse->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    mTHnSparse->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    mTHnSparse->GetAxis(axistocut3)->SetRange(binlow3, binup3);
    TH1* h1 = (TH1*)mTHnSparse->Projection(axistoproj, option);
    h1->SetName(hname.data());
    h1->SetDirectory(0);
    return h1;
}
