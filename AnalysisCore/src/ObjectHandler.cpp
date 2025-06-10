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


ObjectHandler::ObjectHandler(const char* filename, const std::vector<std::string>& requiredKeys) : fRequiredKeys(requiredKeys) 
{ 
    std::map<std::string, std::string> meta;
    ReadFromJSON(meta, filename);
    ObjectAcquisition(meta);
    std::cout << "ObjectHandler initialized with file: " << filename << std::endl; 
}
ObjectHandler::~ObjectHandler()
{
    if (mTHnSparse)
        delete mTHnSparse;

}

/*std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> ObjectHandler::GetSetHisto2D(int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
{
    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> setHisto2D;
    for (int i = 1; i <= nbin_deltay; ++i) {
        for (int j = 0; j < nbin_mult; ++j) {
            std::vector<TH2*> histos;
            for (int k = 0; k < nbin_pT; ++k) {
                std::cout << "Creating histogram for delta y bin: " << i << ", mult bin: " << j << ", pT bin: " << k << std::endl;
                std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                histos.push_back(Project2D(hName.data(), {{0, i+1, i+1}, {1, j+1, j+1}, {2, k+1, k+1}}, axixtoproject));
            }
            setHisto2D[i-1][j] = histos;
        }
    }
    return setHisto2D;
}

std::array<std::vector<TH2*>, nbin_deltay> ObjectHandler::GetSetHistoMultInt2D(int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
{
    std::array<std::vector<TH2*>, nbin_deltay> setHistoMultInt2D;
    for (int i = 1; i <= nbin_deltay; ++i) {
        std::vector<TH2*> histos;
        for (int k = 0; k < nbin_pT; ++k) {
            std::cout << "Creating histogram for delta y bin: " << i << ", pT bin: " << k << std::endl;
            std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(k);
            histos.push_back(Project2D(hName.data(), {{0, i+1, i+1}, {1, 1, nbin_mult}, {2, k+1, k+1}}, axixtoproject));
        }
        setHistoMultInt2D[i-1] = histos;
    }
    return setHistoMultInt2D;
}*/

void ObjectHandler::ExportProjections(const char* filename, int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
{
    std::map<std::string, std::string> meta;
    TFile* file = TFile::Open(filename, "RECREATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening output file: " << filename << std::endl;
        return;
    }

    for (int i = 1; i <= nbin_deltay; ++i) {
        for (int j = 0; j < nbin_mult; ++j) {
            for (int k = 0; k < nbin_pT; ++k) {
                std::cout << "Writing histogram for delta y bin: " << i << ", mult bin: " << j << ", pT bin: " << k << std::endl;
                std::string hName = hSetName + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                TH2* h2 = Project2D(hName.data(), {{0, i+1, i+1}, {1, j+1, j+1}, {2, k+1, k+1}}, axixtoproject);
                h2->Write();
                delete h2;
            }
        }
    }

    file->Close();
}

void ObjectHandler::ReadFromJSON(std::map<std::string, std::string>& meta, const char* filename)
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
    }
}

void ObjectHandler::ObjectAcquisition(std::map<std::string, std::string>& meta)
{
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

TH2* ObjectHandler::Project2D(const char* hname, const std::vector<std::tuple<Int_t, Int_t, Int_t>>& axistocut, const std::pair<Int_t, Int_t>& axixtoproject, Option_t* option) 
{ 
    for (const auto [naxis, binlow, binup] : axistocut)
        mTHnSparse->GetAxis(naxis)->SetRange(binlow, binup);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject.first, axixtoproject.second, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}

TH1* ObjectHandler::Project1D(const char* hname, const std::vector<std::tuple<Int_t, Int_t, Int_t>>& axistocut, Int_t axixtoproject, Option_t* option) 
{ 
    for (const auto [naxis, binlow, binup] : axistocut)
        mTHnSparse->GetAxis(naxis)->SetRange(binlow, binup);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}

void ObjectHandler::CheckValidMembers()
{
    if (!mTHnSparse) {
        std::cerr << "Error: mTHnSparse is not initialized." << std::endl;
    } else {
        std::cout << "mTHnSparse is valid." << std::endl;
    
    if (mNEvents <= 0) {
        std::cerr << "Error: mNEvents is not set or invalid." << std::endl;
    } else {
        std::cout << "mNEvents is valid: " << mNEvents << std::endl;
    }

    if (mOutPath.empty()) {
        std::cerr << "Error: mOutPath is not set." << std::endl;
    } else {
        std::cout << "mOutPath is valid: " << mOutPath << std::endl;
    }

    if (mOutFileName.empty()) {
        std::cerr << "Error: mOutFileName is not set." << std::endl;
    } else {
        std::cout << "mOutFileName is valid: " << mOutFileName << std::endl;
    }
    }
}
