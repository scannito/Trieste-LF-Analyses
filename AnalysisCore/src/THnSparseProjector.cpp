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
#include <map>

#include <rapidjson/document.h>

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/Plot.h"

#include "JSONReader.h"
#include "THnSparseProjector.h"

THnSparseProjector::THnSparseProjector(const char* filename, const std::vector<std::string>& requiredKeys)
{ 
    JSONReader jsonReader(filename, requiredKeys);
    ObjectAcquisition(jsonReader.GetMeta());
    std::cout << "THnSparseProjector initialized with file: " << filename << std::endl; 
}

/*std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> THnSparseProjector::GetSetHisto2D(int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
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

std::array<std::vector<TH2*>, nbin_deltay> THnSparseProjector::GetSetHistoMultInt2D(int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
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

void THnSparseProjector::ExportProjections(int nbin_pT, const std::string& hSetName, const std::pair<Int_t, Int_t>& axixtoproject)
{
    std::map<std::string, std::string> meta;
    TFile* outputFile = TFile::Open(mOutputFileName.data(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error opening output file: " << mOutputFileName << std::endl;
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

    outputFile->Close();
    delete outputFile;
}

void THnSparseProjector::ObjectAcquisition(const std::map<std::string, std::string>& meta)
{
    TFile* inputFile = TFile::Open(meta.at("inputFile").data());
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << meta.at("inputFile") << std::endl;
        return;
    }

    THnSparseF* hnPhiAssoc = (THnSparseF*)inputFile->Get(meta.at("objectPath").data());
    if (!hnPhiAssoc) {
        std::cerr << "Error retrieving TObject: " << meta.at("objectPath") << std::endl;
        return;
    }

    mTHnSparse = std::unique_ptr<THnSparse>((THnSparseF*)hnPhiAssoc->Clone("hnPhiAssocClone"));

    /*TH1F* hEventSelection = (TH1F*)eventHist->Get(meta.at("eventHistName").data());
    hEventSelection->SetDirectory(0);
    Int_t binNumber = hEventSelection->GetXaxis()->FindBin(meta.at("binEventHistName").data());
    Double_t nEventsPhi = hEventSelection->GetBinContent(binNumber);*/

    mOutputFileName = meta.at("outputFile");

    delete hnPhiAssoc;

    inputFile->Close();
    delete inputFile;
}

TH2* THnSparseProjector::Project2D(const char* hname, const std::vector<AxisToCut>& axistocut, const std::pair<Int_t, Int_t>& axixtoproject, Option_t* option) 
{ 
    for (const auto [naxis, binlow, binup] : axistocut)
        mTHnSparse->GetAxis(naxis)->SetRange(binlow, binup);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject.first, axixtoproject.second, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}

TH1* THnSparseProjector::Project1D(const char* hname, const std::vector<AxisToCut>& axistocut, Int_t axixtoproject, Option_t* option) 
{ 
    for (const auto [naxis, binlow, binup] : axistocut)
        mTHnSparse->GetAxis(naxis)->SetRange(binlow, binup);

    TH2* h2 = (TH2*)mTHnSparse->Projection(axixtoproject, option);
    h2->SetName(hname);
    h2->SetDirectory(0);
    return h2;
}

void THnSparseProjector::CheckValidMembers()
{
    if (!mTHnSparse)
        std::cerr << "Error: mTHnSparse is not initialized." << std::endl;
    else
        std::cout << "mTHnSparse is valid." << std::endl;

    if (mOutputFileName.empty())
        std::cerr << "Error: mOutFileName is not set." << std::endl;
    else
        std::cout << "mOutFileName is valid: " << mOutputFileName << std::endl;
}
