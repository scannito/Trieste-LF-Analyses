#include "Riostream.h"

#include "TFile.h"

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>

#include "JSONReader.h"
#include "EfficiencyHandler.h"

EfficiencyHandler::EfficiencyHandler(const char* filename, const std::vector<std::string>& requiredKeys)
{
    JSONReader jsonReader(filename, requiredKeys);
    HistogramsAcquisition(jsonReader.GetMeta());
    std::cout << "EfficiencyHandler initialized with file: " << filename << std::endl;
}

void EfficiencyHandler::HistogramsAcquisition(const std::map<std::string, std::string>& meta)
{
    std::unique_ptr<TFile> inputFile = std::unique_ptr<TFile>(TFile::Open(meta.at("inputFile").c_str()));
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << meta.at("inputFile") << std::endl;
        return;
    }

    TH3* tempRecoHistogram = (TH3*)inputFile->Get(meta.at("recoPath").c_str());
    if (!tempRecoHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("recoPath") << std::endl;
        return;
    }
    TH3* cloneReco = (TH3*)tempRecoHistogram->Clone("CloneReco");
    cloneReco->SetDirectory(0);
    mRecoHistogram.reset(cloneReco);
    std::cout << "Directory of cloned reco histogram: " << mRecoHistogram->GetDirectory() << std::endl;

    TH3* tempGenHistogram = (TH3*)inputFile->Get(meta.at("genPath").c_str());
    if (!tempGenHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genPath") << std::endl;
        return;
    }
    TH3* cloneGen = (TH3*)tempGenHistogram->Clone("CloneGen");
    cloneGen->SetDirectory(0);
    mGenHistogram.reset(cloneGen);
    std::cout << "Directory of cloned gen histogram: " << mGenHistogram->GetDirectory() << std::endl;

    TH3* tempGenAssocRecoHistogram = (TH3*)inputFile->Get(meta.at("genAssocRecoPath").c_str());
    if (!tempGenAssocRecoHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genAssocRecoPath") << std::endl;
        return;
    }
    TH3* cloneGenAssocReco = (TH3*)tempGenAssocRecoHistogram->Clone("CloneGenAssocReco");
    cloneGenAssocReco->SetDirectory(0);
    mGenAssocRecoHistogram.reset(cloneGenAssocReco);
    std::cout << "Directory of cloned gen assoc reco histogram: " << mGenAssocRecoHistogram->GetDirectory() << std::endl;

    mOutputFileName = meta.at("outputFile");
}

void EfficiencyHandler::ExportCorrections()
{
    TFile* outputFile = TFile::Open(mOutputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error opening output file: " << mOutputFileName << std::endl;
        return;
    }

    TH3* efficiencyHist = GetEfficiencyHistogram();
    if (efficiencyHist) {
        efficiencyHist->Write();
        delete efficiencyHist;
    }

    TH3* signalLossHist = GetSignalLossHistogram();
    if (signalLossHist) {
        signalLossHist->Write();
        delete signalLossHist;
    }

    TH3* combinedEffHist = GetCombinedEfficiencyHistogram();
    if (combinedEffHist) {
        combinedEffHist->Write();
        delete combinedEffHist;
    }

    outputFile->Close();
    delete outputFile;
}

void EfficiencyHandler::ExportCorrectionsForCCDB()
{
    TFile* outputFile = TFile::Open(mOutputFileName.c_str(), "UPDATE");
    if (!outputFile || outputFile->IsZombie()) {
        outputFile = TFile::Open(mOutputFileName.c_str(), "RECREATE");
    }

    TList* outputList = (TList*)outputFile->Get("ccdb_object");
    if (!outputList) {
        outputList = new TList();
        outputList->SetName("ccdb_object");
    }

    TH3* combinedEffHist = GetCombinedEfficiencyHistogram();
    if (combinedEffHist) {
        outputList->Add(combinedEffHist);
    }

    outputFile->cd();
    outputList->Write("ccdb_object", TObject::kOverwrite | TObject::kSingleKey);
    outputFile->Close();
    delete outputList;
    delete outputFile;
}

TH3* EfficiencyHandler::GetEfficiencyHistogram()
{
    std::cout << "Calculating efficiency histogram..." << std::endl;
    if (!mRecoHistogram || !mGenAssocRecoHistogram) {
        std::cerr << "Error: Reco or GenAssocReco histogram is not initialized." << std::endl;
        return nullptr;
    } else {
        std::cout << "Reco and GenAssocReco histograms are initialized." << std::endl;
        std::cout << "Cloning histogram: " << mRecoHistogram->GetName() << std::endl;
    }
    TH3* efficiencyHist = (TH3*)mRecoHistogram->Clone("EfficiencyHistogram");
    std::cout << "Efficiency histogram cloned." << std::endl;
    efficiencyHist->Divide(mRecoHistogram.get(), mGenAssocRecoHistogram.get(), 1, 1, "B");

    return efficiencyHist;

}

TH3* EfficiencyHandler::GetSignalLossHistogram()
{
    TH3* signalLossHist = (TH3*)mGenAssocRecoHistogram->Clone("SignalLossHistogram");
    signalLossHist->Divide(mGenAssocRecoHistogram.get(), mGenHistogram.get(), 1, 1, "B");

    return signalLossHist;
}

TH3* EfficiencyHandler::GetCombinedEfficiencyHistogram()
{
    TH3* combinedHist = GetEfficiencyHistogram();
    combinedHist->Multiply(GetSignalLossHistogram());

    return combinedHist;
}