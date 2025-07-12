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
    TFile* inputFile = TFile::Open(meta.at("inputFile").data());
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << meta.at("inputFile") << std::endl;
        return;
    }

    TH3* tempRecoHistogram = (TH3*)inputFile->Get(meta.at("recoPath").data());
    if (!tempRecoHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("recoPath") << std::endl;
        return;
    }
    mRecoHistogram = std::unique_ptr<TH3>((TH3*)tempRecoHistogram->Clone("CloneReco"));

    TH3* tempGenHistogram = (TH3*)inputFile->Get(meta.at("genPath").data());
    if (!tempGenHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genPath") << std::endl;
        return;
    }
    mGenHistogram = std::unique_ptr<TH3>((TH3*)tempGenHistogram->Clone("CloneGen"));

    TH3* tempGenAssocRecoHistogram = (TH3*)inputFile->Get(meta.at("genAssocRecoPath").data());
    if (!tempGenAssocRecoHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genAssocRecoPath") << std::endl;
        return;
    }
    mGenAssocRecoHistogram = std::unique_ptr<TH3>((TH3*)tempGenAssocRecoHistogram->Clone("CloneGenAssocReco"));

    mOutputFileName = meta.at("outputFile");

    delete tempRecoHistogram;
    delete tempGenHistogram;
    delete tempGenAssocRecoHistogram;

    inputFile->Close();
    delete inputFile;
}

void EfficiencyHandler::ExportCorrection()
{
    TFile* outputFile = TFile::Open(mOutputFileName.data(), "RECREATE");
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

    TH3* combinedHist = GetCombinedEfficiencyHistogram();
    if (combinedHist) {
        combinedHist->Write();
        delete combinedHist;
    }

    outputFile->Close();
    delete outputFile;
}

TH3* EfficiencyHandler::GetEfficiencyHistogram()
{
    TH3* efficiencyHist = (TH3*)mRecoHistogram->Clone("EfficiencyHistogram");
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