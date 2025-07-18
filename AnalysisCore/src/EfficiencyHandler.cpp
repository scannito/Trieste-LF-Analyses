#include "Riostream.h"

#include "TFile.h"

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <set>
#include <filesystem>

#include "JSONReader.h"
#include "EfficiencyHandler.h"

EfficiencyHandler::EfficiencyHandler(const std::string& part, const char* filename, const std::vector<std::string>& requiredKeys) : 
                                     mParticleType(StringToPart(part))
{
    JSONReader jsonReader(filename, requiredKeys);
    HistogramsAcquisition(jsonReader.GetMeta());
    std::cout << "EfficiencyHandler initialized with file: " << filename << std::endl;
    if (mParticleType == ParticleType::Unknown) {
        std::cerr << "Error: Unknown associated particle type. Please check the input." << std::endl;
        return;
    }
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
    delete tempRecoHistogram; // Clean up the original pointer

    if (mParticleType == ParticleType::Pion) {
        TH3* tempReco2Histogram = (TH3*)inputFile->Get(meta.at("recoPath2").c_str());
        if (!tempReco2Histogram) {
            std::cerr << "Error retrieving TObject: " << meta.at("recoPath2") << std::endl;
            return;
        }
        TH3* cloneReco2 = (TH3*)tempReco2Histogram->Clone("CloneReco2");
        cloneReco2->SetDirectory(0);
        mReco2Histogram.reset(cloneReco2);
        std::cout << "Directory of cloned reco2 histogram: " << mReco2Histogram->GetDirectory() << std::endl;
        delete tempReco2Histogram; // Clean up the original pointer
    }

    TH3* tempGenHistogram = (TH3*)inputFile->Get(meta.at("genPath").c_str());
    if (!tempGenHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genPath") << std::endl;
        return;
    }
    TH3* cloneGen = (TH3*)tempGenHistogram->Clone("CloneGen");
    cloneGen->SetDirectory(0);
    mGenHistogram.reset(cloneGen);
    std::cout << "Directory of cloned gen histogram: " << mGenHistogram->GetDirectory() << std::endl;
    delete tempGenHistogram; // Clean up the original pointer

    TH3* tempGenAssocRecoHistogram = (TH3*)inputFile->Get(meta.at("genAssocRecoPath").c_str());
    if (!tempGenAssocRecoHistogram) {
        std::cerr << "Error retrieving TObject: " << meta.at("genAssocRecoPath") << std::endl;
        return;
    }
    TH3* cloneGenAssocReco = (TH3*)tempGenAssocRecoHistogram->Clone("CloneGenAssocReco");
    cloneGenAssocReco->SetDirectory(0);
    mGenAssocRecoHistogram.reset(cloneGenAssocReco);
    std::cout << "Directory of cloned gen assoc reco histogram: " << mGenAssocRecoHistogram->GetDirectory() << std::endl;
    delete tempGenAssocRecoHistogram; // Clean up the original pointer

    mOutputFileName = meta.at("outputFile");
}

void EfficiencyHandler::ExportCorrections()
{
    std::unique_ptr<TFile> outputFile =  std::unique_ptr<TFile>(TFile::Open(mOutputFileName.c_str(), "UPDATE"));
    if (!outputFile || outputFile->IsZombie()) {
        outputFile = std::unique_ptr<TFile>(TFile::Open(mOutputFileName.c_str(), "RECREATE"));
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

    if (mParticleType == ParticleType::Pion) {
        TH3* efficiencyHist2 = GetEfficiencyHistogram2();
        if (efficiencyHist2) {
            efficiencyHist2->Write();
            delete efficiencyHist2;
        }

        TH3* combinedEffHist2 = GetCombinedEfficiencyHistogram2();
        if (combinedEffHist2) {
            combinedEffHist2->Write();
            delete combinedEffHist2;
        }
    }
}

void EfficiencyHandler::ExportCorrectionsForCCDB()
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

    TList* outputList = (TList*)outputFile->Get("ccdb_object");
    if (!outputList) {
        outputList = new TList();
        outputList->SetName("ccdb_object");
    }

    TH3* combinedEffHist = GetCombinedEfficiencyHistogram();
    if (combinedEffHist) {
        outputList->Add(combinedEffHist);
    }

    if (mParticleType == ParticleType::Pion) {
        TH3* combinedEffHist2 = GetCombinedEfficiencyHistogram2();
        if (combinedEffHist2) {
            outputList->Add(combinedEffHist2);
        }
    }

    outputFile->cd();
    outputList->Write("ccdb_object", TObject::kOverwrite | TObject::kSingleKey);
    delete outputList;
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
    TH3* efficiencyHist = (TH3*)mRecoHistogram->Clone(Form("EfficiencyHistogram%s", PartToString(mParticleType).c_str()));
    std::cout << "Efficiency histogram cloned." << std::endl;
    efficiencyHist->Divide(mRecoHistogram.get(), mGenAssocRecoHistogram.get(), 1, 1, "B");

    return efficiencyHist;

}

TH3* EfficiencyHandler::GetEfficiencyHistogram2()
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    if (!mReco2Histogram || !mRecoHistogram) {
        std::cerr << "Error: Reco2 or Reco histogram is not initialized." << std::endl;
        return nullptr;
    }
    TH3* efficiencyHist2 = (TH3*)mReco2Histogram->Clone(Form("EfficiencyHistogram2%s", PartToString(mParticleType).c_str()));
    efficiencyHist2->Divide(mReco2Histogram.get(), mRecoHistogram.get(), 1, 1, "B");

    return efficiencyHist2;
}

TH3* EfficiencyHandler::GetSignalLossHistogram()
{
    TH3* signalLossHist = (TH3*)mGenAssocRecoHistogram->Clone(Form("SignalLossHistogram%s", PartToString(mParticleType).c_str()));
    signalLossHist->Divide(mGenAssocRecoHistogram.get(), mGenHistogram.get(), 1, 1, "B");

    return signalLossHist;
}

TH3* EfficiencyHandler::GetCombinedEfficiencyHistogram()
{
    TH3* combinedHist = GetEfficiencyHistogram();
    combinedHist->Multiply(GetSignalLossHistogram());
    if (mParticleType != ParticleType::Pion)
        combinedHist->SetName(Form("h3Efficiency%s", PartToString(mParticleType).c_str()));
    else
        combinedHist->SetName(Form("h3Efficiency%sTPC", PartToString(mParticleType).c_str()));

    return combinedHist;
}

TH3* EfficiencyHandler::GetCombinedEfficiencyHistogram2()
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH3* combinedHist2 = GetEfficiencyHistogram2();
    combinedHist2->Multiply(GetSignalLossHistogram());
    combinedHist2->SetName(Form("h3Efficiency%sTPCTOF", PartToString(mParticleType).c_str()));

    return combinedHist2;
}