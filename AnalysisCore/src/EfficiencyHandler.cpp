#include "Riostream.h"

#include "TFile.h"
#include "TH3.h"

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <set>
#include <filesystem>

#include "JSONReader.h"
#include "EfficiencyHandler.h"

#include "AnalysisUtils/Parameters.h"

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

    TH3* efficiencyHist = GetEfficiency();
    if (efficiencyHist) {
        efficiencyHist->Write();
        delete efficiencyHist;
    }

    TH3* signalLossHist = GetSignalLoss();
    if (signalLossHist) {
        signalLossHist->Write();
        delete signalLossHist;
    }

    TH3* combinedEffHist = GetEffXSigLoss();
    if (combinedEffHist) {
        combinedEffHist->Write();
        delete combinedEffHist;
    }

    if (mParticleType == ParticleType::Pion) {
        TH3* efficiencyHist2 = GetEfficiency2();
        if (efficiencyHist2) {
            efficiencyHist2->Write();
            delete efficiencyHist2;
        }

        TH3* combinedEffHist2 = GetCombEffXSigLoss();
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

    TH3* combinedEffHist = GetEffXSigLoss();
    if (combinedEffHist) {
        outputList->Add(combinedEffHist);
    }

    if (mParticleType == ParticleType::Pion) {
        TH3* combinedEffHist2 = GetCombEffXSigLoss();
        if (combinedEffHist2) {
            outputList->Add(combinedEffHist2);
        }
    }

    outputFile->cd();
    outputList->Write("ccdb_object", TObject::kOverwrite | TObject::kSingleKey);
    delete outputList;
}

TH3* EfficiencyHandler::GetEfficiency()
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

TH3* EfficiencyHandler::GetEfficiency2()
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

TH3* EfficiencyHandler::GetCombinedEfficiency()
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH3* combinedHist = GetEfficiency();
    combinedHist->Multiply(GetEfficiency2());
    combinedHist->SetName(Form("h3CombinedEfficiency%s", PartToString(mParticleType).c_str()));

    return combinedHist;
}

TH3* EfficiencyHandler::GetSignalLoss()
{
    TH3* signalLossHist = (TH3*)mGenAssocRecoHistogram->Clone(Form("SignalLossHistogram%s", PartToString(mParticleType).c_str()));
    signalLossHist->Divide(mGenAssocRecoHistogram.get(), mGenHistogram.get(), 1, 1, "B");

    return signalLossHist;
}

TH3* EfficiencyHandler::GetEffXSigLoss()
{
    TH3* combinedHist = GetEfficiency();
    combinedHist->Multiply(GetSignalLoss());
    if (mParticleType != ParticleType::Pion)
        combinedHist->SetName(Form("h3Efficiency%s", PartToString(mParticleType).c_str()));
    else
        combinedHist->SetName(Form("h3Efficiency%sTPC", PartToString(mParticleType).c_str()));

    return combinedHist;
}

TH3* EfficiencyHandler::GetCombEffXSigLoss()
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH3* combinedHist = GetCombinedEfficiency();
    combinedHist->Multiply(GetSignalLoss());
    combinedHist->SetName(Form("h3Efficiency%sTPCTOF", PartToString(mParticleType).c_str()));

    return combinedHist;
}

TH1* EfficiencyHandler::GetEfficiencySpectrum(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    TH1* projReco = mRecoHistogram->ProjectionY(Form("ProjEfficiency%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projReco) {
        std::cerr << "Error: Failed to project reco histogram." << std::endl;
        return nullptr;
    }
    TH1* projGenAssocReco = mGenAssocRecoHistogram->ProjectionY(Form("ProjSignalLoss%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projGenAssocReco) {
        std::cerr << "Error: Failed to project GenAssocReco histogram." << std::endl;
        return nullptr;
    }

    // Optional rebinning
    if (rebin && rebinnedpTAxis.size() > 1) {
        const int nBinsNew = rebinnedpTAxis.size() - 1;
        const double* binEdges = rebinnedpTAxis.data();

        TH1* rebinnedReco = projReco->Rebin(nBinsNew, Form("ProjEfficiencyRebinned%s", PartToString(mParticleType).c_str()), binEdges);
        TH1* rebinnedGenAssoc = projGenAssocReco->Rebin(nBinsNew, Form("ProjSignalLossRebinned%s", PartToString(mParticleType).c_str()), binEdges);

        delete projReco;
        delete projGenAssocReco;

        projReco = rebinnedReco;
        projGenAssocReco = rebinnedGenAssoc;
    }
    
    TH1* efficiencySpectrum = (TH1*)projReco->Clone(Form("EfficiencySpectrum%s", PartToString(mParticleType).c_str()));
    efficiencySpectrum->Divide(projReco, projGenAssocReco, 1, 1, "B");

    delete projReco; // Clean up the projected histograms
    delete projGenAssocReco; // Clean up the projected histograms

    return efficiencySpectrum;
}

TH1* EfficiencyHandler::GetEfficiencySpectrum2(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH1* projReco2 = mReco2Histogram->ProjectionY(Form("ProjEfficiency2%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projReco2) {
        std::cerr << "Error: Failed to project reco2 histogram." << std::endl;
        return nullptr;
    }
    TH1* projReco = mRecoHistogram->ProjectionY(Form("ProjEfficiency%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projReco) {
        std::cerr << "Error: Failed to project reco histogram." << std::endl;
        return nullptr;
    }

    // Optional rebinning
    if (rebin && rebinnedpTAxis.size() > 1) {
        const int nBinsNew = rebinnedpTAxis.size() - 1;
        const double* binEdges = rebinnedpTAxis.data();

        TH1* rebinnedReco2 = projReco2->Rebin(nBinsNew, Form("ProjEfficiency2Rebinned%s", PartToString(mParticleType).c_str()), binEdges);
        TH1* rebinnedReco = projReco->Rebin(nBinsNew, Form("ProjEfficiencyRebinned%s", PartToString(mParticleType).c_str()), binEdges);

        delete projReco2;
        delete projReco;

        projReco2 = rebinnedReco2;
        projReco = rebinnedReco;
    }

    TH1* efficiencySpectrum2 = (TH1*)projReco2->Clone(Form("EfficiencySpectrum2%s", PartToString(mParticleType).c_str()));
    efficiencySpectrum2->Divide(projReco2, projReco, 1, 1, "B");

    delete projReco2; // Clean up the projected histograms
    delete projReco; // Clean up the projected histograms

    return efficiencySpectrum2;
}

TH1* EfficiencyHandler::GetCombinedEfficiencySpectrum(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH1* efficiencySpectrum = GetEfficiencySpectrum(binMult, rebin, rebinnedpTAxis);
    TH1* efficiencySpectrum2 = GetEfficiencySpectrum2(binMult, rebin, rebinnedpTAxis);
    if (!efficiencySpectrum || !efficiencySpectrum2) {
        std::cerr << "Error: Failed to get efficiency spectrum histograms." << std::endl;
        return nullptr;
    }
    efficiencySpectrum->Multiply(efficiencySpectrum2);
    efficiencySpectrum->SetName(Form("h1CombinedEfficiency%s", PartToString(mParticleType).c_str()));

    delete efficiencySpectrum2; // Clean up the second efficiency spectrum histogram

    return efficiencySpectrum;
}

TH1* EfficiencyHandler::GetSignalLossSpectrum(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    TH1* projGenAssocReco = mGenAssocRecoHistogram->ProjectionY(Form("ProjSignalLoss%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projGenAssocReco) {
        std::cerr << "Error: Failed to project GenAssocReco histogram." << std::endl;
        return nullptr;
    }
    TH1* projGen = mGenHistogram->ProjectionY(Form("ProjGen%s", PartToString(mParticleType).c_str()), binMult + 1, binMult + 1, 1, nbin_deltay);
    if (!projGen) {
        std::cerr << "Error: Failed to project Gen histogram." << std::endl;
        return nullptr;
    }

    // Optional rebinning
    if (rebin && rebinnedpTAxis.size() > 1) {
        const int nBinsNew = rebinnedpTAxis.size() - 1;
        const double* binEdges = rebinnedpTAxis.data();

        TH1* rebinnedGenAssocReco = projGenAssocReco->Rebin(nBinsNew, Form("ProjSignalLossRebinned%s", PartToString(mParticleType).c_str()), binEdges);
        TH1* rebinnedGen = projGen->Rebin(nBinsNew, Form("ProjGenRebinned%s", PartToString(mParticleType).c_str()), binEdges);

        delete projGenAssocReco;
        delete projGen;

        projGenAssocReco = rebinnedGenAssocReco;
        projGen = rebinnedGen;
    }

    TH1* signalLossSpectrum = (TH1*)projGenAssocReco->Clone(Form("SignalLossSpectrum%s", PartToString(mParticleType).c_str()));
    signalLossSpectrum->Divide(projGenAssocReco, projGen, 1, 1, "B");

    delete projGenAssocReco; // Clean up the projected histograms
    delete projGen; // Clean up the projected histograms

    return signalLossSpectrum;
}

TH1* EfficiencyHandler::GetEffXSigLossSpectrum(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    TH1* efficiencySpectrum = GetEfficiencySpectrum(binMult, rebin, rebinnedpTAxis);
    TH1* signalLossSpectrum = GetSignalLossSpectrum(binMult, rebin, rebinnedpTAxis);
    if (!efficiencySpectrum || !signalLossSpectrum) {
        std::cerr << "Error: Failed to get efficiency or signal loss spectrum histogram." << std::endl;
        return nullptr;
    }
    efficiencySpectrum->Multiply(signalLossSpectrum);
    efficiencySpectrum->SetName(Form("h1Efficiency%s", PartToString(mParticleType).c_str()));

    delete signalLossSpectrum; // Clean up the signal loss spectrum histogram

    return efficiencySpectrum;
}

TH1* EfficiencyHandler::GetCombEffXSigLossSpectrum(int binMult, bool rebin, const std::vector<double>& rebinnedpTAxis)
{
    if (mParticleType != ParticleType::Pion) {
        std::cerr << "Error: This method is only applicable for Pion particle type." << std::endl;
        return nullptr;
    }
    TH1* efficiencySpectrum = GetCombinedEfficiencySpectrum(binMult, rebin, rebinnedpTAxis);
    TH1* signalLossSpectrum = GetSignalLossSpectrum(binMult, rebin, rebinnedpTAxis);
    if (!efficiencySpectrum || !signalLossSpectrum) {
        std::cerr << "Error: Failed to get efficiency or signal loss spectrum histogram." << std::endl;
        return nullptr;
    }

    efficiencySpectrum->Multiply(signalLossSpectrum);
    efficiencySpectrum->SetName(Form("h1Efficiency%sTPCTOF", PartToString(mParticleType).c_str()));

    delete signalLossSpectrum; // Clean up the signal loss spectrum histogram

    return efficiencySpectrum;
}