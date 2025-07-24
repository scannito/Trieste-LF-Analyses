#pragma once

#include "TH3.h"

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>

#include "AnalysisUtils/ParticleTypes.h"
 
class EfficiencyHandler
{
public:
    EfficiencyHandler() = default;
    EfficiencyHandler(const std::string& part, const char* filename, const std::vector<std::string>& requiredKeys);
    ~EfficiencyHandler() = default;

    void ExportCorrections();
    void ExportCorrectionsForCCDB();

    // To be generalized
    TH1* GetEfficiencySpectrum(int binMult);
    TH1* GetEfficiencySpectrum2(int binMult); // Only for Pi (at the moment, in future may be used for other charged particles)
    TH1* GetCombinedEfficiencySpectrum(int binMult); // Only for Pi (at the moment, in future may be used for other charged particles)
    TH1* GetSignalLossSpectrum(int binMult);
    TH1* GetEffXSigLossSpectrum(int binMult);
    TH1* GetCombEffXSigLossSpectrum(int binMult); // Only for Pi (at the moment, in future may be used for other charged particles)

private:
    ParticleType mParticleType = ParticleType::Unknown;

    std::unique_ptr<TH3> mRecoHistogram;
    std::unique_ptr<TH3> mReco2Histogram; // Only for Pi (at the moment, in future may be used for other charged particles)
    std::unique_ptr<TH3> mGenHistogram;
    std::unique_ptr<TH3> mGenAssocRecoHistogram;

    std::string mOutputFileName;

    void HistogramsAcquisition(const std::map<std::string, std::string>& meta);

    TH3* GetEfficiency();
    TH3* GetEfficiency2(); // Only for Pi (at the moment, in future may be used for other charged particles)
    TH3* GetCombinedEfficiency(); // Only for Pi (at the moment, in future may be used for other charged particles)
    TH3* GetSignalLoss();
    TH3* GetEffXSigLoss();
    TH3* GetCombEffXSigLoss(); // Only for Pi (at the moment, in future may be used for other charged particles)
};