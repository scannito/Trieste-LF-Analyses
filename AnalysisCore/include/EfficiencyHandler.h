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
    EfficiencyHandler(const std::string& assoc, const char* filename, const std::vector<std::string>& requiredKeys);
    ~EfficiencyHandler() = default;

    void ExportCorrections();
    void ExportCorrectionsForCCDB();

private:
    ParticleType mParticleType = ParticleType::Unknown;

    std::unique_ptr<TH3> mRecoHistogram;
    std::unique_ptr<TH3> mReco2Histogram; // Only for Pi(at the moment, in future may be used for other charged particles)
    std::unique_ptr<TH3> mGenHistogram;
    std::unique_ptr<TH3> mGenAssocRecoHistogram;

    std::string mOutputFileName;

    void HistogramsAcquisition(const std::map<std::string, std::string>& meta);

    TH3* GetEfficiencyHistogram();
    TH3* GetEfficiencyHistogram2(); // Only for Pi(at the moment, in future may be used for other charged particles)
    TH3* GetSignalLossHistogram();
    TH3* GetCombinedEfficiencyHistogram();
    TH3* GetCombinedEfficiencyHistogram2(); // Only for Pi(at the moment, in future may be used for other charged particles)
};