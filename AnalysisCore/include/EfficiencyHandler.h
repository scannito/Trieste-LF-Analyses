#pragma once

#include "TH3.h"

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>

class EfficiencyHandler
{
public:
    EfficiencyHandler() = default;
    EfficiencyHandler(const char* filename, const std::vector<std::string>& requiredKeys);
    ~EfficiencyHandler() = default;

    void ExportCorrections();
    void ExportCorrectionsForCCDB();

private:
    std::unique_ptr<TH3> mRecoHistogram;
    std::unique_ptr<TH3> mGenHistogram;
    std::unique_ptr<TH3> mGenAssocRecoHistogram;

    std::string mOutputFileName;

    void HistogramsAcquisition(const std::map<std::string, std::string>& meta);

    TH3* GetEfficiencyHistogram();
    TH3* GetSignalLossHistogram();
    TH3* GetCombinedEfficiencyHistogram();
};