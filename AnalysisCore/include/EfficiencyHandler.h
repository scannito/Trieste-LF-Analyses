#pragma once 

#include <vector>
#include <string>

class EfficiencyHandler
{
public:
    EfficiencyHandler() = default;
    EfficiencyHandler(const char* filename, const std::vector<std::string>& requiredKeys);
    ~EfficiencyHandler() = default;
};
