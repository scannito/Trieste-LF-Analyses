#pragma once

#include <string>
#include <map>

enum class ParticleType { Phi, K0S, Pi, PhiK0S, PhiPi, Unknown };

inline ParticleType StringToPart(const std::string& assoc) 
{
    std::map<std::string, ParticleType> map = {
        {"Phi", ParticleType::Phi},
        {"K0S", ParticleType::K0S},
        {"Pi", ParticleType::Pi},
        {"Phi-K0S", ParticleType::PhiK0S},
        {"Phi-Pi", ParticleType::PhiPi}
    };

    std::string key(assoc);

    auto it = map.find(key);
    return it != map.end() ? it->second : ParticleType::Unknown;
}

inline std::string PartToString(ParticleType assoc) 
{
    switch (assoc) {
        case ParticleType::Phi: return "Phi";
        case ParticleType::K0S: return "K0S";
        case ParticleType::Pi: return "Pi";
        case ParticleType::PhiK0S: return "Phi-K0S";
        case ParticleType::PhiPi: return "Phi-Pi";
        default: return "Unknown";
    }
}

inline std::string PartToSymbol(ParticleType assoc) 
{
    switch (assoc) {
        case ParticleType::Phi: return "#phi";
        case ParticleType::K0S: return "K^{0}_{S}";
        case ParticleType::Pi: return "#pi";
        case ParticleType::PhiK0S: return "#phi-K^{0}_{S}";
        case ParticleType::PhiPi: return "#phi-#pi";
        default: return "Unknown";
    }
}