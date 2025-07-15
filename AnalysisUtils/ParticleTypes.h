#pragma once

#include <string>
#include <map>

enum class ParticleType { Phi, K0S, Pion, PhiK0S, PhiPion, Unknown };

inline ParticleType StringToPart(const std::string& part) 
{
    std::map<std::string, ParticleType> map = {
        {"Phi", ParticleType::Phi},
        {"K0S", ParticleType::K0S},
        {"Pion", ParticleType::Pion},
        {"Phi-K0S", ParticleType::PhiK0S},
        {"Phi-Pi", ParticleType::PhiPion}
    };

    std::string key(part);

    auto it = map.find(key);
    return it != map.end() ? it->second : ParticleType::Unknown;
}

inline std::string PartToString(ParticleType part) 
{
    switch (part) {
        case ParticleType::Phi: return "Phi";
        case ParticleType::K0S: return "K0S";
        case ParticleType::Pion: return "Pion";
        case ParticleType::PhiK0S: return "Phi-K0S";
        case ParticleType::PhiPion: return "Phi-Pion";
        default: return "Unknown";
    }
}

inline std::string PartToSymbol(ParticleType part) 
{
    switch (part) {
        case ParticleType::Phi: return "#phi";
        case ParticleType::K0S: return "K^{0}_{S}";
        case ParticleType::Pion: return "#pi^{+}+#pi^{#minus}";
        case ParticleType::PhiK0S: return "#phi-K^{0}_{S}";
        case ParticleType::PhiPion: return "#phi-(#pi^{+}+#pi^{#minus})";
        default: return "Unknown";
    }
}